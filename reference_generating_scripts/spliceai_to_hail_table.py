#!/usr/bin/env python3

"""
Build the SpliceAI v1.3 GRCh38 Hail Table from the masked score VCFs.

Inputs are the spliceai_resources masked SNV and indel VCFs already in the references
bucket. SpliceAI scores a variant once per overlapping gene, one VCF record each, with
INFO/SpliceAI = ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL. Output is
one row per (locus, alleles) holding every gene's scores in a dict keyed by gene symbol,
plus ds_max across genes and score types. Only chr1-22, X, Y, M are kept: the files carry
records on 18 alt and random contigs that their headers do not declare, and those are
dropped at import. Partitioned on the hail_intervals_hg38 variant-balanced intervals.

Paths default to the deployed references config (cpg_utils reference_path), so the
spliceai_v1-3_ht key must be deployed before the default --out resolves. Runs on Hail
Query-on-Batch, writing straight to cpg-common-main, e.g.

    analysis-runner --dataset agdd --access-level full \
        --output-dir references/spliceai \
        --description "SpliceAI v1.3 Hail Table" \
        python3 reference_generating_scripts/spliceai_to_hail_table.py

An existing --out is refused unless --overwrite is given. Each VCF is sorted, so imported
one at a time there is no shuffle: import, collect records per key, union the two keyed
tables, write.
"""

import re
from argparse import ArgumentParser

import hail as hl
from cpg_utils.config import reference_path
from cpg_utils.hail_batch import init_batch
from hail_reference_utils import (
    CONTIG_RECODING,
    read_intervals,
    refuse_existing,
    write_on_intervals,
)

SCORE_FIELDS = ['ds_ag', 'ds_al', 'ds_dg', 'ds_dl']
POSITION_FIELDS = ['dp_ag', 'dp_al', 'dp_dg', 'dp_dl']
# Which variants each masked file may contain, by the name used on the command line.
VARIANT_KIND = {'snv': hl.is_snp, 'indel': hl.is_indel}
# import_vcf drops any line this matches: a record whose contig is not chr1-22, X, Y, M,
# with or without the chr prefix. The remaining records stay in GRCh38 contig order, so
# the import needs no sort. Header lines are left alone.
DROP_OTHER_CONTIGS = (
    '^(?!#|('
    + '|'.join(re.escape(c) for c in {*CONTIG_RECODING, *CONTIG_RECODING.values()})
    + r')\t)'
)


def parse_scores(parts: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """One split INFO/SpliceAI entry into (gene symbol, typed scores)."""
    return hl.struct(
        symbol=parts[1],
        scores=hl.struct(
            **{name: hl.float32(parts[i + 2]) for i, name in enumerate(SCORE_FIELDS)},
            **{name: hl.int32(parts[i + 6]) for i, name in enumerate(POSITION_FIELDS)},
        ),
    )


def collect_by_variant(vcf_path: str, kind: str) -> hl.Table:
    """
    One masked SpliceAI VCF as a keyed table with one row per (locus, alleles).

    Args:
        vcf_path: the masked VCF
        kind: 'snv' or 'indel'; a record of the other kind fails the run, since the
            union of the two tables relies on their keys never colliding.
    """
    ht = hl.import_vcf(
        vcf_path,
        reference_genome='GRCh38',
        contig_recoding=CONTIG_RECODING,
        force_bgz=True,
        filter=DROP_OTHER_CONTIGS,
    ).rows()
    # The field is Number=., so it arrives as an array: one element per gene in the record
    # (SpliceAI's own files write one record per gene, but nothing here relies on that).
    # explode gives one row per element and keeps the key order; a record without the field
    # would vanish in the explode, so it fails the run instead (a missing array makes the
    # case condition missing, which hl.case treats as false-and-missing, hence the coalesce).
    entries = (
        hl.case()
        .when(VARIANT_KIND[kind](ht.alleles[0], ht.alleles[1]), ht.info.SpliceAI)
        .or_error(f'not a {kind} at ' + hl.str(ht.locus))
    )
    ht = ht.select(
        entries=hl.case()
        .when(hl.coalesce(hl.len(entries), 0) > 0, entries)
        .or_error('record without a SpliceAI INFO value at ' + hl.str(ht.locus)),
    )
    ht = ht.explode('entries')
    # The ALLELE field must be the record's ALT, or the scores would be filed under the
    # wrong key with no error.
    parts = ht.entries.split('\\|')
    ht = ht.select(
        gene=hl.case()
        .when(parts[0] == ht.alleles[1], parse_scores(parts))
        .or_error('SpliceAI ALLELE differs from ALT at ' + hl.str(ht.locus)),
    )
    # Consecutive records share a key, so this groups without a shuffle.
    ht = ht.collect_by_key()
    # hl.dict keeps one entry per key, so a repeated symbol at one variant would drop a
    # record from by_gene while ds_max still counted it. Fail instead.
    by_gene = hl.dict(ht.values.map(lambda v: (v.gene.symbol, v.gene.scores)))
    return ht.select(
        by_gene=hl.case()
        .when(hl.len(by_gene) == hl.len(ht.values), by_gene)
        .or_error('repeated gene symbol at ' + hl.str(ht.locus)),
        ds_max=hl.max(
            ht.values.flatmap(
                lambda v: hl.array([v.gene.scores[f] for f in SCORE_FIELDS])
            )
        ),
    )


def main(snvs: str, indels: str, intervals_bed: str, out: str, overwrite: bool):
    refuse_existing(out, overwrite)
    # Same worker size as cadd_to_hail_table: a sort step holding a whole partition in
    # memory killed the default 1-core 3.75 GB worker JVM there. 4 highmem cores is 26 GB.
    init_batch(
        driver_cores=2,
        driver_memory='highmem',
        worker_cores=4,
        worker_memory='highmem',
    )
    intervals = read_intervals(intervals_bed)

    # Each file is sorted on its own, so imported on its own Hail keeps that order and
    # collect_by_key stays partition-local. One import over both files would interleave two
    # whole-genome partition sets and send ~14B rows through a distributed sort. The union of
    # two keyed tables is a range merge, not a shuffle; collect_by_variant checks each file
    # holds only its own kind of variant, so the keys never collide.
    ht = collect_by_variant(snvs, 'snv').union(collect_by_variant(indels, 'indel'))

    write_on_intervals(
        ht,
        intervals,
        out,
        overwrite,
        source=dict(
            version='SpliceAI v1.3 GRCh38 masked',
            snvs=snvs,
            indels=indels,
            intervals=intervals_bed,
        ),
    )


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        '--snvs', help="default reference_path('spliceai_resources/splice_ai_snvs')"
    )
    parser.add_argument(
        '--indels', help="default reference_path('spliceai_resources/splice_ai_indels')"
    )
    parser.add_argument(
        '--intervals-bed',
        help=(
            "default reference_path('hail_intervals_hg38/"
            "gnomad_v4_1_variants_balanced_intervals_bed')"
        ),
    )
    parser.add_argument('--out', help="default reference_path('spliceai_v1-3_ht')")
    parser.add_argument(
        '--overwrite', action='store_true', help='replace an existing --out'
    )
    args = parser.parse_args()
    main(
        args.snvs or reference_path('spliceai_resources/splice_ai_snvs'),
        args.indels or reference_path('spliceai_resources/splice_ai_indels'),
        args.intervals_bed
        or reference_path(
            'hail_intervals_hg38/gnomad_v4_1_variants_balanced_intervals_bed'
        ),
        args.out or reference_path('spliceai_v1-3_ht'),
        args.overwrite,
    )


if __name__ == '__main__':
    cli_main()
