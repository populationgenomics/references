#!/usr/bin/env python3

"""
Build the SpliceAI v1.3 GRCh38 Hail Table from the masked score VCFs.

Inputs are the spliceai_resources masked SNV and indel VCFs already in the references
bucket. SpliceAI scores a variant once per overlapping gene, one VCF record each, with
INFO/SpliceAI = ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL. Output is
one row per (locus, alleles) holding every gene's scores in a dict keyed by gene symbol,
plus ds_max across genes and score types. Every scored variant is kept. Partitioned on
the hail_intervals_hg38 variant-balanced intervals.

Runs on Hail Query-on-Batch, writing straight to cpg-common-main, e.g.

    analysis-runner --dataset agdd --access-level full \
        --output-dir references/spliceai \
        --description "SpliceAI v1.3 Hail Table" \
        python3 reference_generating_scripts/spliceai_to_hail_table.py

The VCFs are sorted, so no shuffle: import, collect records per key, write.
"""

from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch

REFERENCES = 'gs://cpg-common-main/references'
SPLICEAI = f'{REFERENCES}/ourdna_browser/v0/spliceai-resources/v1-3'
DEFAULT_SNVS = f'{SPLICEAI}/spliceai_scores.masked.snv.hg38.vcf.gz'
DEFAULT_INDELS = f'{SPLICEAI}/spliceai_scores.masked.indel.hg38.vcf.gz'
DEFAULT_INTERVALS = (
    f'{REFERENCES}/hail_intervals/hg38/gnomad_v4.1_variants_balanced_intervals.bed.gz'
)
DEFAULT_OUT = f'{SPLICEAI}/spliceai_v1-3.ht'

# The hg38 files name contigs without the chr prefix.
CONTIG_RECODING = {
    **{str(c): f'chr{c}' for c in [*range(1, 23), 'X', 'Y']},
    'MT': 'chrM',
    'M': 'chrM',
}
SCORE_FIELDS = ['ds_ag', 'ds_al', 'ds_dg', 'ds_dl']
POSITION_FIELDS = ['dp_ag', 'dp_al', 'dp_dg', 'dp_dl']


def parse_scores(record: hl.expr.StringExpression) -> hl.expr.StructExpression:
    """One INFO/SpliceAI string into (gene symbol, typed scores)."""
    parts = record.split('\\|')
    return hl.struct(
        symbol=parts[1],
        scores=hl.struct(
            **{name: hl.float32(parts[i + 2]) for i, name in enumerate(SCORE_FIELDS)},
            **{name: hl.int32(parts[i + 6]) for i, name in enumerate(POSITION_FIELDS)},
        ),
    )


def read_intervals(bed_path: str) -> list[hl.Interval]:
    """The variant-balanced intervals as a Python list, for a partitioned read."""
    return hl.import_bed(bed_path, reference_genome='GRCh38').interval.collect()


def main(snvs: str, indels: str, intervals_bed: str, out: str, min_partitions: int):
    init_batch(driver_cores=2, driver_memory='highmem')

    ht = hl.import_vcf(
        [snvs, indels],
        reference_genome='GRCh38',
        contig_recoding=CONTIG_RECODING,
        force_bgz=True,
        min_partitions=min_partitions,
        skip_invalid_loci=True,
    ).rows()
    # One record per gene; the field is Number=. so it arrives as a one-element array.
    ht = ht.select(gene=parse_scores(ht.info.SpliceAI[0]))
    # Consecutive records share a key, so this groups without a shuffle.
    ht = ht.collect_by_key()
    ht = ht.select(
        by_gene=hl.dict(ht.values.map(lambda v: (v.gene.symbol, v.gene.scores))),
        ds_max=hl.max(
            ht.values.flatmap(
                lambda v: hl.array([v.gene.scores[f] for f in SCORE_FIELDS])
            )
        ),
    )

    tmp = hl.utils.new_temp_file('spliceai_by_key', 'ht')
    ht.checkpoint(tmp)
    ht = hl.read_table(tmp, _intervals=read_intervals(intervals_bed))
    ht.write(out, overwrite=True)
    ht = hl.read_table(out)
    ht.describe()
    print(f'{ht.count():,} rows in {ht.n_partitions()} partitions at {out}')


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--snvs', default=DEFAULT_SNVS)
    parser.add_argument('--indels', default=DEFAULT_INDELS)
    parser.add_argument('--intervals-bed', default=DEFAULT_INTERVALS)
    parser.add_argument('--out', default=DEFAULT_OUT)
    parser.add_argument('--min-partitions', type=int, default=10000)
    args = parser.parse_args()
    main(args.snvs, args.indels, args.intervals_bed, args.out, args.min_partitions)


if __name__ == '__main__':
    cli_main()
