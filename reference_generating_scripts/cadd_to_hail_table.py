#!/usr/bin/env python3

"""
Build the CADD v1.7 GRCh38 Hail Table from the raw score TSVs.

Inputs are the CADD_v1.7_snvs and CADD_v1.7_indels files already in the references
bucket (~8.8 billion SNV rows plus the gnomAD-observed indels). Output is one table keyed
by (locus, alleles) with raw_score and phred as float32, partitioned on the
hail_intervals_hg38 variant-balanced intervals so downstream joins against tables read on
the same intervals are co-partitioned merges. The read itself splits by bgzip block; the
intervals decide the layout that is written.

Paths default to the deployed references config (cpg_utils reference_path), so the
CADD_v1.7_ht key must be deployed before the default --out resolves. Runs on Hail
Query-on-Batch, writing straight to cpg-common-main, e.g.

    analysis-runner --dataset agdd --access-level full \
        --output-dir references/cadd \
        --description "CADD v1.7 Hail Table" \
        python3 reference_generating_scripts/cadd_to_hail_table.py

An existing --out is refused unless --overwrite is given. Each file is sorted by
position, so keyed on its own it should pass Hail's sortedness check without sort rounds;
if it does not, the SNV sort is a few hundred core-hours.
"""

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


def import_cadd(path: str) -> hl.Table:
    """
    Read one CADD score TSV (#Chrom Pos Ref Alt RawScore PHRED) as a keyed table.

    Args:
        path: bgzipped CADD TSV
    """
    ht = hl.import_table(
        path,
        force_bgz=True,
        comment='#',
        no_header=True,
        types={
            'f0': hl.tstr,
            'f1': hl.tint32,
            'f2': hl.tstr,
            'f3': hl.tstr,
            'f4': hl.tfloat32,
            'f5': hl.tfloat32,
        },
    )
    ht = ht.rename(
        {
            'f0': 'chrom',
            'f1': 'pos',
            'f2': 'ref',
            'f3': 'alt',
            'f4': 'raw_score',
            'f5': 'phred',
        }
    )
    contig = hl.literal(CONTIG_RECODING).get(ht.chrom, ht.chrom)
    ht = ht.transmute(
        locus=hl.locus(contig, ht.pos, reference_genome='GRCh38'),
        alleles=[ht.ref, ht.alt],
    )
    return ht.key_by('locus', 'alleles')


def main(snvs: str, indels: str, intervals_bed: str, out: str, overwrite: bool):
    refuse_existing(out, overwrite)
    init_batch(driver_cores=2, driver_memory='highmem')
    intervals = read_intervals(intervals_bed)

    # Imported one file at a time so each key_by sees a position-sorted input. One import
    # over both files would lay the indel partitions after the SNV ones, covering the
    # genome twice, and force a distributed sort of every row. The union of two keyed
    # tables is a range merge, not a shuffle; the SNV and indel keys never collide.
    ht = import_cadd(snvs).union(import_cadd(indels))

    write_on_intervals(
        ht,
        intervals,
        out,
        overwrite,
        source=dict(
            version='CADD v1.7 GRCh38', snvs=snvs, indels=indels, intervals=intervals_bed
        ),
    )


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--snvs', help="default reference_path('CADD_v1.7_snvs')")
    parser.add_argument('--indels', help="default reference_path('CADD_v1.7_indels')")
    parser.add_argument(
        '--intervals-bed',
        help=(
            "default reference_path('hail_intervals_hg38/"
            "gnomad_v4_1_variants_balanced_intervals_bed')"
        ),
    )
    parser.add_argument('--out', help="default reference_path('CADD_v1.7_ht')")
    parser.add_argument('--overwrite', action='store_true', help='replace an existing --out')
    args = parser.parse_args()
    main(
        args.snvs or reference_path('CADD_v1.7_snvs'),
        args.indels or reference_path('CADD_v1.7_indels'),
        args.intervals_bed
        or reference_path('hail_intervals_hg38/gnomad_v4_1_variants_balanced_intervals_bed'),
        args.out or reference_path('CADD_v1.7_ht'),
        args.overwrite,
    )


if __name__ == '__main__':
    cli_main()
