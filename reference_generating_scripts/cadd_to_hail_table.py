#!/usr/bin/env python3

"""
Build the CADD v1.7 GRCh38 Hail Table from the raw score TSVs.

Inputs are the CADD_v1.7_snvs and CADD_v1.7_indels files already in the references
bucket (~8.8 billion SNV rows plus the gnomAD-observed indels). Output is one table keyed
by (locus, alleles) with raw_score and phred as float32, partitioned on the
hail_intervals_hg38 variant-balanced intervals so downstream joins against tables read on
the same intervals are co-partitioned merges. The read itself splits by bgzip block; the
intervals decide the layout that is written.

Runs on Hail Query-on-Batch, writing straight to cpg-common-main, e.g.

    analysis-runner --dataset agdd --access-level full \
        --output-dir references/cadd \
        --description "CADD v1.7 Hail Table" \
        python3 reference_generating_scripts/cadd_to_hail_table.py

The key_by is a full sort shuffle of the SNV file; expect a few hundred core-hours.
"""

import gzip
from argparse import ArgumentParser

import hail as hl
from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch

REFERENCES = 'gs://cpg-common-main/references'
CADD = f'{REFERENCES}/CADD/v1.7/GRCh38'
DEFAULT_SNVS = f'{CADD}/whole_genome_SNVs.tsv.gz'
DEFAULT_INDELS = f'{CADD}/gnomad.genomes.r4.0.indel.tsv.gz'
DEFAULT_INTERVALS = (
    f'{REFERENCES}/hail_intervals/hg38/gnomad_v4.1_variants_balanced_intervals.bed.gz'
)
DEFAULT_OUT = f'{CADD}/cadd_v1.7.ht'

# CADD names contigs without the chr prefix; MT is chrM in GRCh38.
CONTIG_RECODING = {
    **{str(c): f'chr{c}' for c in [*range(1, 23), 'X', 'Y']},
    'MT': 'chrM',
    'M': 'chrM',
}


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


def read_intervals(bed_path: str) -> list[hl.Interval]:
    """
    The variant-balanced intervals as Python Interval objects, for a partitioned read.

    Parsed here rather than with hl.import_bed, which cannot open a .gz path; the same
    parser ourdna_genomic_atlas uses for this file. BED is 0-based half-open, Hail loci
    are 1-based, so start + 1 with both ends included.
    """
    with to_path(bed_path).open('rb') as raw, gzip.open(raw, 'rt') as bed:
        rows = (line.split('\t') for line in bed if line.strip() and line[0] not in '#t')
        return [
            hl.Interval(
                hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38'),
                hl.Locus(chrom, int(end), reference_genome='GRCh38'),
                includes_start=True,
                includes_end=True,
            )
            for chrom, start, end, *_ in rows
        ]


def main(snvs: str, indels: str, intervals_bed: str, out: str):
    init_batch(driver_cores=2, driver_memory='highmem')
    intervals = read_intervals(intervals_bed)

    ht = import_cadd(snvs).union(import_cadd(indels))

    # The shuffle decides its own layout; re-read on the shared intervals before writing.
    tmp = hl.utils.new_temp_file('cadd_keyed', 'ht')
    ht.checkpoint(tmp)
    # _intervals is a private Hail argument, the standard idiom for read-time partitioning.
    ht = hl.read_table(tmp, _intervals=intervals)
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
    args = parser.parse_args()
    main(args.snvs, args.indels, args.intervals_bed, args.out)


if __name__ == '__main__':
    cli_main()
