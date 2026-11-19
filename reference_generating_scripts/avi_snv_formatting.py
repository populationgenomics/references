#! /usr/bin/env python3

"""
Expects the AlphaGenome Variant Impact (AVI) SNV table, available from:
- https://deepmind.google.com/science/alphagenome/downloads

Parses the block-gzipped (tabix-ready) TSV of AVI scores directly into a Hail Table,
keyed on locus & alleles, for annotating variant MTs.

Expected columns (header line starts with '#'):
#CHROM  POS  REF  ALT  raw_score  PHRED
"""

from argparse import ArgumentParser

import hail as hl

# positional column order in the source TSV -> names in the output Table
COLUMNS = ['chrom', 'pos', 'ref', 'alt', 'raw_score', 'phred']


def avi_tsv_to_hail_table(tsv_path: str, ht_out: str, min_partitions: int = 200):
    """
    Import the AVI TSV as a Hail Table, and write it out

    Args:
        tsv_path (str): path to the bgzipped AVI tsv.gz
        ht_out (str): path to write the Hail Table to
        min_partitions (int): minimum number of partitions to split the import into
    """

    # force_bgz: the file is BGZF (tabix-compatible) despite the .gz suffix, so Hail can split it and read in parallel
    # comment='#' drops the header (and any preamble), so columns are assigned by position
    ht = hl.import_table(
        tsv_path,
        no_header=True,
        comment='#',
        force_bgz=True,
        delimiter='\t',
        missing='',
        min_partitions=min_partitions,
        types={'f1': hl.tint32, 'f4': hl.tfloat32, 'f5': hl.tfloat32},
    )
    ht = ht.rename({f'f{index}': name for index, name in enumerate(COLUMNS)})

    ht = ht.transmute(
        locus=hl.locus(ht.chrom, ht.pos, reference_genome='GRCh38'),
        alleles=[ht.ref, ht.alt],
    )
    ht = ht.key_by('locus', 'alleles')
    ht.write(ht_out, overwrite=True)
    ht.describe()


def main(avi_tsv: str, ht_path: str):
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)
    avi_tsv_to_hail_table(avi_tsv, ht_out=ht_path)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--avi_tsv', help='path to the AVI SNV tsv.gz file', required=True)
    parser.add_argument('--ht_out', help='path to write a new Hail Table', required=True)
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(unknown)
    main(avi_tsv=args.avi_tsv, ht_path=args.ht_out)


if __name__ == '__main__':
    cli_main()
