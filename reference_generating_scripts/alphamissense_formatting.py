#! /usr/bin/env python3

"""
Expects the from-source AlphaMissense table, available from:
- https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_hg38.tsv.gz
or
- https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz

Script for parsing the compressed tsv file of AlphaMissense results into a Hail Table
This Hail Table can be used for integrating AM scores into VCFs annotated with VEP or BCFtools

Process:
1. read through the compressed data, skip non-pathogenic entries
2. write the pathogenic entries back out to a reduced schema
3. parse that data as a Hail Table using specific variant types
4. write the Hail Table

Currently this skips over everything AM deems non-pathogenic, and only writes pathogenic variants
We may want to adjust that behaviour in the future.
"""


import gzip
import json
import requests

from cloudpathlib import AnyPath

from cpg_utils.hail_batch import init_batch

import hail as hl


AM_ZENODO = 'https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz'
DESTINATION = 'gs://cpg-common-test-tmp/references/alphamissense/AlphaMissense_hg38.tsv.gz'
TABLE_DESTINATION = 'gs://cpg-common-test-tmp/references/alphamissense/AlphaMissense_hg38.ht'


def process_header(final_header_line: str) -> dict[str, int]:
    """
    the TSV format is a little different in the all-isoforms vs. main transcript only
    this method determines the column indexes to use based on the header content
    we could determine preset columns by filename/flag, but... that's more burden on operator

    Args:
        final_header_line (str): the header line from the TSV file
    """
    # remove newline and hash, lowercase, split into a list
    broken_line = final_header_line.rstrip().replace('#', '').lower().split()

    return {
        'chrom': broken_line.index('chrom'),
        'pos': broken_line.index('pos'),
        'ref': broken_line.index('ref'),
        'alt': broken_line.index('alt'),
        'transcript': broken_line.index('transcript_id'),
        'am_pathogenicity': broken_line.index('am_pathogenicity'),
        'am_class': broken_line.index('am_class'),
    }


def filter_for_pathogenic_am(input_file: str, intermediate_file: str):
    """
    read the tsv file, skim for pathogenic entries, then write out to a new file

    Args:
        input_file ():
        intermediate_file ():
    """

    headers = ['chrom', 'pos', 'ref', 'alt', 'transcript', 'am_pathogenicity', 'am_class']

    # empty dictionary to contain the target indexes
    header_indexes: dict[str, int] = {}
    with gzip.open(input_file, 'rt') as read_handle, open(intermediate_file, 'wt') as write_handle:
        for line in read_handle:
            # skip over the headers
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # set the indexes (isoform and main-only have different columns)
                    header_indexes = process_header(line)
                continue

            if not header_indexes:
                raise ValueError('No header line was identified, columns are a mystery')

            # skip over everything except pathogenic
            if 'pathogenic' not in line:
                continue

            content = line.rstrip().split()

            # grab all the content
            content_dict: dict[str, str | float] = {key: content[header_indexes[key]] for key in headers}

            # trim transcripts
            assert isinstance(content_dict['transcript'], str)
            content_dict['transcript'] = content_dict['transcript'].split('.')[0]

            # convert the AM score to a float, and pos to an int
            content_dict['pos'] = int(content_dict['pos'])
            content_dict['am_pathogenicity'] = float(content_dict['am_pathogenicity'])
            write_handle.write(f'{json.dumps(content_dict)}\n')


def json_to_hail_table(json_file: str):
    """
    take a previously created JSON file and ingest it as a Hail Table
    requires an initiated Hail context

    Args:
        json_file ():
    """

    # define the schema for each written line
    schema = hl.dtype(
        'struct{'
        'chrom:str,'
        'pos:int32,'
        'ref:str,'
        'alt:str,'
        'transcript:str,'
        'am_pathogenicity:float64,'
        'am_class:str'
        '}',
    )

    # import as a hail table, force=True as this isn't Block-Zipped so all read on one core
    # We also provide a full attribute schema
    ht = hl.import_table(json_file, types={'f0': schema}, force=True, no_header=True)
    ht = ht.transmute(**ht.f0)

    # combine the two alleles into a single list
    ht = ht.transmute(locus=hl.locus(contig=ht.chrom, pos=ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.key_by('locus', 'alleles')
    ht.write(TABLE_DESTINATION)
    ht.describe()


def main():
    """
    takes the path to an AlphaMissense TSV, reorganises it into a Hail Table
    """

    init_batch()

    # get the AM file using Hail's hadoop open to read/write
    r = requests.get(AM_ZENODO, stream=True)
    with open('temp.tsv.gz', 'wb') as f:
        for chunk in r.raw.stream(1024, decode_content=False):
            if chunk:
                f.write(chunk)

    all_data_in_bytes = open('temp.tsv.gz', 'rb').read()

    # if it doesn't exist in GCP, push it there
    if not AnyPath(DESTINATION).exists():
        with hl.hadoop_open(DESTINATION, 'wb') as f:
            f.write(all_data_in_bytes)

    # generate a new tsv of just pathogenic entries
    filter_for_pathogenic_am('temp.tsv.gz', 'temp.json')

    # now ingest as HT and re-jig some fields
    json_to_hail_table('temp.json')


if __name__ == '__main__':
    main()
