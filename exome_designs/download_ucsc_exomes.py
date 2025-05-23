#!/usr/bin/env python3

import re
import unicodedata
import urllib.parse as up
from pathlib import Path

import bbi
import polars as pl
import requests

# CONSTs
API_URL = 'https://api.genome.ucsc.edu'
DL_URL = 'https://hgdownload.soe.ucsc.edu'
GENOME_ENDPOINT = '/list/ucscGenomes'
TRACK_ENDPOINT = '/list/tracks'
GENOME = 'hg38'
CANONICAL_CHRS = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'M']]
SOURCE_NAME = 'exome_probesets'
CPG_DEST = 'exome-probsets/hg38/'


def slug_and_under(line: str):
    """
    Slugify and replace hyphen with underscore string.

    Example:
    >>> slug_and_under(u'Héllø W.1')
    'hello_w_1'
    """

    line = unicodedata.normalize('NFKD', line).encode('ascii', 'ignore').decode()
    line = line.strip().lower()
    line = re.sub(
        r'[\s.]+',
        '-',
        line,
    )
    line = re.sub(r'[\s-]+', '_', line)
    line = line.strip('_')
    return line


def reformat_endpoint_dict(
    end_point_dict: dict[str, dict],
) -> list[tuple[str, str, str]]:
    """
    Given a dict from json response from UCSC track endpoint, returns
    a list of tuples ('long name', 'short name', 'bigDataUrl download link')
    """

    bbs_list = []
    for k in end_point_dict.keys():
        # keep only those with a download link
        if 'bigDataUrl' in end_point_dict[k]:
            bbs_list.append(
                (
                    slug_and_under(end_point_dict[k]['longLabel']),
                    slug_and_under(end_point_dict[k]['shortLabel']),
                    end_point_dict[k]['bigDataUrl'],
                )
            )
    return bbs_list


def get_exome_bbs(api_url: str, track_endpoint: str, genome: str) -> list[tuple[str, str, str]]:
    """
    Construct a query for the UCSC API. Submit it, and then parse the response
    with call to reformat_endpoint_dict()
    """

    params = {'genome': genome}
    response = requests.get(up.urljoin(api_url, track_endpoint), params)
    bbs = reformat_endpoint_dict(response.json()[genome]['exomeProbesets'])
    return bbs


def make_beds_from_bbs(bbs: list[tuple[str, str, str]], dl_url: str, canonical_chrs: list[str]) -> None:
    """
    Iteate over list of retreived BigBed file names + locations.
    Construct full dl link, rerieve using bbi, convert to bed and write to file.
    """

    for bb in bbs:
        tmp = bbi.open(up.urljoin(dl_url, bb[2]))
        bed_df = pl.DataFrame(schema={'chrom': pl.String, 'start': pl.Int64, 'end': pl.Int64})
        for c, s in tmp.chromsizes.items():
            if c in canonical_chrs:
                c_bed = pl.from_pandas(tmp.fetch_intervals(chrom=c, start=1, end=s)[['chrom', 'start', 'end']])
                bed_df = pl.concat([bed_df, c_bed])
        bed_df.write_csv(file=Path(bb[2]).stem + '.bed', include_header=False, separator='\t')


def make_reference_source(source_name: str, destination_path: str, bbs: list[tuple[str, str, str]]) -> None:
    """
    Create a 'Source' string to be copied into references.py, and write it to file
    """

    sp_tab = '    '
    open_list = [
        'Source(',
        f'{sp_tab}\'{source_name}\',',
        f'{sp_tab}# exome probset defintions (bed file and interval_list) format\n{sp_tab}# downloaded from UCSC with the download_ucsc_exomes.py script',
        f'{sp_tab}\'dst={destination_path}\',',
        ]
    source_open = '\n'.join(open_list)
    source_close = '\n),'
    files = build_files_list(bbs)
    smp_source_str = source_open + files + source_close
    with open('ucsc_exome_probsets_source_entry.txt', 'w') as write_handle:
        write_handle.write(smp_source_str)
    

def build_files_list(bbs: list[tuple[str, str, str]]) -> str:
    """
    Iterate over list of retreived BigBed file names + locations.
    Construct a txt str dict of keys (human readable names) and values ( actual bed / intervals )
    
    example:
    'files=dict(
            vcf='twist_exome_benchmark_truth.vcf.gz',
            index='twist_exome_benchmark_truth.vcf.gz.tbi',
            bed='Twist_Exome_Core_Covered_Targets_hg38.bed',
    ),'
    
    """
    sp_tab = '    '
    files_open = f'\n{sp_tab*2}files=dict(\n'
    files_close = f',\n{sp_tab*2}),'
    files_body = []

    for bb in bbs:
        files_body.append(f'{sp_tab*3}{Path(bb[0]).stem + "_bed"}=\'{Path(bb[2]).stem + ".bed"}\'')
        files_body.append(f'{sp_tab*3}{Path(bb[0]).stem + "_interval_list"}=\'{Path(bb[2]).stem + ".interval_list"}\'')
    return (files_open + ',\n'.join(files_body) + files_close)
    

def main():
    """
    Retrieves names and links for BigBed formatted Exome probesets from UCSC
    Converts from BigBed to Bed and writes to file.
    Bed output keeps same prefix as BigBed original.
    """

    exome_dl = get_exome_bbs(API_URL, TRACK_ENDPOINT, GENOME)
    make_beds_from_bbs(exome_dl, DL_URL, CANONICAL_CHRS)
    make_reference_source(SOURCE_NAME, CPG_DEST, exome_dl)


if __name__ == '__main__':
    main()
