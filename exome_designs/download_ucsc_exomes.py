#!/usr/bin/env python3


import urllib.parse as up
from pathlib import Path

import bbi
import polars as pl
import requests
from slugify import slugify

# CONSTs
API_URL = 'https://api.genome.ucsc.edu'
DL_URL = 'https://hgdownload.soe.ucsc.edu'
GENOME_ENDPOINT = '/list/ucscGenomes'
TRACK_ENDPOINT = '/list/tracks'
GENOME = 'hg38'
CANONICAL_CHRS = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'M']]


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
                    slugify(end_point_dict[k]['longLabel']),
                    slugify(end_point_dict[k]['shortLabel']),
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


def main():
    """
    Retrieves names and links for BigBed formatted Exome probesets from UCSC
    Converts from BigBed to Bed and writes to file.
    Bed output keeps same prefix as BigBed original.
    """

    exome_dl = get_exome_bbs(API_URL, TRACK_ENDPOINT, GENOME)
    make_beds_from_bbs(exome_dl, DL_URL, CANONICAL_CHRS)


if __name__ == '__main__':
    main()
