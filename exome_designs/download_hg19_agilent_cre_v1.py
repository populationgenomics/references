#!/usr/bin/env python3

"""
Download the Agilent SureSelect Clinical Research Exome v1 (S06588914)
Regions BED from UCSC at hg19 coordinates. Output:

  - S06588914_Regions_hg19.bed

CRE v1 is built on the SureSelectXT All Exon V5 backbone, which was specified
at exon resolution; Agilent ships Regions = Covered for this design by
construction (the portal Regions track header literally reads "This is same as
Covered.bed"), so only the Regions BED is downloaded. CREv2 (S30409818), a
ground-up redesign with sub-exon target resolution, keeps both files via
download_ucsc_exomes.py.

The hg19 BED is an intermediate artifact: upload it to cpg-common-test and
lift over to hg38 with liftover_and_convert_hg19_bedfiles.py.

Manufacturer-portal baseline (hg19, also distributed by Agilent in hg19 only)
for sanity-check comparison against the UCSC track downloaded here:
    /Users/jossch/Downloads/S06588914/S06588914_Regions.bed
"""

import urllib.parse as up
from pathlib import Path

import bbi
import polars as pl
import requests

API_URL = 'https://api.genome.ucsc.edu'
DL_URL = 'https://hgdownload.soe.ucsc.edu'
TRACK_ENDPOINT = '/list/tracks'
GENOME = 'hg19'
CANONICAL_CHRS = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'M']]
TARGET_DESIGN_ID = 'S06588914'
OUTPUT_FILENAMES = {
    'regions': 'S06588914_Regions_hg19.bed',
}


def get_exome_probesets(
    api_url: str, track_endpoint: str, genome: str,
) -> dict[str, dict]:
    """
    Fetch the UCSC ``exomeProbesets`` super-track listing for ``genome``.
    Returns the raw dict keyed by track name; values carry ``shortLabel``,
    ``longLabel`` and ``bigDataUrl`` among other fields.
    """

    params = {'genome': genome}
    response = requests.get(up.urljoin(api_url, track_endpoint), params)
    return response.json()[genome]['exomeProbesets']


def find_design_tracks(
    probesets: dict[str, dict], design_id: str,
) -> dict[str, str]:
    """
    Locate the Regions BigBed URL for the given Agilent design id
    (e.g. ``S06588914``) within the UCSC exomeProbesets track listing. Matches
    by BigBed filename stem (``<design>_Regions``) and cross-checks the track
    ``shortLabel`` for the design id.
    """

    matches: dict[str, str] = {}
    for track in probesets.values():
        # UCSC nests parent composite-track metadata (compositeContainer,
        # shortLabel, etc.) as plain strings alongside the per-track dicts;
        # skip those.
        if not isinstance(track, dict):
            continue
        url = track.get('bigDataUrl', '')
        if not url:
            continue
        stem = Path(url).stem
        in_url = design_id in stem
        in_label = design_id in track.get('shortLabel', '')
        if not (in_url or in_label):
            continue
        if stem.endswith('_Regions'):
            matches['regions'] = url
    missing = {'regions'} - matches.keys()
    if missing:
        raise RuntimeError(
            f'Could not locate {design_id} hg19 entries for: {sorted(missing)}',
        )
    return matches


def write_bed_from_bigbed(
    big_data_url: str, dl_url: str, output_path: Path,
    canonical_chrs: list[str],
) -> None:
    """
    Stream the BigBed at ``big_data_url`` (relative to ``dl_url``), keep only
    ``canonical_chrs``, and write to ``output_path`` as 3-column BED.
    """

    track = bbi.open(up.urljoin(dl_url, big_data_url))
    bed_df = pl.DataFrame(
        schema={'chrom': pl.String, 'start': pl.Int64, 'end': pl.Int64},
    )
    for chrom, size in track.chromsizes.items():
        if chrom not in canonical_chrs:
            continue
        chrom_bed = pl.from_pandas(
            track.fetch_intervals(chrom=chrom, start=1, end=size)[
                ['chrom', 'start', 'end']
            ],
        )
        bed_df = pl.concat([bed_df, chrom_bed])
    bed_df.write_csv(file=output_path, include_header=False, separator='\t')


def main() -> None:
    """
    Download the Regions hg19 BED for Agilent CRE v1 (S06588914).
    """

    probesets = get_exome_probesets(API_URL, TRACK_ENDPOINT, GENOME)
    urls = find_design_tracks(probesets, TARGET_DESIGN_ID)
    for kind, filename in OUTPUT_FILENAMES.items():
        write_bed_from_bigbed(
            urls[kind], DL_URL, Path(filename), CANONICAL_CHRS,
        )


if __name__ == '__main__':
    main()
