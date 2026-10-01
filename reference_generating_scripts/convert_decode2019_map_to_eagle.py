#!/usr/bin/env python3
"""
Convert the deCODE 2019 sex-averaged GRCh38 genetic map to an Eagle-style map.

Provenance for decode2019_sexavg_GRCh38.eagle.txt.gz, staged under
gs://cpg-common-*/references/genetic_maps/decode2019/ via the
`genetic_maps_decode2019` Source in references.py.

Input: aau1043_datas3 (Data S3) from Halldorsson et al. (2019), "Characterizing
mutagenic effects of recombination through a sequence-level genetic map",
Science 363, eaau1043, doi:10.1126/science.aau1043. It is the average of the
paternal and maternal maps, in GRCh38 coordinates, as intervals:
    Chr  Begin  End  cMperMb  cM
where cM is the genetic position at End. chrX carries the maternal rate (its
total is ~176 cM, the female X map length), which is the right scale for X IBD.

Output: one gzipped whitespace-separated file in the Eagle map layout that
`plink2 --cm-map` accepts in place of per-chromosome SHAPEIT maps:
    chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
with one row per interval start (genetic position back-calculated from the
interval's end position and rate) plus a closing row at the last interval end.
Chromosome codes drop the 'chr' prefix; plink2 matches them to 'chr'-prefixed
datasets. Rows are in natural chromosome order (1..22, X), since plink2 rejects
a map whose chromosomes are out of order ("Chromosome ... is split").

Fails on non-contiguous intervals or a decreasing genetic position.

Usage:
    python convert_decode2019_map_to_eagle.py aau1043_datas3 \
        decode2019_sexavg_GRCh38.eagle.txt.gz
"""

import argparse
import gzip
import itertools

CHROM_ORDER = [str(i) for i in range(1, 23)] + ['X']
HEADER = 'chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n'
# Back-calculated start positions can dip below the previous end by float error.
CM_TOLERANCE = 1e-6


def read_intervals(path: str) -> list[tuple[str, int, int, float, float]]:
    """
    Read the deCODE interval map, skipping comment lines and the header row.

    Args:
        path: Path to aau1043_datas3 (plain or gzipped).

    Returns:
        (chrom without 'chr', begin, end, cM/Mb, cM at end) per interval, in file order.
    """
    opener = gzip.open if path.endswith('.gz') else open
    intervals = []
    with opener(path, 'rt') as f:
        lines = (line for line in f if not line.startswith('#'))
        header = next(lines).split()
        if header != ['Chr', 'Begin', 'End', 'cMperMb', 'cM']:
            raise ValueError(f'Unexpected header in {path}: {header}')
        for line in lines:
            chrom, begin, end, rate, cm = line.split()
            intervals.append(
                (chrom.removeprefix('chr'), int(begin), int(end), float(rate), float(cm))
            )
    return intervals


def to_eagle_rows(
    intervals: list[tuple[str, int, int, float, float]],
) -> list[tuple[str, int, float, float]]:
    """
    Turn intervals into Eagle map rows, validating contiguity and monotonic cM.

    Args:
        intervals: Output of read_intervals.

    Returns:
        (chrom, position, cM/Mb, cM) rows in natural chromosome order.

    Raises:
        ValueError: On an unknown chromosome, a gap or overlap between consecutive
            intervals, or a genetic position that decreases along a chromosome.
    """
    by_chrom: dict[str, list[tuple[str, int, int, float, float]]] = {}
    for chrom, group in itertools.groupby(intervals, key=lambda iv: iv[0]):
        if chrom not in CHROM_ORDER:
            raise ValueError(f'Unexpected chromosome {chrom!r}')
        if chrom in by_chrom:
            raise ValueError(f'Chromosome {chrom} appears in more than one block')
        by_chrom[chrom] = list(group)

    rows: list[tuple[str, int, float, float]] = []
    for chrom in CHROM_ORDER:
        if chrom not in by_chrom:
            raise ValueError(f'Chromosome {chrom} missing from the map')
        prev_end, prev_cm = None, 0.0
        for _, begin, end, rate, cm in by_chrom[chrom]:
            if prev_end is not None and begin != prev_end:
                raise ValueError(
                    f'chr{chrom}: interval starting {begin} does not follow {prev_end}'
                )
            start_cm = cm - rate * (end - begin) / 1e6
            if start_cm < prev_cm - CM_TOLERANCE or cm < start_cm - CM_TOLERANCE:
                raise ValueError(f'chr{chrom}: genetic position decreases at {begin}')
            rows.append((chrom, begin, rate, max(start_cm, prev_cm)))
            prev_end, prev_cm = end, cm
        rows.append((chrom, prev_end, 0.0, prev_cm))
    return rows


def main() -> None:
    """Parse arguments, convert the map, and write the gzipped Eagle-style file."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('input', help='aau1043_datas3 from the paper (plain or .gz)')
    parser.add_argument('output', help='Output path, ending .txt.gz')
    args = parser.parse_args()

    rows = to_eagle_rows(read_intervals(args.input))
    with gzip.open(args.output, 'wt') as out:
        out.write(HEADER)
        for chrom, pos, rate, cm in rows:
            out.write(f'{chrom} {pos} {rate:.10g} {cm:.10g}\n')
    print(f'Wrote {len(rows)} rows for {len(CHROM_ORDER)} chromosomes to {args.output}')


if __name__ == '__main__':
    main()
