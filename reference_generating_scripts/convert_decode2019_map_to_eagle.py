#!/usr/bin/env python3
"""
Convert the deCODE 2019 sex-averaged GRCh38 genetic map to an Eagle-style map.

Provenance for decode2019_sexavg_GRCh38.eagle.txt.gz, staged under
gs://cpg-common-*/references/genetic_maps/decode2019/ via the
`genetic_maps_decode2019` Source in references.py.

Input: aau1043_datas3 (Data S3) from Halldorsson et al. (2019), "Characterizing
mutagenic effects of recombination through a sequence-level genetic map",
Science 363, eaau1043, doi:10.1126/science.aau1043, downloaded from the
paper's Supplementary Materials (https://doi.org/10.1126/science.aau1043). The
copy used is staged alongside the output as aau1043_datas3.gz (gzipped as
published). It is the average of the paternal and maternal maps, in GRCh38
coordinates, as intervals:
    Chr  Begin  End  cMperMb  cM
where cM is the genetic position at End. chrX carries the maternal rate (its
total is ~176 cM, the female X map length), which is the right scale for X IBD.

Output: one gzipped whitespace-separated file in the Eagle map layout that
`plink2 --cm-map` accepts in place of per-chromosome SHAPEIT maps:
    chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
with one row per interval start, carrying the published cM at the previous
interval's end (back-calculated from rate for a chromosome's first interval),
plus a closing row at the last interval end. Chromosome codes drop the 'chr'
prefix; plink2 matches them to 'chr'-prefixed datasets. Rows are in natural
chromosome order (1..22, X), since plink2 rejects a map whose chromosomes are
out of order ("Chromosome ... is split").

Using the map with plink2 --cm-map:
  - plink2 also reports "Chromosome ... is split" when the map lists chromosomes
    the dataset lacks, so a run restricted to some chromosomes (--chr, or one
    chromosome per job) must subset the map to those first, e.g.
    awk 'NR==1 || $1=="1"'.
  - Variants before a chromosome's first map row get an extrapolated cM, which
    can be slightly negative; variants after its last row are clamped to the
    end value.

Fails on a missing header, a malformed row, an empty or inverted interval,
non-contiguous intervals, a decreasing genetic position, or a rate-derived
interval start that disagrees with the previous interval's published cM.

Usage:
    python convert_decode2019_map_to_eagle.py aau1043_datas3 \
        decode2019_sexavg_GRCh38.eagle.txt.gz
"""

import argparse
import gzip
import itertools

CHROM_ORDER = [str(i) for i in range(1, 23)] + ['X']
HEADER = 'chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n'
# Allowed float error between a rate-derived interval start and the previous
# interval's published end cM (the published file agrees to within ~3e-14).
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
        lines = ((n, line) for n, line in enumerate(f, 1) if not line.startswith('#'))
        try:
            _, header_line = next(lines)
        except StopIteration:
            raise ValueError(f'No header row in {path}') from None
        header = header_line.split()
        if header != ['Chr', 'Begin', 'End', 'cMperMb', 'cM']:
            raise ValueError(f'Unexpected header in {path}: {header}')
        for lineno, line in lines:
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f'{path}:{lineno}: expected 5 fields, got {fields!r}')
            chrom, begin, end, rate, cm = fields
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
        ValueError: On an unknown chromosome, an empty or inverted interval, a gap
            or overlap between consecutive intervals, a decreasing genetic position,
            or an interval whose rate-derived start cM disagrees with the previous
            interval's published end cM.
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
            if end <= begin:
                raise ValueError(f'chr{chrom}: empty or inverted interval at {begin}')
            start_cm = cm - rate * (end - begin) / 1e6
            if cm < start_cm - CM_TOLERANCE:
                raise ValueError(f'chr{chrom}: genetic position decreases at {begin}')
            if prev_end is None:
                # The first interval may start at a nonzero cM; it only must not be negative.
                if start_cm < -CM_TOLERANCE:
                    raise ValueError(f'chr{chrom}: negative genetic position at {begin}')
                row_cm = max(start_cm, 0.0)
            else:
                if begin != prev_end:
                    raise ValueError(
                        f'chr{chrom}: interval starting {begin} does not follow {prev_end}'
                    )
                # cM is cumulative, so the start derived from this interval's rate must
                # match the published cM at the previous interval's end, either way.
                if abs(start_cm - prev_cm) > CM_TOLERANCE:
                    raise ValueError(
                        f'chr{chrom}: derived start cM {start_cm} at {begin} disagrees '
                        f'with the published {prev_cm}'
                    )
                row_cm = prev_cm
            rows.append((chrom, begin, rate, row_cm))
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
