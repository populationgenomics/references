#!/usr/bin/env python3
"""
Convert the deCODE 2019 sex-averaged GRCh38 genetic map to the two layouts that
`plink2 --cm-map` reads.

Provenance for decode2019_sexavg_GRCh38.eagle.txt.gz and
decode2019_sexavg_GRCh38_chr{1..22,X}.txt.gz, staged under
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

Output: the same rows in two gzipped, whitespace-separated layouts. Each has
one row per interval start, carrying the published cM at the previous
interval's end (back-calculated from rate for a chromosome's first interval),
plus a closing row at the last interval end.
  - decode2019_sexavg_GRCh38.eagle.txt.gz: one genome-wide file with a leading
    chromosome column (plink2's "Eagle-style" map):
        chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    Chromosome codes drop the 'chr' prefix; plink2 matches them to
    'chr'-prefixed datasets. Rows are in natural chromosome order (1..22, X),
    since plink2 rejects a map whose chromosomes are out of order.
  - decode2019_sexavg_GRCh38_chr{1..22,X}.txt.gz: one file per chromosome, with
    no chromosome column (the 3-column layout plink2 calls "SHAPEIT-format",
    from the SHAPEIT2/IMPUTE2 maps; SHAPEIT4/5 read a different layout):
        pposition rrate gposition

Using the map with plink2 --cm-map:
  - Whole-genome dataset: either layout; they give identical positions.
  - Dataset restricted to some chromosomes (--chr, or one chromosome per job):
    use the per-chromosome files through plink2's '@' pattern,
    --cm-map <dir>/decode2019_sexavg_GRCh38_chr@.txt.gz. plink2 only opens the
    files for chromosomes present. The genome-wide file fails there with
    "Chromosome ... is split", because it lists chromosomes the dataset lacks.
  - Variants before a chromosome's first map row get an extrapolated cM, which
    can be slightly negative; variants after its last row are clamped to the
    end value.

Gzip headers carry no timestamp, so re-running on the same input reproduces the
staged files byte for byte.

Fails on a missing header, a malformed row, an empty or inverted interval,
non-contiguous intervals, a decreasing genetic position, or a rate-derived
interval start that disagrees with the previous interval's published cM.

Usage:
    python convert_decode2019_map_for_plink2.py aau1043_datas3.gz <output dir>
"""

import argparse
import gzip
import io
import itertools
import os

CHROM_ORDER = [str(i) for i in range(1, 23)] + ['X']
OUTPUT_PREFIX = 'decode2019_sexavg_GRCh38'
EAGLE_HEADER = 'chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n'
PER_CHROM_HEADER = 'pposition rrate gposition\n'
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


def write_gzip_text(path: str, lines: list[str]) -> None:
    """
    Write text lines to a gzip file with no timestamp or name in its header.

    Args:
        path: Output path.
        lines: Lines to write, each ending in a newline.
    """
    with (
        open(path, 'wb') as raw,
        gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz,
        io.TextIOWrapper(gz, encoding='utf-8', newline='\n') as out,
    ):
        out.writelines(lines)


def main() -> None:
    """Parse arguments, convert the map, and write both gzipped layouts."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('input', help='aau1043_datas3 from the paper (plain or .gz)')
    parser.add_argument('output_dir', help='Existing directory to write the maps to')
    args = parser.parse_args()

    rows = to_eagle_rows(read_intervals(args.input))

    eagle_path = os.path.join(args.output_dir, f'{OUTPUT_PREFIX}.eagle.txt.gz')
    write_gzip_text(
        eagle_path,
        [EAGLE_HEADER]
        + [f'{chrom} {pos} {rate:.10g} {cm:.10g}\n' for chrom, pos, rate, cm in rows],
    )
    print(f'Wrote {len(rows)} rows for {len(CHROM_ORDER)} chromosomes to {eagle_path}')

    for chrom, chrom_rows in itertools.groupby(rows, key=lambda row: row[0]):
        chrom_path = os.path.join(args.output_dir, f'{OUTPUT_PREFIX}_chr{chrom}.txt.gz')
        write_gzip_text(
            chrom_path,
            [PER_CHROM_HEADER]
            + [f'{pos} {rate:.10g} {cm:.10g}\n' for _, pos, rate, cm in chrom_rows],
        )
    print(f'Wrote {len(CHROM_ORDER)} per-chromosome maps to {args.output_dir}')


if __name__ == '__main__':
    main()
