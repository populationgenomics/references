#!/usr/bin/env python3

"""
Pull the associations worth analysing out of the formatted GWAS summary statistics.

`format.py` writes every row a study published, up to ~200 million per file and 264 GB
in total. Almost all of it is noise: an analysis that ranks variants by effect only
needs the ones that reached significance. Reading all 336 files to find them takes
about a day in one process, which is why that step does not belong in a notebook.

This script does it once, one Batch job per file, and writes the result next to the
formatted files. Downstream work then reads a few hundred megabytes in seconds.

Two outputs, both under `v1/`:

    hits/<file_id>.parquet      associations with p <= THRESHOLD, all columns kept
    chromosome_coverage.csv     which chromosomes each file contains, and whether the
                                source published chromosome X

The coverage table matters as much as the hits. 95 of the 336 files have no
chromosome X, including every Million Veteran Program file, because the authors
never published it. Without the table a biomarker sourced only from those files
looks cleanly negative for X-linked causes such as G6PD, when in truth it was never
tested. Coverage is read from each file's tabix index, so it needs no data pass.

X rows alone do not prove the source had X. Liftover from GRCh37 or NCBI36 moves a
few autosome and Y positions onto GRCh38 X (e.g. chr7:100614152), and format.py
keeps them. Those can only land in the stretches of X that the chain file maps from
other chromosomes: about 9.6 kb, in 20 of X's ~9,500 16-kb windows, for GRCh37. So
a file counts as having X only if it has X data in a window outside those. The index
records which windows hold data, so this too needs no data pass.

Filtering uses neg_log_10_p_value where the file has it, falling back to parsing
p_value. The published P value is kept as text by format.py precisely because the
smallest ones (e.g. 8.51e-3553) are below the float range, and those are the
strongest associations: reading them as numbers would silently turn them into zero.
A published P of exactly 0, or an infinite published -log10 P, is NA in both columns
of the formatted file (format.py counts it as p_value_zero), yet those are the
strongest hits of all, e.g. APOE in Timsina 2026. So a row with no readable P value
is kept when its z (or beta / standard_error) reaches the same threshold; its P
columns stay NA. A row with no P value and no effect size cannot be recovered.

The threshold and row counts are stored in each parquet's metadata. A file is
redone when its formatted source is newer than its hits, so re-formatting a file
and running this again refreshes it.

Every column is written as text, as format.py wrote it, with NA as null rather than
the string 'NA', so files with and without hits share one schema and the folder
reads back as a single table. Cast on reading.

Both commands write to main, so run them from `main` after merge, once every
format.py job has finished. On Hail Batch, one job per file; the analysis-runner job
only submits them, and prints the batch link to follow:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Extract significant GWAS associations" \
        python3 gwas_sumstats/extract.py

The coverage table needs only the formatted files, not the hits, and is quick
enough to run in the analysis-runner job itself:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "GWAS chromosome coverage" \
        python3 gwas_sumstats/extract.py --coverage

Locally, on a file formatted locally (for testing; see the README for setup):

    python3 gwas_sumstats/extract.py --local --formatted ~/gwas_sumstats_local/v1 \
        --out ~/gwas_sumstats_local/hits \
        --only 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951
"""

import argparse
import csv
import gzip
import io
import json
import math
import shlex
import struct
from pathlib import Path
from statistics import NormalDist

FILES_CSV = Path(__file__).with_name('files.csv')
FORMATTED = 'gs://cpg-common-main/references/gwas_sumstats/v1'
DEFAULT_OUT = f'{FORMATTED}/hits'
DEFAULT_COVERAGE = f'{FORMATTED}/chromosome_coverage.csv'

THRESHOLD = 1e-5
"""
Deliberately looser than genome-wide significance (5e-8).

Whoever uses these files chooses their own threshold; this one only has to be loose
enough that they never need to come back and re-read the originals. Under the null it
keeps 200 times more rows than 5e-8 while still discarding over 99% of each file.
"""

TABIX_PSEUDO_BIN = 37450
"""Bin in which htslib records each contig's row count in a tabix index."""

WINDOW_SHIFT = 14
"""Tabix's linear index has one entry per 2**14 = 16,384 bp window."""

LIFTED_BUILDS = {'GRCh37': 'liftover_37_to_38', 'NCBI36': 'liftover_36_to_38'}
"""Source builds format.py lifts, and the references key of each chain file."""

CHAIN_SOURCES_OFF_X = {f'chr{i}' for i in range(1, 23)} | {'chrY'}
"""Source contigs format.py lifts (its CHAIN_CONTIGS) other than X."""

CHUNK_ROWS = 2_000_000
PIP_PACKAGES = 'pandas==2.2.3 pyarrow==18.1.0'


def formatted_name(file_id: str) -> str:
    """Name format.py gave this file: <file_id>_GRCh38_formatted.tsv.gz.

    The manifest's own `formatted_file` column cannot be used: every row of it reads
    'out.tsv.gz', the name each Batch job used before its output was copied out.
    """
    return f'{file_id}_GRCh38_formatted.tsv.gz'


def passing(p_text, neg_log_10, cutoff: float, z=None, beta=None, se=None):
    """
    Which associations are at or below the cutoff, one boolean per row.

    Prefers neg_log_10_p_value, which format.py computed from the published text and
    which therefore survives P values below the float range; inf, where a study
    published -log10 P as inf, passes. Rows without a readable
    neg_log_10_p_value fall back to parsing p_value, where text such as 8.51e-3553
    reads as 0 and so also passes.

    A row with neither is what format.py writes for a published P of exactly 0, the
    strongest hits of all. It is kept when |z|, or |beta / standard_error| where z
    is NA, reaches the two-sided z for the cutoff (4.42 for 1e-5). A row with no P
    value and no effect size is dropped: it cannot be ranked.

    Args:
        p_text: the p_value column as written, a pandas Series of text
        neg_log_10: the neg_log_10_p_value column, text, possibly empty or NA
        cutoff: the P value threshold
        z, beta, se: the z, beta and standard_error columns, text; None if absent
    """
    import pandas as pd

    def number(column):
        if column is None:
            return pd.Series(float('nan'), index=p_text.index)
        return pd.to_numeric(column, errors='coerce')

    neg_log = number(neg_log_10)
    p_value = number(p_text)
    z_score = number(z)
    z_score = z_score.where(z_score.notna(), number(beta) / number(se))
    by_z = z_score.abs() >= NormalDist().inv_cdf(1 - cutoff / 2)
    by_p_value = (p_value <= cutoff).where(p_value.notna(), by_z)
    by_neg_log = neg_log >= -math.log10(cutoff)
    return by_neg_log.where(neg_log.notna(), by_p_value).fillna(False).astype(bool)


def extract_one(source: Path, dest: Path, file_id: str, cutoff: float) -> dict:
    """
    Write the associations at or below `cutoff` from one formatted file.

    Read in chunks so memory follows the chunk, not the file: the largest of these
    holds ~200 million rows. The threshold and counts go into the parquet's metadata
    under b'gwas_sumstats', so each file records how it was made.

    Args:
        source: a formatted .tsv.gz
        dest: parquet to write
        file_id: stamped onto every row, so files can be concatenated
        cutoff: the P value threshold

    Returns:
        counts for the run record
    """
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    kept: list[pd.DataFrame] = []
    columns: list[str] = []
    rows_in = 0
    with gzip.open(source, 'rt') as handle:
        for chunk in pd.read_csv(handle, sep='\t', dtype=str, chunksize=CHUNK_ROWS):
            rows_in += len(chunk)
            columns = list(chunk.columns)
            keep = passing(
                chunk['p_value'],
                chunk.get('neg_log_10_p_value'),
                cutoff,
                z=chunk.get('z'),
                beta=chunk.get('beta'),
                se=chunk.get('standard_error'),
            )
            if keep.any():
                kept.append(chunk[keep])

    # A file with no hits still gets a parquet, with the right columns and no rows,
    # so that reading the whole folder back never has to special-case it. The explicit
    # schema keeps those columns text: left to infer, pyarrow types an empty column as
    # null, and a folder mixing null and text columns fails to read.
    hits = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=columns)
    hits.insert(0, 'file_id', file_id)
    schema = pa.schema([(column, pa.string()) for column in hits.columns])
    counts = {
        'file_id': file_id,
        'rows_in': rows_in,
        'rows_kept': len(hits),
        'rows_kept_without_p': int(hits['p_value'].isna().sum()),
    }
    table = pa.Table.from_pandas(hits, schema=schema, preserve_index=False)
    record = json.dumps({**counts, 'threshold': cutoff}).encode()
    table = table.replace_schema_metadata(
        {**(table.schema.metadata or {}), b'gwas_sumstats': record}
    )
    pq.write_table(table, dest)
    return counts


def read_tabix_index(handle) -> dict[str, dict]:
    """
    Rows and occupied 16-kb windows on each contig, read from a tabix (.tbi) index.

    The index is BGZF, which gzip reads. Its header names every contig with data.
    For each contig htslib records the row count in a pseudo-bin, and a linear index
    holding, per window, the file offset of its first row; an empty window repeats
    the offset of the next one, so the windows with data are those whose offset
    differs from the next, plus the last. Parsed here rather than with pysam so the
    coverage run needs no extra package, and no access token for htslib to open
    gs:// paths.

    Args:
        handle: the .tbi opened in binary mode

    Returns:
        contig -> {'rows': row count, 'windows': set of window indices with data}
    """
    with gzip.GzipFile(fileobj=handle) as index:
        data = index.read()
    if data[:4] != b'TBI\x01':
        raise ValueError(f'not a tabix index: starts {data[:4]!r}')
    # n_ref, format, col_seq, col_beg, col_end, meta, skip, l_nm
    header = struct.unpack_from('<8i', data, 4)
    names = [name.decode() for name in data[36 : 36 + header[7]].split(b'\0') if name]
    offset = 36 + header[7]
    contigs = {}
    for name in names:
        rows = None
        (n_bin,) = struct.unpack_from('<i', data, offset)
        offset += 4
        for _ in range(n_bin):
            bin_number, n_chunk = struct.unpack_from('<Ii', data, offset)
            offset += 8
            if bin_number == TABIX_PSEUDO_BIN:
                # (file start, file end), then (rows mapped, rows unmapped)
                rows = struct.unpack_from('<4Q', data, offset)[2]
            offset += 16 * n_chunk
        (n_intv,) = struct.unpack_from('<i', data, offset)
        starts = struct.unpack_from(f'<{n_intv}Q', data, offset + 4)
        offset += 4 + 8 * n_intv
        if rows is None:
            raise ValueError(f'tabix index has no row count for {name}')
        windows = {
            i for i in range(n_intv) if i == n_intv - 1 or starts[i] != starts[i + 1]
        }
        contigs[name] = {'rows': rows, 'windows': windows}
    return contigs


def chromosomes_in(path: str) -> dict[str, dict]:
    """
    Rows and occupied windows per chromosome in a formatted file, from its index.

    Chromosomes are GWAS-SSF numeric, so X is '23'.

    Args:
        path: a formatted .tsv.gz (its .tbi must sit beside it)
    """
    from cpg_utils import to_path

    with to_path(f'{path}.tbi').open('rb') as handle:
        return read_tabix_index(handle)


def liftover_windows_on_x(chain_text) -> set[int]:
    """
    GRCh38 X windows that liftover can fill from source chromosomes other than X.

    A UCSC chain maps source (t) blocks to GRCh38 (q) blocks. Every block from a
    lifted source contig other than X that lands on chrX marks where a stray row can
    appear; rows anywhere else on X came from the source's own X.

    Args:
        chain_text: lines of a .over.chain file, decompressed
    """
    windows: set[int] = set()
    source = target = strand = None
    target_size = position = 0
    for line in chain_text:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'chain':
            # score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            source, target, strand = fields[2], fields[7], fields[9]
            target_size, position = int(fields[8]), int(fields[10])
            continue
        size = int(fields[0])
        if target == 'chrX' and source in CHAIN_SOURCES_OFF_X:
            start, end = position, position + size
            if strand == '-':
                start, end = target_size - end, target_size - start
            windows.update(
                range(start >> WINDOW_SHIFT, ((end - 1) >> WINDOW_SHIFT) + 1)
            )
        if len(fields) == 3:
            position += size + int(fields[2])
    return windows


def coverage_record(
    row: dict, contigs: dict[str, dict] | None, liftover_windows: set[int]
) -> dict:
    """
    One row of the coverage table.

    Args:
        row: the files.csv row
        contigs: read_tabix_index of the formatted file; None if not formatted yet
        liftover_windows: X windows liftover can fill from other chromosomes, for
            this file's source build; empty for GRCh38 sources, which are not lifted
    """
    record = {
        'file_id': row['file_id'],
        'biomarker': row['biomarker'],
        'ancestry': row['ancestry'],
        'source_build': row['source_build'],
        'formatted': contigs is not None,
    }
    if contigs is None:
        empty = ['has_x', 'x_rows', 'x_windows', 'x_windows_from_source_x']
        return {**record, **dict.fromkeys(empty, ''), 'chromosomes': ''}
    x = contigs.get('23', {'rows': 0, 'windows': set()})
    from_source_x = len(x['windows'] - liftover_windows)
    return {
        **record,
        'has_x': from_source_x > 0,
        'x_rows': x['rows'],
        'x_windows': len(x['windows']),
        'x_windows_from_source_x': from_source_x,
        'chromosomes': ' '.join(contigs),
    }


def run_coverage(
    rows: list[dict], out: str, formatted: str, chains: dict[str, str]
) -> None:
    """
    Write one row per file saying which chromosomes it has, and whether the source
    published X.

    Files not formatted yet are listed with formatted=False, so a partial run is
    visible in the table itself.

    Args:
        rows: files.csv rows
        out: destination CSV
        formatted: folder holding the formatted files
        chains: source build -> chain file format.py lifted it with
    """
    from cpg_utils import to_path

    liftover_windows = {'GRCh38': set()}
    for build, chain in chains.items():
        with gzip.open(to_path(chain).open('rb'), 'rt') as handle:
            liftover_windows[build] = liftover_windows_on_x(handle)

    records = []
    for row in rows:
        path = f'{formatted}/{formatted_name(row["file_id"])}'
        contigs = chromosomes_in(path) if to_path(path).exists() else None
        if contigs is None:
            print(f'{row["file_id"]}: not formatted')
        windows = liftover_windows[row['source_build']]
        records.append(coverage_record(row, contigs, windows))
    done = [record for record in records if record['formatted']]
    if not done:
        raise SystemExit(
            f'none of the {len(rows)} selected files is formatted in {formatted}'
        )

    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
    with to_path(out).open('w') as handle:
        handle.write(text.getvalue())

    with_x = sum(record['has_x'] for record in done)
    print(
        f'{len(done)} files: {with_x} with chromosome X, {len(done) - with_x} without; '
        f'{len(records) - len(done)} not formatted'
    )
    print(f'written to {out}')


def run_local(
    rows: list[dict], out: Path, formatted: str, cutoff: float, force: bool = False
) -> None:
    """
    Extract from local or remote files one after another, into a local folder.

    Local formatted files are read in place. Remote ones are downloaded to a folder
    beside `out`, not into it, so `out` holds only parquets and reads as one table.
    """
    out.mkdir(parents=True, exist_ok=True)
    downloads = out.parent / f'{out.name}_formatted'
    for row in rows:
        dest = out / f'{row["file_id"]}.parquet'
        if dest.exists() and not force:
            print(f'{row["file_id"]}: exists, skipped')
            continue
        source = f'{formatted}/{formatted_name(row["file_id"])}'
        if source.startswith('gs://'):
            from cpg_utils import to_path

            local = downloads / formatted_name(row['file_id'])
            if not local.exists():
                downloads.mkdir(parents=True, exist_ok=True)
                to_path(source).download_to(local)
        else:
            local = Path(source).expanduser()
        counts = extract_one(local, dest, row['file_id'], cutoff)
        print(
            f'{counts["file_id"]}: {counts["rows_kept"]:,} of {counts["rows_in"]:,} rows'
        )


def is_current(dest, source) -> bool:
    """
    Whether hits already written are newer than the formatted file they came from.

    After format.py --only <id> --force, the formatted file is newer than its hits,
    so the hits are redone rather than skipped.

    Args:
        dest: the hits parquet
        source: the formatted file
    """
    return dest.exists() and dest.stat().st_mtime >= source.stat().st_mtime


def run_batch(rows: list[dict], args) -> None:
    """
    One Hail Batch job per formatted file not yet extracted, or re-formatted since.

    Each job reads one file and writes one parquet, so 336 files take as long as the
    largest single file rather than the sum of all of them.
    """
    from cpg_utils import to_path
    from cpg_utils.config import config_retrieve
    from cpg_utils.hail_batch import get_batch

    batch = get_batch(name='Extract significant GWAS associations')
    script = Path(__file__).read_text()
    submitted, missing = 0, []
    for row in rows:
        source = to_path(f'{args.formatted}/{formatted_name(row["file_id"])}')
        dest = to_path(f'{args.out}/{row["file_id"]}.parquet')
        if not source.exists():
            missing.append(row['file_id'])
            continue
        if not args.force and is_current(dest, source):
            print(f'{row["file_id"]}: exists, skipped')
            continue

        size_gib = int(row['size_bytes'] or 0) / 2**30
        job = batch.new_bash_job(f'extract {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2)
        job.storage(f'{int(size_gib * 3) + 10}Gi')
        job.command(
            f'python3 -m pip install --quiet {PIP_PACKAGES}\n'
            f"cat > extract.py <<'GWAS_SUMSTATS_SCRIPT'\n{script}\nGWAS_SUMSTATS_SCRIPT\n"
            f'python3 extract.py --file-id {shlex.quote(row["file_id"])} '
            f'--source {batch.read_input(str(source))} --dest {job.out} '
            f'--threshold {args.threshold}'
        )
        batch.write_output(job.out, str(dest))
        submitted += 1

    print(f'{submitted} extract jobs submitted')
    if missing:
        print(
            f'{len(missing)} not formatted yet, so skipped:\n  ' + '\n  '.join(missing)
        )
    if submitted:
        batch.run(wait=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument(
        '--files', type=Path, default=FILES_CSV, help='files.csv to read'
    )
    parser.add_argument('--formatted', default=FORMATTED, help='formatted files folder')
    parser.add_argument('--out', default=DEFAULT_OUT, help='destination folder')
    parser.add_argument(
        '--coverage-out',
        help=f'coverage CSV to write (default {DEFAULT_COVERAGE}); required with --only',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=THRESHOLD,
        help=f'P value (not -log10 P) to keep at or below (default {THRESHOLD}); '
        'anything else needs its own --out',
    )
    parser.add_argument('--only', nargs='+', help='file_ids to extract')
    for build, key in LIFTED_BUILDS.items():
        parser.add_argument(
            f'--chain-{build.lower()}',
            help=f'coverage: chain format.py used; default references {key}',
        )
    parser.add_argument('--local', action='store_true', help='run here, not on Batch')
    parser.add_argument(
        '--coverage', action='store_true', help='write the chromosome coverage table'
    )
    parser.add_argument(
        '--force', action='store_true', help='redo files already written; needs --only'
    )
    parser.add_argument('--file-id', help='job mode: file_id stamped on every row')
    parser.add_argument('--source', type=Path, help='job mode: input file')
    parser.add_argument('--dest', type=Path, help='job mode: output parquet')
    args = parser.parse_args()

    if args.file_id:
        counts = extract_one(args.source, args.dest, args.file_id, args.threshold)
        print(
            f'{counts["file_id"]}: {counts["rows_kept"]:,} of {counts["rows_in"]:,} '
            f'rows, {counts["rows_kept_without_p"]:,} of them kept on z'
        )
        return

    if args.force and not args.only:
        parser.error('--force needs --only, so every file is not redone by mistake')
    if args.coverage and args.only and not args.coverage_out:
        parser.error(
            '--coverage with --only would replace chromosome_coverage.csv with a '
            'subset; give an output path: --coverage-out subset.csv'
        )
    if args.threshold != THRESHOLD and args.out == DEFAULT_OUT:
        parser.error(
            f'--threshold {args.threshold} would mix with the {THRESHOLD} hits in '
            f'{DEFAULT_OUT}; give another --out'
        )

    from download import read_files_csv

    rows = read_files_csv(args.files, args.only)
    if args.coverage:
        from cpg_utils.config import reference_path

        chains = {
            build: getattr(args, f'chain_{build.lower()}') or str(reference_path(key))
            for build, key in LIFTED_BUILDS.items()
        }
        run_coverage(
            rows, args.coverage_out or DEFAULT_COVERAGE, args.formatted, chains
        )
    elif args.local:
        if args.out.startswith('gs://'):
            parser.error(f'--local needs a local --out folder, got {args.out}')
        run_local(rows, Path(args.out), args.formatted, args.threshold, args.force)
    else:
        run_batch(rows, args)


if __name__ == '__main__':
    main()
