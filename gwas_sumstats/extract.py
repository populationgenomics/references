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
    chromosome_coverage.csv     which chromosomes each file contains

The coverage table matters as much as the hits. 95 of the 336 files have no
chromosome X, including every Million Veteran Program file, because the authors
never published it. Without the table a biomarker sourced only from those files
looks cleanly negative for X-linked causes such as G6PD, when in truth it was never
tested. Coverage is read from each file's tabix index, so it needs no data pass.

Filtering uses neg_log_10_p_value where the file has it, falling back to parsing
p_value. The published P value is kept as text by format.py precisely because the
smallest ones (e.g. 8.51e-3553) are below the float range, and those are the
strongest associations: reading them as numbers would silently turn them into zero.

On Hail Batch, one job per file:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Extract significant GWAS associations" \
        python3 gwas_sumstats/extract.py

Then the coverage table, which is quick enough to run in the driver:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "GWAS chromosome coverage" \
        python3 gwas_sumstats/extract.py --coverage

Locally, on a few files (for testing):

    python3 gwas_sumstats/extract.py --local --out ./hits \
        --only 2017_Wheeler_PLoSMed_HbA1c_AFR_GCST007952
"""

import argparse
import csv
import gzip
import io
import json
import math
import shlex
from pathlib import Path

FILES_CSV = Path(__file__).with_name('files.csv')
FORMATTED = 'gs://cpg-common-main/references/gwas_sumstats/v1'
DEFAULT_OUT = f'{FORMATTED}/hits'
DEFAULT_COVERAGE = f'{FORMATTED}/chromosome_coverage.csv'

THRESHOLD = 1e-5
"""
Deliberately looser than genome-wide significance (5e-8).

Whoever uses these files chooses their own threshold; this one only has to be loose
enough that they never need to come back and re-read the originals. It keeps roughly
a thousand times more rows than 5e-8 while still discarding over 99% of each file.
"""

CHUNK_ROWS = 2_000_000
PIP_PACKAGES = 'pandas==2.2.3 pyarrow==18.1.0 pysam==0.22.1'


def formatted_name(file_id: str) -> str:
    """Name format.py gave this file: <file_id>_GRCh38_formatted.tsv.gz.

    The manifest's own `formatted_file` column cannot be used: every row of it reads
    'out.tsv.gz', the name each Batch job used before its output was copied out.
    """
    return f'{file_id}_GRCh38_formatted.tsv.gz'


def read_files_csv(path: Path, only: list[str] | None = None) -> list[dict]:
    """
    Rows of files.csv, optionally restricted to the given file_ids.

    Args:
        path: files.csv
        only: file_ids to keep; all rows if None
    """
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    if only:
        missing = set(only) - {row['file_id'] for row in rows}
        if missing:
            raise ValueError(f'not in {path.name}: {sorted(missing)}')
        rows = [row for row in rows if row['file_id'] in only]
    return rows


def passes(p_text: str, neg_log_10: str, cutoff: float) -> bool:
    """
    Is this association at or below the cutoff?

    Prefers neg_log_10_p_value, which format.py computed from the published text and
    which therefore survives P values below the float range. Falls back to parsing
    p_value. A row whose P value cannot be read either way is dropped: it cannot be
    ranked, and format.py already fails any file where more than 1% are unreadable.

    A p_value of exactly 0 is kept. It means the source underflowed, so the true
    value is smaller than a float can hold, which makes it one of the strongest
    associations in the file rather than one to discard.

    Args:
        p_text: the p_value field as written
        neg_log_10: the neg_log_10_p_value field, possibly empty or NA
        cutoff: the P value threshold
    """
    if neg_log_10 and neg_log_10 != 'NA':
        try:
            return float(neg_log_10) >= -math.log10(cutoff)
        except ValueError:
            pass
    try:
        return float(p_text) <= cutoff
    except ValueError:
        return False


def extract_one(source: Path, dest: Path, file_id: str, cutoff: float) -> dict:
    """
    Write the associations at or below `cutoff` from one formatted file.

    Read in chunks so memory follows the chunk, not the file: the largest of these
    holds ~200 million rows.

    Args:
        source: a formatted .tsv.gz
        dest: parquet to write
        file_id: stamped onto every row, so files can be concatenated
        cutoff: the P value threshold

    Returns:
        counts for the run record
    """
    import pandas as pd

    kept: list[pd.DataFrame] = []
    columns: list[str] = []
    rows_in = 0
    with gzip.open(source, 'rt') as handle:
        for chunk in pd.read_csv(handle, sep='\t', dtype=str, chunksize=CHUNK_ROWS):
            rows_in += len(chunk)
            columns = list(chunk.columns)
            neg_log = chunk.get('neg_log_10_p_value', pd.Series('', index=chunk.index))
            keep = [
                passes(p, n, cutoff)
                for p, n in zip(chunk['p_value'].fillna(''), neg_log.fillna(''))
            ]
            if any(keep):
                kept.append(chunk[keep])

    # A file with no hits still gets a parquet, with the right columns and no rows,
    # so that reading the whole folder back never has to special-case it.
    hits = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=columns)
    hits.insert(0, 'file_id', file_id)
    hits.to_parquet(dest, index=False)
    return {'file_id': file_id, 'rows_in': rows_in, 'rows_kept': len(hits)}


def chromosomes_in(path: str) -> list[str]:
    """
    Chromosomes a formatted file contains, read from its tabix index.

    The index lists its own contigs, so this costs one small read rather than a pass
    over the data. Chromosomes are GWAS-SSF numeric, so X is '23'.

    Args:
        path: a formatted .tsv.gz (its .tbi must sit beside it)
    """
    import pysam

    with pysam.TabixFile(path) as tabix:
        return list(tabix.contigs)


def run_coverage(rows: list[dict], out: str, formatted: str) -> None:
    """
    Write one row per file saying which chromosomes it has, and whether X is one.

    Args:
        rows: files.csv rows
        out: destination CSV
        formatted: folder holding the formatted files
    """
    from cpg_utils import to_path

    records = []
    for row in rows:
        path = f'{formatted}/{formatted_name(row["file_id"])}'
        if not to_path(path).exists():
            print(f'{row["file_id"]}: not formatted, skipped')
            continue
        contigs = chromosomes_in(path)
        records.append(
            {
                'file_id': row['file_id'],
                'biomarker': row['biomarker'],
                'ancestry': row['ancestry'],
                'has_x': '23' in contigs,
                'n_chromosomes': len(contigs),
                'chromosomes': ' '.join(contigs),
            }
        )

    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
    with to_path(out).open('w') as handle:
        handle.write(text.getvalue())

    with_x = sum(record['has_x'] for record in records)
    print(
        f'{len(records)} files: {with_x} with chromosome X, {len(records) - with_x} without'
    )
    print(f'written to {out}')


def run_local(rows: list[dict], out: Path, formatted: str, cutoff: float) -> None:
    """Extract from local or remote files one after another, into a local folder."""
    from cpg_utils import to_path

    out.mkdir(parents=True, exist_ok=True)
    for row in rows:
        dest = out / f'{row["file_id"]}.parquet'
        if dest.exists():
            print(f'{row["file_id"]}: exists, skipped')
            continue
        source = f'{formatted}/{formatted_name(row["file_id"])}'
        local = out / formatted_name(row['file_id'])
        if not local.exists():
            to_path(source).copy(local)
        counts = extract_one(local, dest, row['file_id'], cutoff)
        print(
            f'{counts["file_id"]}: {counts["rows_kept"]:,} of {counts["rows_in"]:,} rows'
        )


def run_batch(rows: list[dict], args) -> None:
    """
    One Hail Batch job per formatted file not yet extracted.

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
        source = f'{args.formatted}/{formatted_name(row["file_id"])}'
        dest = f'{args.out}/{row["file_id"]}.parquet'
        if to_path(dest).exists() and not args.force:
            print(f'{row["file_id"]}: exists, skipped')
            continue
        if not to_path(source).exists():
            missing.append(row['file_id'])
            continue

        size_gib = int(row['size_bytes'] or 0) / 2**30
        job = batch.new_bash_job(f'extract {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2)
        job.storage(f'{int(size_gib * 3) + 10}Gi')
        job.command(
            f'python3 -m pip install --quiet {PIP_PACKAGES}\n'
            f"cat > extract.py <<'GWAS_SUMSTATS_SCRIPT'\n{script}\nGWAS_SUMSTATS_SCRIPT\n"
            f'python3 extract.py --one {shlex.quote(json.dumps(row))} '
            f'--source {batch.read_input(source)} --dest {job.out} '
            f'--threshold {args.threshold}'
        )
        batch.write_output(job.out, dest)
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
    parser.add_argument('--files', type=Path, default=FILES_CSV)
    parser.add_argument('--formatted', default=FORMATTED, help='formatted files folder')
    parser.add_argument('--out', default=DEFAULT_OUT, help='destination folder')
    parser.add_argument('--coverage-out', default=DEFAULT_COVERAGE)
    parser.add_argument('--threshold', type=float, default=THRESHOLD)
    parser.add_argument('--only', nargs='+', help='file_ids to extract')
    parser.add_argument('--local', action='store_true', help='run here, not on Batch')
    parser.add_argument(
        '--coverage', action='store_true', help='write the chromosome coverage table'
    )
    parser.add_argument(
        '--force', action='store_true', help='redo files already written'
    )
    parser.add_argument('--one', help='a single files.csv row as JSON (job mode)')
    parser.add_argument('--source', type=Path, help='job mode: input file')
    parser.add_argument('--dest', type=Path, help='job mode: output parquet')
    args = parser.parse_args()

    if args.one:
        row = json.loads(args.one)
        counts = extract_one(args.source, args.dest, row['file_id'], args.threshold)
        print(
            f'{counts["file_id"]}: {counts["rows_kept"]:,} of {counts["rows_in"]:,} rows'
        )
        return

    rows = read_files_csv(args.files, args.only)
    if args.coverage:
        run_coverage(rows, args.coverage_out, args.formatted)
    elif args.local:
        if args.out.startswith('gs://'):
            parser.error(f'--local needs a local --out folder, got {args.out}')
        run_local(rows, Path(args.out), args.formatted, args.threshold)
    else:
        run_batch(rows, args)


if __name__ == '__main__':
    main()
