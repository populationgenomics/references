#!/usr/bin/env python3

"""
Extract the significant associations from the formatted GWAS summary statistics.

Writes, under `v1/`:

    hits/p_<threshold>/<file_id>.parquet
                                associations with P <= the threshold (default 1e-5),
                                every column kept, typed
    chromosomes_reported.csv    which chromosomes each file reports, and whether the
                                source published chromosome X

Why this exists, how P values and X are handled, and the column types are in
README.md, under Extracted hits.

Both commands write to main, so run them from `main` after merge. Once every
format.py job has finished, the hits, one Hail Batch job per file; the
analysis-runner job only submits them and prints the batch link to follow:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Extract significant GWAS associations" \
        python3 gwas_sumstats/extract.py

Then, once every extract job has finished, the chromosome table, which runs in the
analysis-runner job itself:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "GWAS chromosomes reported" \
        python3 gwas_sumstats/extract.py --chromosomes

On a laptop, one file at a time, with no Batch, gs:// or cpg_utils (see the README
for setup); this is the command each Batch job runs:

    python3 gwas_sumstats/extract.py \
        --file-id 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951 \
        --source ~/gwas_sumstats_local/v1/2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951_GRCh38_formatted.tsv.gz \
        --dest ~/gwas_sumstats_local/hits/p_1e-05/2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951.parquet
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import shlex
from collections import Counter
from pathlib import Path
from statistics import NormalDist
from typing import TYPE_CHECKING

# pandas and pyarrow are imported inside the functions that use them, and cpg_utils
# inside the ones that reach gs://: each Batch job installs only PIP_PACKAGES, and
# the tests and the single-file command run without cpg_utils.
if TYPE_CHECKING:
    import pandas as pd
    import pyarrow as pa

FILES_CSV = Path(__file__).with_name('files.csv')
FORMATTED = 'gs://cpg-common-main/references/gwas_sumstats/v1'
CHROMOSOMES_NAME = 'chromosomes_reported.csv'

THRESHOLD = 1e-5
"""
Deliberately looser than genome-wide significance (5e-8).

Whoever uses these files chooses their own threshold; this one only has to be loose
enough that they never need to come back and re-read the originals. Under the null it
keeps 200 times more rows than 5e-8 while still discarding over 99% of each file.
"""

HIT_TYPES = {
    'file_id': 'string',
    'chromosome': 'int64',
    'base_pair_location': 'int64',
    'effect_allele': 'string',
    'other_allele': 'string',
    'beta': 'float64',
    'standard_error': 'float64',
    'effect_allele_frequency': 'float64',
    'p_value': 'string',
    'neg_log_10_p_value': 'float64',
    'rsid': 'string',
    'n': 'int64',
    'z': 'float64',
    'effect_allele_is_alt': 'bool',
}
"""
Arrow type of each hits column: file_id, then format.py's OUTPUT_COLUMNS in order.

Kept here rather than imported from format.py because each Batch job runs this
script on its own; a test checks the two agree.
"""

CHUNK_ROWS = 2_000_000
PIP_PACKAGES = 'pandas==2.2.3 pyarrow==18.1.0'


def hits_folder(out: str, threshold: float) -> str:
    """
    Return the folder for hits at this threshold: <out>/p_<threshold>.

    Every threshold gets its own subfolder, so hits made at two thresholds can never
    share a folder and be read back as one table. Named from P rather than -log10 P,
    which gives names such as 7.30103 for non-round cutoffs.

    Args:
        out: the parent folder, local or gs://
        threshold: the P value threshold

    Returns:
        e.g. gs://.../v1/hits/p_1e-05 for 1e-5
    """
    return f'{out.rstrip("/")}/p_{threshold:g}'


def formatted_name(file_id: str) -> str:
    """
    Return the name format.py gave this file.

    The manifest's own `formatted_file` column cannot be used: every row of it reads
    'out.tsv.gz', the name each Batch job used before its output was copied out.

    Args:
        file_id: the files.csv file_id

    Returns:
        <file_id>_GRCh38_formatted.tsv.gz
    """
    return f'{file_id}_GRCh38_formatted.tsv.gz'


def passing(
    p_text: pd.Series,
    neg_log_10: pd.Series | None,
    cutoff: float,
    z: pd.Series | None = None,
    beta: pd.Series | None = None,
    se: pd.Series | None = None,
) -> pd.Series:
    """
    Flag the associations at or below the cutoff.

    Uses neg_log_10_p_value where readable, then p_value, then, for rows with
    neither (a published P of 0), |z| or |beta / standard_error| against the
    two-sided z for the cutoff. A standard error of 0 or less is not used.

    Args:
        p_text: the p_value column, text
        neg_log_10: the neg_log_10_p_value column, text; None if absent
        cutoff: the P value threshold
        z: the z column, text; None if absent
        beta: the beta column, text; None if absent
        se: the standard_error column, text; None if absent

    Returns:
        one boolean per row, True to keep
    """
    import pandas as pd

    def number(column: pd.Series | None) -> pd.Series:
        if column is None:
            return pd.Series(float('nan'), index=p_text.index)
        return pd.to_numeric(column, errors='coerce')

    neg_log = number(neg_log_10)
    p_value = number(p_text)
    z_score = number(z)
    standard_error = number(se)
    z_score = z_score.where(
        z_score.notna(), number(beta) / standard_error.where(standard_error > 0)
    )
    # From the lower tail, which stays exact for cutoffs far below 1e-16, where
    # 1 - cutoff / 2 rounds to 1.
    by_z = z_score.abs() >= -NormalDist().inv_cdf(cutoff / 2)
    by_p_value = (p_value <= cutoff).where(p_value.notna(), by_z)
    # format.py rounds neg_log_10_p_value to 6 decimals, so the cutoff is rounded
    # the same way, or a P exactly at a cutoff such as 2e-5 would be dropped.
    by_neg_log = neg_log >= round(-math.log10(cutoff), 6)
    return by_neg_log.where(neg_log.notna(), by_p_value).fillna(False).astype(bool)


def fell_back_on_zero_se(chunk: pd.DataFrame) -> pd.Series:
    """
    Flag rows with no P value and no z whose standard error is 0 or less.

    passing() would have ranked these on beta / standard_error, so they are dropped
    for a source error rather than for being weak, and are counted apart.

    Args:
        chunk: formatted rows, every column text

    Returns:
        one boolean per row
    """
    import pandas as pd

    def missing(column: str) -> pd.Series:
        if column not in chunk:
            return pd.Series(True, index=chunk.index)
        return pd.to_numeric(chunk[column], errors='coerce').isna()

    if 'standard_error' not in chunk:
        return pd.Series(False, index=chunk.index)
    standard_error = pd.to_numeric(chunk['standard_error'], errors='coerce')
    return (
        missing('neg_log_10_p_value')
        & missing('p_value')
        & missing('z')
        & (standard_error <= 0)
    )


def typed(hits: pd.DataFrame) -> pd.DataFrame:
    """
    Convert hits read as text to the types in HIT_TYPES.

    Strict: a value that does not convert raises rather than becoming null, since
    format.py writes only numbers, true/false or NA in these columns.

    Args:
        hits: formatted rows, every column text, NA as missing

    Returns:
        the same rows with typed columns

    Raises:
        ValueError: if a value does not convert to its column's type
    """
    import pandas as pd

    hits = hits.copy()
    for column, kind in HIT_TYPES.items():
        if column not in hits or kind == 'string':
            continue
        if kind == 'bool':
            values = hits[column].map({'true': True, 'false': False})
            unknown = hits[column].notna() & values.isna()
            if unknown.any():
                raise ValueError(
                    f'{column}: not true/false: {hits.loc[unknown, column].iloc[0]!r}'
                )
            hits[column] = values.astype('boolean')
        else:
            numbers = pd.to_numeric(hits[column])
            hits[column] = numbers.astype('Int64' if kind == 'int64' else 'float64')
    return hits


def hits_schema(columns: list[str]) -> pa.Schema:
    """
    Return the schema every hits parquet is written with.

    One explicit schema, rather than types inferred per file, means files with and
    without hits read back as one table.

    Args:
        columns: the formatted file's header

    Returns:
        file_id, then each column, typed as in HIT_TYPES (text if not listed)
    """
    import pyarrow as pa

    return pa.schema(
        [
            (column, pa.type_for_alias(HIT_TYPES.get(column, 'string')))
            for column in ['file_id', *columns]
        ]
    )


def hits_in(chunk: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    """
    Return the rows of a chunk at or below the cutoff.

    Args:
        chunk: formatted rows, every column text
        cutoff: the P value threshold

    Returns:
        the passing rows, still text
    """
    return chunk[
        passing(
            chunk['p_value'],
            chunk.get('neg_log_10_p_value'),
            cutoff,
            z=chunk.get('z'),
            beta=chunk.get('beta'),
            se=chunk.get('standard_error'),
        )
    ]


def extract_one(source: Path, dest: Path, file_id: str, cutoff: float) -> dict:
    """
    Write the associations at or below `cutoff` from one formatted file.

    Reads in chunks so memory follows the chunk, not the file: the largest holds
    ~200 million rows. Every row is counted per chromosome for the chromosome table.
    The threshold and counts go into the parquet's file metadata under
    b'gwas_sumstats', so each file records how it was made.

    Args:
        source: a formatted .tsv.gz
        dest: parquet to write
        file_id: stamped onto every row, so files can be concatenated
        cutoff: the P value threshold

    Returns:
        file_id, rows_in, rows_kept, rows_kept_without_p, rows_dropped_se_zero and
        rows_per_chromosome

    Raises:
        ValueError: if a hit's value does not fit its column's type
    """
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    # Each filtered chunk is written as it comes, so memory follows the chunk rather
    # than the number of hits. The header is read first so the schema is known even
    # for a file with no rows.
    rows_in = rows_kept = rows_kept_without_p = rows_dropped_se_zero = 0
    rows_per_chromosome: Counter[str] = Counter()
    with gzip.open(source, 'rt') as handle:
        columns = handle.readline().rstrip('\n').split('\t')
        schema = hits_schema(columns)
        with pq.ParquetWriter(dest, schema) as writer:
            for chunk in pd.read_csv(
                handle,
                sep='\t',
                dtype=str,
                header=None,
                names=columns,
                chunksize=CHUNK_ROWS,
            ):
                rows_in += len(chunk)
                rows_per_chromosome.update(chunk['chromosome'].value_counts().to_dict())
                rows_dropped_se_zero += int(fell_back_on_zero_se(chunk).sum())
                hits = hits_in(chunk, cutoff)
                if hits.empty:
                    continue
                hits.insert(0, 'file_id', file_id)
                rows_kept += len(hits)
                rows_kept_without_p += int(hits['p_value'].isna().sum())
                writer.write_table(
                    pa.Table.from_pandas(
                        typed(hits), schema=schema, preserve_index=False
                    )
                )
            counts = {
                'file_id': file_id,
                'rows_in': rows_in,
                'rows_kept': rows_kept,
                'rows_kept_without_p': rows_kept_without_p,
                'rows_dropped_se_zero': rows_dropped_se_zero,
                'rows_per_chromosome': dict(
                    sorted(rows_per_chromosome.items(), key=lambda item: int(item[0]))
                ),
            }
            record = json.dumps({**counts, 'threshold': cutoff})
            writer.add_key_value_metadata({'gwas_sumstats': record})
    return counts


def chromosomes_record(row: dict, record: dict | None, stats: dict | None) -> dict:
    """
    Return one row of the chromosome table.

    has_x needs more rows on X than liftover moved between chromosomes in total,
    the most that could have landed on X from elsewhere. A GRCh38 source is not
    lifted, so its count is 0 and any X row is its own.

    Args:
        row: the files.csv row
        record: the hits parquet's gwas_sumstats metadata; None if not extracted
        stats: format.py's stats for the file; None if not extracted

    Returns:
        file details, extracted, has_x, x_rows, liftover_chromosome_changed and
        chromosomes; the last four empty when not extracted
    """
    details = {
        'file_id': row['file_id'],
        'biomarker': row['biomarker'],
        'ancestry': row['ancestry'],
        'source_build': row['source_build'],
        'extracted': record is not None,
    }
    if record is None:
        empty = ['has_x', 'x_rows', 'liftover_chromosome_changed', 'chromosomes']
        return {**details, **dict.fromkeys(empty, '')}
    per_chromosome = record['rows_per_chromosome']
    x_rows = per_chromosome.get('23', 0)
    changed = stats.get('liftover_chromosome_changed', 0)
    return {
        **details,
        'has_x': x_rows > changed,
        'x_rows': x_rows,
        'liftover_chromosome_changed': changed,
        'chromosomes': ' '.join(per_chromosome),
    }


def run_chromosomes(rows: list[dict], out: str, hits: str, stats: str) -> None:
    """
    Write one row per file saying which chromosomes it reports, and whether the
    source published X.

    Files not extracted yet are listed with extracted=False, so a partial run is
    visible in the table itself.

    Args:
        rows: files.csv rows
        out: destination CSV
        hits: folder of hits parquets, any threshold: every row was counted
        stats: format.py's per-file stats folder

    Raises:
        ValueError: if none of the files has been extracted, or one has no stats
    """
    import pyarrow.parquet as pq
    from format import as_path

    records = []
    for row in rows:
        parquet = as_path(f'{hits}/{row["file_id"]}.parquet')
        record = file_stats = None
        if parquet.exists():
            with parquet.open('rb') as handle:
                record = json.loads(pq.read_metadata(handle).metadata[b'gwas_sumstats'])
            stats_file = as_path(f'{stats}/{row["file_id"]}.json')
            if not stats_file.exists():
                raise ValueError(f'{row["file_id"]}: no format.py stats at {stats_file}')
            file_stats = json.loads(stats_file.read_text())
        else:
            print(f'{row["file_id"]}: not extracted')
        records.append(chromosomes_record(row, record, file_stats))
    done = [record for record in records if record['extracted']]
    if not done:
        raise ValueError(f'none of the {len(rows)} selected files is extracted in {hits}')

    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
    with as_path(out).open('w') as handle:
        handle.write(text.getvalue())

    with_x = sum(record['has_x'] for record in done)
    print(
        f'{len(done)} files: {with_x} with chromosome X, {len(done) - with_x} without; '
        f'{len(records) - len(done)} not extracted'
    )
    print(f'written to {out}')


def is_current(dest: Path, source: Path) -> bool:
    """
    Return whether hits already written are newer than their formatted file.

    After format.py --only <id> --force, the formatted file is newer than its hits,
    so the hits are redone rather than skipped.

    Args:
        dest: the hits parquet, local or a cloud path
        source: the formatted file, local or a cloud path

    Returns:
        True if dest exists and is at least as new as source
    """
    return dest.exists() and dest.stat().st_mtime >= source.stat().st_mtime


def job_storage(source_bytes: int, file_id: str) -> str:
    """
    Return the disk to request for a job that reads this formatted file.

    Sized from the file the job actually localises, with room for the parquet it
    writes. Refuses a size of 0 or less, which would mean the file is broken,
    rather than falling back to the floor.

    Args:
        source_bytes: size of the formatted .tsv.gz
        file_id: for the error message

    Returns:
        a Hail Batch storage string, e.g. '13Gi'

    Raises:
        ValueError: if the size is 0 or less
    """
    if source_bytes <= 0:
        raise ValueError(f'{file_id}: formatted file has size {source_bytes}')
    return f'{int(source_bytes / 2**30 * 3) + 10}Gi'


def run_batch(rows: list[dict], args: argparse.Namespace, out: str) -> None:
    """
    Submit one Hail Batch job per formatted file not yet extracted, or re-formatted
    since.

    Each job reads one file and writes one parquet, so 336 files take as long as the
    largest single file rather than the sum of all of them.

    Args:
        rows: files.csv rows
        args: the parsed command line
        out: folder for this threshold's hits, from hits_folder
    """
    from cpg_utils import to_path
    from cpg_utils.config import config_retrieve
    from cpg_utils.hail_batch import get_batch

    batch = get_batch(name='Extract significant GWAS associations')
    script = Path(__file__).read_text()
    submitted, missing = 0, []
    for row in rows:
        source = to_path(f'{args.formatted}/{formatted_name(row["file_id"])}')
        dest = to_path(f'{out}/{row["file_id"]}.parquet')
        if not source.exists():
            missing.append(row['file_id'])
            continue
        if not args.force and is_current(dest, source):
            print(f'{row["file_id"]}: exists, skipped')
            continue

        job = batch.new_bash_job(f'extract {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2)
        job.storage(job_storage(source.stat().st_size, row['file_id']))
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


def build_parser() -> argparse.ArgumentParser:
    """
    Return the command-line parser.

    Returns:
        the parser for every mode: batch extraction, the chromosome table and one file
    """
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument(
        '--files', type=Path, default=FILES_CSV, help='files.csv to read'
    )
    parser.add_argument(
        '--formatted',
        default=FORMATTED,
        help='formatted files folder; the other paths default to places in it',
    )
    parser.add_argument(
        '--out',
        help='hits parent folder (default hits/ in --formatted); hits go to a '
        'subfolder per threshold, e.g. p_1e-05',
    )
    parser.add_argument(
        '--stats', help="format.py's per-file stats (default stats/ in --formatted)"
    )
    parser.add_argument(
        '--chromosomes-out',
        help=f'chromosome table to write (default {CHROMOSOMES_NAME} in --formatted); '
        'required with --only or another --files',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=THRESHOLD,
        help=f'P value (not -log10 P) to keep at or below (default {THRESHOLD})',
    )
    parser.add_argument('--only', nargs='+', help='file_ids to extract')
    parser.add_argument(
        '--chromosomes',
        action='store_true',
        help='write the table of chromosomes each file reports; run after extracting',
    )
    parser.add_argument(
        '--force', action='store_true', help='redo files already written; needs --only'
    )
    parser.add_argument(
        '--file-id', help='one file (a Batch job, or a laptop): file_id for every row'
    )
    parser.add_argument('--source', type=Path, help='one file: formatted .tsv.gz')
    parser.add_argument('--dest', type=Path, help='one file: parquet to write')
    return parser


def check_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """
    Refuse argument combinations that would waste a run or overwrite shared output.

    Args:
        parser: reports the error and exits
        args: the parsed command line
    """
    if not 0 < args.threshold < 1:
        parser.error(f'--threshold is a P value, so between 0 and 1; got {args.threshold}')
    if args.file_id:
        return
    if args.force and not args.only:
        parser.error('--force needs --only, so every file is not redone by mistake')
    subset = args.only or args.files.resolve() != FILES_CSV.resolve()
    if args.chromosomes and subset and not args.chromosomes_out:
        parser.error(
            '--chromosomes on a subset of files (--only or another --files) would '
            f'replace {CHROMOSOMES_NAME} with a partial table; give an output path: '
            '--chromosomes-out subset.csv'
        )


def run_one_file(args: argparse.Namespace) -> None:
    """
    Extract one file and print its counts: what each Batch job runs.

    Args:
        args: the parsed command line, with file_id, source and dest
    """
    counts = extract_one(args.source, args.dest, args.file_id, args.threshold)
    print(
        f'{counts["file_id"]}: {counts["rows_kept"]:,} of {counts["rows_in"]:,} '
        f'rows, {counts["rows_kept_without_p"]:,} of them kept on z; '
        f'{counts["rows_dropped_se_zero"]:,} with no P value dropped for a '
        'standard error of 0'
    )


def main() -> None:
    """Run one file, the chromosome table, or the batch extraction."""
    parser = build_parser()
    args = parser.parse_args()
    check_args(parser, args)
    if args.file_id:
        run_one_file(args)
        return

    from download import read_files_csv

    formatted = args.formatted.rstrip('/')
    hits = hits_folder(args.out or f'{formatted}/hits', args.threshold)
    rows = read_files_csv(args.files, args.only)
    if args.chromosomes:
        try:
            run_chromosomes(
                rows,
                args.chromosomes_out or f'{formatted}/{CHROMOSOMES_NAME}',
                hits,
                args.stats or f'{formatted}/stats',
            )
        except ValueError as error:
            parser.error(str(error))
    else:
        run_batch(rows, args, hits)


if __name__ == '__main__':
    main()
