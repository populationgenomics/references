#!/usr/bin/env python3

"""
Download the GWAS summary statistics listed in files.csv, unchanged, and check each
one against its published MD5 where the source publishes one.

Each file lands as <file_id>_<original|harmonised><suffix>, next to a <same name>.json
record of where it came from (URL, MD5, size, download time). The name leaves out the
genome build: that is files.csv's claim, which format.py checks and may prove wrong,
and correcting it must not orphan the download. A file whose destination already
exists is skipped, so a rerun resumes where the last stopped.

On Hail Batch, one job per file, writing to the tmp bucket (objects there are deleted
after 8 days, so run format.py within the week):

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Download biomarker GWAS summary statistics" \
        python3 gwas_sumstats/download.py

Locally, sequentially, into a folder (for testing on a few files):

    python3 gwas_sumstats/download.py --local --out ./original \
        --only 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951
"""

import argparse
import csv
import hashlib
import http.client
import json
import shlex
import shutil
import time
import urllib.request
import zlib
from datetime import datetime, timezone
from pathlib import Path

FILES_CSV = Path(__file__).with_name('files.csv')
DEFAULT_OUT = 'gs://cpg-common-main-tmp/gwas_sumstats/original'
USER_AGENT = 'Mozilla/5.0 (compatible; cpg-references-gwas-sumstats)'
CHUNK = 1 << 20
RETRIES = 4


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


def original_name(row: dict) -> str:
    """File name of the downloaded original, e.g. ..._original.tsv.gz"""
    kind = 'harmonised' if row['source_kind'] == 'catalog_harmonised' else 'original'
    return f'{row["file_id"]}_{kind}{row["suffix"]}'


def open_url(url: str, headers: dict | None = None):
    """Open a URL with a user agent; some hosts refuse urllib's default one."""
    headers = {'User-Agent': USER_AGENT} | (headers or {})
    request = urllib.request.Request(url, headers=headers)
    return urllib.request.urlopen(request, timeout=300)  # noqa: S310


def fetch_head(url: str, n_bytes: int = 200_000) -> str:
    """
    The first complete lines of a (possibly gzip/bgzip) file, from its first
    n_bytes. A server that ignores the Range header is read only that far.
    """
    with open_url(url, {'Range': f'bytes=0-{n_bytes - 1}'}) as response:
        data = response.read(n_bytes)
    text = b''
    while data[:2] == b'\x1f\x8b':
        member = zlib.decompressobj(31)
        text += member.decompress(data)
        data = member.unused_data
    lines = (text or data).decode('utf-8', 'replace').splitlines(keepends=True)
    return ''.join(lines[:-1])


def published_md5(row: dict) -> str | None:
    """
    The MD5 the source publishes for this file, if any.

    files.csv carries it directly (md5) or points at a GWAS Catalog md5sum.txt
    (md5_url), whose lines are '<md5> <file name>' with one or two spaces.

    Args:
        row: one files.csv row
    """
    if row['md5']:
        return row['md5']
    if not row['md5_url']:
        return None
    file_name = row['source_url'].rsplit('/', 1)[1]
    with open_url(row['md5_url']) as response:
        for line in response.read().decode().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1] == file_name:
                return parts[0]
    raise ValueError(f'{file_name} is not listed in {row["md5_url"]}')


def fetch(url: str, dest: Path, expected_md5: str | None = None) -> tuple[str, int]:
    """
    Stream a URL to a local file, retrying with backoff until it arrives whole.

    A body shorter than its Content-Length counts as a failed attempt: urllib
    returns it without raising (http.client only raises IncompleteRead for chunked
    responses). So does an MD5 that differs from expected_md5. After the last
    failed attempt the partial file is deleted and the error raised.

    Args:
        url: source URL
        dest: local file to write
        expected_md5: the publisher's MD5, if any

    Returns:
        (MD5 hex digest, size in bytes) of what was written
    """
    for attempt in range(1, RETRIES + 1):
        try:
            md5 = hashlib.md5()  # noqa: S324
            size = 0
            with open_url(url) as response, dest.open('wb') as out:
                length = response.headers.get('Content-Length')
                while chunk := response.read(CHUNK):
                    md5.update(chunk)
                    size += len(chunk)
                    out.write(chunk)
            if length is not None and size != int(length):
                raise OSError(f'got {size:,} of {int(length):,} bytes')
            if expected_md5 and md5.hexdigest() != expected_md5:
                raise OSError(f'MD5 {md5.hexdigest()} != published {expected_md5}')
            return md5.hexdigest(), size
        except (OSError, http.client.HTTPException) as error:
            if attempt == RETRIES:
                dest.unlink(missing_ok=True)
                if isinstance(error, http.client.HTTPException):
                    raise OSError(f'{url}: {error!r}') from error
                raise
            wait = 30 * attempt
            print(f'{url}: {error!r}; retry {attempt}/{RETRIES - 1} in {wait}s')
            time.sleep(wait)
    raise AssertionError('unreachable')


def download_one(row: dict, data_path: Path, record_path: Path) -> None:
    """
    Download one file, refuse it if its MD5 differs from the published one, and
    write the provenance record next to it.

    Args:
        row: one files.csv row
        data_path: local path for the downloaded file
        record_path: local path for the JSON record
    """
    expected = published_md5(row)
    md5, size = fetch(row['source_url'], data_path, expected)
    record = {
        'file_id': row['file_id'],
        'file_name': original_name(row),
        'source_url': row['source_url'],
        'md5': md5,
        'md5_published': expected or '',
        'md5_checked': bool(expected),
        'size_bytes': size,
        'downloaded_at': datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }
    record_path.write_text(json.dumps(record, indent=1) + '\n')
    print(f'{row["file_id"]}: {size:,} bytes, MD5 checked: {bool(expected)}')


def run_local(rows: list[dict], out: Path) -> None:
    """Download rows one after another into a local folder."""
    out.mkdir(parents=True, exist_ok=True)
    for row in rows:
        data_path = out / original_name(row)
        if data_path.exists():
            print(f'{row["file_id"]}: exists, skipped')
            continue
        partial = data_path.with_name(data_path.name + '.partial')
        download_one(row, partial, out / f'{data_path.name}.json')
        shutil.move(partial, data_path)


def run_batch(rows: list[dict], out: str) -> None:
    """
    One Hail Batch job per file not yet in `out`. The job runs this script on
    one row (passed as JSON) and Batch copies its two outputs to `out` only when
    the job succeeds, so a failed MD5 check leaves nothing behind.
    """
    from cpg_utils import to_path
    from cpg_utils.config import config_retrieve
    from cpg_utils.hail_batch import get_batch

    batch = get_batch(name='Download biomarker GWAS summary statistics')
    script = Path(__file__).read_text()
    submitted = 0
    for row in rows:
        dest = f'{out}/{original_name(row)}'
        if to_path(dest).exists():
            print(f'{row["file_id"]}: exists, skipped')
            continue
        job = batch.new_bash_job(f'download {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        size_gib = int(row['size_bytes'] or 0) / 2**30
        job.storage(f'{max(10, int(size_gib * 1.5) + 5)}Gi')
        job.command(
            "cat > download.py <<'GWAS_SUMSTATS_SCRIPT'\n"
            f'{script}\n'
            'GWAS_SUMSTATS_SCRIPT\n'
            f'python3 download.py --one {shlex.quote(json.dumps(row))} '
            f'--data {job.data} --record {job.record}'
        )
        batch.write_output(job.data, dest)
        batch.write_output(job.record, f'{dest}.json')
        submitted += 1
    print(f'{submitted} download jobs submitted')
    if submitted:
        batch.run(wait=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--files', type=Path, default=FILES_CSV)
    parser.add_argument('--out', default=DEFAULT_OUT, help='destination folder')
    parser.add_argument('--only', nargs='+', help='file_ids to download')
    parser.add_argument('--local', action='store_true', help='run here, not on Batch')
    parser.add_argument('--one', help='a single files.csv row as JSON (job mode)')
    parser.add_argument('--data', type=Path, help='job mode: output data path')
    parser.add_argument('--record', type=Path, help='job mode: output JSON path')
    args = parser.parse_args()

    if args.one:
        download_one(json.loads(args.one), args.data, args.record)
        return
    rows = read_files_csv(args.files, args.only)
    if args.local:
        if args.out.startswith('gs://'):
            parser.error(f'--local needs a local --out folder, got {args.out}')
        run_local(rows, Path(args.out))
    else:
        run_batch(rows, args.out.rstrip('/'))


if __name__ == '__main__':
    main()
