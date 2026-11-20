#!/usr/bin/env python3

"""
The steps of alphagenome_echtvar.py that need no cloud library: reading one of
DeepMind's AlphaGenome score tables by byte range, merging per-chromosome echtvar
archives into one, and writing the receipt of a finished archive.

alphagenome_echtvar.py imports this module to plan a run, and copies the file into
every job. There it runs in the annotation pipeline's tools image, which has Python
and no cloud SDK, so the module uses the standard library alone:

    python3 alphagenome_echtvar_steps.py rows --plan plan.json --chromosome chr1
    python3 alphagenome_echtvar_steps.py sample --plan plan.json
    python3 alphagenome_echtvar_steps.py merge --out merged.zip chr1.zip chr2.zip
    python3 alphagenome_echtvar_steps.py records --plan plan.json ...

Reading by byte range. DeepMind publishes each table as a zip of two members, both
stored without compression: the table, bgzipped and sorted by position, and its tabix
index. A stored member is a run of the zip's own bytes, and the index gives the place
of each chromosome's rows in the table. So a chromosome's rows are one byte range of
the zip, and no job copies or unpacks a zip, the largest of which has 284 GB.
`locate` finds the two members, reads the index and returns the plan: the object and
its generation, and for each chromosome its range, its number of rows and the blocks
that the check reads. `rows` writes one chromosome's rows as text and fails unless it
wrote as many as the index counts. `sample` writes the rows of the check's blocks.

A source is an object in Cloud Storage, gs://<bucket>/<name>, or a local file in a
test. An object is read through the JSON API with a token from the metadata server,
which Hail Batch runs for every job and which answers as the job's service account.
Every read asks for the generation that `locate` found, so all jobs of a run read
the same bytes, and a run fails if the object is replaced under it.

Merging. An echtvar archive is a zip. `echtvar encode` (0.2.4) writes a directory
entry `echtvar/`, then `echtvar/config.json`, then for every 1 Mb block of every
chromosome one `<alias>.bin` per field, `var32.bin` and `too-long-for-var32.enc`
under `echtvar/<chromosome>/<block>/`. `echtvar anno` opens entries by name. A block
belongs to one chromosome, so an archive encoded from one chromosome holds exactly
the entries that a whole-genome encode writes for that chromosome. `merge` copies the
entries of per-chromosome archives into one archive, as bytes and without inflating
them, and writes a new directory of the entries. It refuses inputs whose
`config.json` differ. On sample rows of the three tables, the merged archive was the
same bytes as the archive of one encode over the same BCFs.

The receipt. `records` writes what is published beside a checked archive: the
configuration that echtvar embedded in it, receipt.json, and TERMS.txt for a table
that is not for permissive use.
"""

import argparse
import dataclasses
import datetime
import http.client
import importlib.metadata
import io
import itertools
import json
import os
import re
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib
from pathlib import Path
from typing import Any, BinaryIO, Callable, Iterator, Mapping, Protocol, Sequence

API = 'https://storage.googleapis.com/storage/v1'
METADATA = (
    'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/'
    'default'
)
# One request of a streamed range. A failed request is made again whole.
STEP = 16 * 2**20
TRIES = 6
# After these statuses, and after a request without an answer, another try can end
# differently.
TRY_AGAIN = frozenset({408, 429, 500, 502, 503, 504})
TIMEOUT = 120.0

# A BGZF block as bgzip writes it: a gzip member with one extra field, BC, which
# holds the block's size less one.
BLOCK_START = b'\x1f\x8b\x08\x04'
BLOCK_EXTRA = b'\x06\x00BC\x02\x00'
BLOCK_HEADER = 18
MAX_BLOCK = 2**16
# The empty block that ends a bgzipped file.
END_OF_FILE = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
# The bin of a tabix index that holds no rows: its two chunks are the range of the
# sequence's rows, and the numbers of its rows with and without a position.
SUMMARY_BIN = 37450
# The last bytes of a chromosome's rows are searched for its last row.
TAIL = 4096
# The most bytes of a table that can lie before its first row and still be a header.
HEADER_LIMIT = 2**20

# The fixed parts of a zip's records, from the ZIP application note (APPNOTE.TXT).
LOCAL_HEADER = struct.Struct('<4s5H3I2H')
CENTRAL_HEADER = struct.Struct('<4s6H3I5H2I')
END_RECORD = struct.Struct('<4s4H2IH')
END_RECORD_64 = struct.Struct('<4sQ2H2I4Q')
END_LOCATOR_64 = struct.Struct('<4sIQI')
# A size or an offset from this value up, and a number of entries from 0xFFFF up, is
# written to the zip64 records.
ZIP64_BYTES = 0xFFFFFFFF
ZIP64_ENTRIES = 0xFFFF
ZIP64_EXTRA = 0x0001
ZIP64_VERSION = 45
COPY_STEP = 16 * 2**20

ARCHIVE_ROOT = 'echtvar/'
ARCHIVE_CONFIG = 'echtvar/config.json'
# echtvar numbers the values of a string field in the order that one encode meets
# them, so the numbers of two encodes do not agree.
ARCHIVE_STRINGS = 'echtvar/strings/'
ARCHIVE_ENTRY = re.compile(r'echtvar/[^/]+/[0-9]+/[^/]+\.(?:bin|enc)')


class UnreadableError(Exception):
    """A source could not be read: it is not there, the account may not read it,
    or the server did not answer."""


class RefusedError(Exception):
    """A request that the server answered with a status that another try does not
    change."""

    def __init__(self, status: int, body: str) -> None:
        super().__init__(f'HTTP {status}: {" ".join(body.split())[:300]}')
        self.status = status


Fetch = Callable[[str, Mapping[str, str]], bytes]


def fetch(
    url: str,
    headers: Mapping[str, str],
    *,
    opener: Callable[..., Any] = urllib.request.urlopen,
    pause: Callable[[float], None] = time.sleep,
) -> bytes:
    """
    The body of one GET. A request that fails in a way that another try can change
    is made again, a few times.

    Args:
        url: one of the two addresses this module builds, the JSON API's or the
            metadata server's
        headers: the request's headers
        opener: stands in for urlopen in a test
        pause: stands in for sleep in a test

    Raises:
        RefusedError: for any other status than success
        UnreadableError: if every try failed
    """
    last = ''
    for attempt in range(TRIES):
        if attempt:
            pause(2.0**attempt)
        try:
            request = urllib.request.Request(url, headers=dict(headers))
            with opener(request, timeout=TIMEOUT) as answer:  # noqa: S310
                return answer.read()
        except urllib.error.HTTPError as error:
            body = error.read().decode('utf-8', 'replace')
            if error.code not in TRY_AGAIN:
                raise RefusedError(error.code, body) from error
            last = f'HTTP {error.code}'
        except (OSError, http.client.HTTPException) as error:
            last = f'{type(error).__name__}: {error}'
    raise UnreadableError(
        f'{url.partition("?")[0]} did not answer in {TRIES} tries, the last of which '
        f'ended with {last}'
    )


@dataclasses.dataclass(frozen=True)
class Described:
    """What was read: an object and the generation of it, or a local file, which
    has none."""

    url: str
    generation: str | None
    size: int
    crc32c: str | None
    updated: str | None


class Source(Protocol):
    """Bytes that can be read by range."""

    @property
    def described(self) -> Described: ...

    def read(self, start: int, end: int) -> bytes:
        """The bytes from `start` up to `end`, which lie inside the source."""


class LocalFile:
    """A file on this machine, for tests."""

    def __init__(self, path: Path) -> None:
        if not path.is_file():
            raise UnreadableError(f'{path} is no file')
        self._path = path
        self._described = Described(str(path), None, path.stat().st_size, None, None)

    @property
    def described(self) -> Described:
        return self._described

    def read(self, start: int, end: int) -> bytes:
        with self._path.open('rb') as file:
            file.seek(start)
            return whole(file.read(end - start), start, end, str(self._path))


class GcsObject:
    """One generation of an object in Cloud Storage, read as the service account
    that the metadata server answers for."""

    def __init__(
        self, bucket: str, name: str, generation: str | None, get: Fetch = fetch
    ) -> None:
        """
        Finds the object.

        Args:
            bucket: the object's bucket
            name: the object's name
            generation: the generation to read, or None for the live one
            get: stands in for `fetch` in a test

        Raises:
            UnreadableError: if the object, or that generation of it, is not there
                or the account may not read it
        """
        self._url = f'gs://{bucket}/{name}'
        quoted = [urllib.parse.quote(part, safe='') for part in (bucket, name)]
        self._object = f'{API}/b/{quoted[0]}/o/{quoted[1]}'
        self._get = get
        self._token = ''
        self._token_until = 0.0
        query = 'fields=generation,size,crc32c,updated'
        if generation is not None:
            query += f'&generation={generation}'
        found = json.loads(self._asked(f'{self._object}?{query}', {}))
        self._described = Described(
            self._url,
            str(found['generation']),
            int(found['size']),
            found.get('crc32c'),
            found.get('updated'),
        )

    @property
    def described(self) -> Described:
        return self._described

    def read(self, start: int, end: int) -> bytes:
        url = f'{self._object}?alt=media&generation={self._described.generation}'
        data = self._asked(url, {'Range': f'bytes={start}-{end - 1}'})
        return whole(data, start, end, self._url)

    def _asked(self, url: str, headers: Mapping[str, str]) -> bytes:
        try:
            authorised = {**headers, 'Authorization': f'Bearer {self._access_token()}'}
            return self._get(url, authorised)
        except RefusedError as refused:
            if refused.status in (401, 403):
                raise UnreadableError(
                    f'{self._url}: {self._account()} may not read it ({refused})'
                ) from refused
            if refused.status == 404:
                raise UnreadableError(
                    f'{self._url}: no such object, or the generation of it that the '
                    f'run began with is gone ({refused})'
                ) from refused
            raise UnreadableError(f'{self._url}: {refused}') from refused

    def _access_token(self) -> str:
        if time.monotonic() >= self._token_until:
            flavor = {'Metadata-Flavor': 'Google'}
            try:
                token = json.loads(self._get(f'{METADATA}/token', flavor))
            except RefusedError as refused:
                raise UnreadableError(
                    f'the metadata server gave no token for the account of this job: '
                    f'{refused}'
                ) from refused
            self._token = token['access_token']
            self._token_until = time.monotonic() + float(token['expires_in']) - 120
        return self._token

    def _account(self) -> str:
        try:
            email = self._get(f'{METADATA}/email', {'Metadata-Flavor': 'Google'})
        except (RefusedError, UnreadableError):
            return 'the account that this job runs as'
        return f'the account that this job runs as, {email.decode().strip()},'


def whole(data: bytes, start: int, end: int, where: str) -> bytes:
    """`data` if it is all the bytes that a read of `start` to `end` asked for."""
    if len(data) != end - start:
        raise UnreadableError(
            f'{where}: a read of the bytes {start} to {end} returned {len(data)} bytes'
        )
    return data


def open_source(url: str, generation: str | None = None, get: Fetch = fetch) -> Source:
    """
    The source at gs://<bucket>/<name>, or the local file at an absolute path.

    Raises:
        ValueError: for any other address, and for a generation of a local file
        UnreadableError: as the source's class does
    """
    if url.startswith('gs://'):
        bucket, _, name = url.removeprefix('gs://').partition('/')
        if not bucket or not name:
            raise ValueError(f'{url!r} is not gs://<bucket>/<name>')
        return GcsObject(bucket, name, generation, get)
    if not url.startswith('/'):
        raise ValueError(
            f'{url!r} is neither gs://<bucket>/<name> nor an absolute path'
        )
    if generation is not None:
        raise ValueError(f'{url} is a local file, which has no generation {generation}')
    return LocalFile(Path(url))


def stream(source: Source, start: int, end: int, step: int = STEP) -> Iterator[bytes]:
    """The bytes from `start` up to `end`, one request of at most `step` bytes at a
    time."""
    for at in range(start, end, step):
        yield source.read(at, min(at + step, end))


class SeekableSource(io.RawIOBase):
    """A source as the seekable file that `zipfile` reads a zip's directory from."""

    def __init__(self, source: Source) -> None:
        super().__init__()
        self._source = source
        self._at = 0

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._at

    def seek(self, offset: int, whence: int = io.SEEK_SET, /) -> int:
        size = self._source.described.size
        at = offset + {io.SEEK_SET: 0, io.SEEK_CUR: self._at, io.SEEK_END: size}[whence]
        if at < 0:
            raise OSError(f'a seek to {at}, before the first byte')
        self._at = at
        return at

    def read(self, size: int = -1, /) -> bytes:
        end = self._source.described.size
        if size >= 0:
            end = min(end, self._at + size)
        data = self._source.read(self._at, end) if end > self._at else b''
        self._at += len(data)
        return data


# ---- the table inside DeepMind's zip ----


@dataclasses.dataclass(frozen=True)
class Table:
    """The table inside a zip: its member, where the member's bytes start in the
    zip, and its index, inflated."""

    member: str
    start: int
    size: int
    index_member: str
    index: bytes


def table_of(source: Source) -> Table:
    """
    The table that a zip holds with its index.

    Raises:
        ValueError: if the zip does not hold exactly one <table>.tsv.gz and its
            <table>.tsv.gz.tbi, both stored without compression
    """
    where = source.described.url
    try:
        with zipfile.ZipFile(SeekableSource(source)) as archive:
            members = {info.filename: info for info in archive.infolist()}
            tables = [name for name in members if name.endswith('.tsv.gz')]
            if len(tables) != 1 or set(members) != {tables[0], f'{tables[0]}.tbi'}:
                held = ', '.join(sorted(members)) or 'no member'
                raise ValueError(
                    f'{where} holds {held}, and not one <table>.tsv.gz with its .tbi'
                )
            for info in members.values():
                stored = info.compress_type == zipfile.ZIP_STORED
                if not stored or info.compress_size != info.file_size:
                    raise ValueError(
                        f'{where}: the member {info.filename} is compressed, so its '
                        f'bytes are not a range of the zip'
                    )
            table = members[tables[0]]
            # Reading the member through the zip also compares it with the checksum
            # that the zip records for it.
            index = inflated(archive.read(f'{table.filename}.tbi'))
    except zipfile.BadZipFile as error:
        raise ValueError(f'{where} is no zip: {error}') from error
    header = source.read(table.header_offset, table.header_offset + LOCAL_HEADER.size)
    signature, *_, name_length, extra_length = LOCAL_HEADER.unpack(header)
    if signature != b'PK\x03\x04':
        raise ValueError(
            f'{where}: no member starts at byte {table.header_offset}, where the zip '
            f'lists the table'
        )
    start = table.header_offset + LOCAL_HEADER.size + name_length + extra_length
    return Table(
        table.filename, start, table.file_size, f'{table.filename}.tbi', index
    )


def inflated(data: bytes) -> bytes:
    """The text of a bgzipped file that is held whole."""
    out = []
    while data:
        member = zlib.decompressobj(wbits=31)
        out.append(member.decompress(data))
        data = member.unused_data
    return b''.join(out)


@dataclasses.dataclass(frozen=True)
class Indexed:
    """
    A chromosome in a tabix index.

    An offset is a virtual one: a block's place in the table times 65,536, plus a
    place in the block's text.
    """

    name: str
    start: int  # of its first row
    end: int  # after its last row
    rows: int
    windows: tuple[int, ...]  # the first row of each 16 kb window, or 0


def indexed(index: bytes) -> list[Indexed]:
    """
    The chromosomes of an inflated tabix index, in the table's order.

    Raises:
        ValueError: if the bytes are no tabix index, or a chromosome has no summary
            of its rows that agrees with its bins
    """
    if index[:4] != b'TBI\x01':
        raise ValueError('the index member is no tabix index')
    sequences, *_, names_length = struct.unpack_from('<8i', index, 4)
    names = index[36 : 36 + names_length].split(b'\x00')[:-1]
    if len(names) != sequences:
        raise ValueError(
            f'the index lists {len(names)} sequences and counts {sequences}'
        )
    at = 36 + names_length
    found = []
    for name in names:
        (bins,) = struct.unpack_from('<i', index, at)
        at += 4
        summary = None
        chunks: list[int] = []
        for _ in range(bins):
            number, count = struct.unpack_from('<Ii', index, at)
            offsets = struct.unpack_from(f'<{2 * count}Q', index, at + 8)
            at += 8 + 16 * count
            if number == SUMMARY_BIN:
                summary = offsets
            else:
                chunks.extend(offsets)
        (count,) = struct.unpack_from('<i', index, at)
        windows = struct.unpack_from(f'<{count}Q', index, at + 4)
        at += 4 + 8 * count
        if (
            summary is None
            or len(summary) != 4
            or not chunks
            or (min(chunks), max(chunks)) != summary[:2]
        ):
            raise ValueError(
                f'the index has no summary of the rows of {name.decode()} that agrees '
                f'with its bins, as an index written by tabix has'
            )
        found.append(
            Indexed(name.decode(), summary[0], summary[1], summary[2], windows)
        )
    return found


def block_bytes(data: bytes | memoryview, at: int, where: str) -> int:
    """
    The size of the BGZF block whose header starts at `at`.

    Raises:
        ValueError: if no block as bgzip writes it starts there
    """
    header = bytes(data[at : at + BLOCK_HEADER])
    if (
        len(header) < BLOCK_HEADER
        or header[:4] != BLOCK_START
        or header[10:16] != BLOCK_EXTRA
    ):
        raise ValueError(f'{where}: no BGZF block starts there')
    return struct.unpack_from('<H', header, 16)[0] + 1


def blocks(chunks: Iterator[bytes], where: str) -> Iterator[tuple[int, int, bytes]]:
    """
    Each BGZF block of a stream: its offset in the stream, its size and its text.

    Raises:
        ValueError: if the stream is not whole blocks from its first byte to its last
    """
    held = b''
    offset = 0
    for chunk in chunks:
        held = held + chunk
        view = memoryview(held)
        at = 0
        while len(held) - at >= BLOCK_HEADER:
            size = block_bytes(view, at, f'{where}, byte {offset + at} of the range')
            if at + size > len(held):
                break
            yield offset + at, size, zlib.decompress(view[at : at + size], wbits=31)
            at += size
        offset += at
        held = held[at:]
    if held:
        raise ValueError(
            f'{where}: the range ends inside a block, {len(held)} bytes after byte '
            f'{offset}'
        )


def sample_blocks(chromosome: Indexed, places: int) -> list[int]:
    """
    The blocks of a chromosome that the check reads: the one of its first row, and
    others spread over its windows.

    The block of a chromosome's first row also holds the last rows of the chromosome
    before it, and the last window is always among the ones taken, so a sample holds
    rows of both ends of every chromosome.

    Raises:
        ValueError: if the index puts a window of the chromosome outside the
            chromosome's range
    """
    first, last = chromosome.start >> 16, chromosome.end >> 16
    held = sorted({window >> 16 for window in chromosome.windows if window})
    if held and not first <= held[0] <= held[-1] <= last:
        raise ValueError(
            f'the index puts a window of {chromosome.name} outside the range of its '
            f'rows'
        )
    if len(held) > places:
        step = (len(held) - 1) / (places - 1)
        held = [held[round(place * step)] for place in range(places)]
    return sorted({first, *held})


def locate(source: Source, places: int) -> dict[str, Any]:
    """
    The plan of one table: what `rows` and `sample` read, and what the receipt says
    of the source.

    Args:
        source: DeepMind's zip
        places: how many of a chromosome's 16 kb windows the check reads rows from,
            beside the chromosome's first rows

    Raises:
        ValueError: as `table_of` and `indexed` do, and if the chromosomes' rows are
            not the whole table after its header, one chromosome after the other
    """
    if places < 2:
        raise ValueError(
            f'a check reads at least 2 places of a chromosome, not {places}'
        )
    where = source.described.url
    table = table_of(source)
    chromosomes = indexed(table.index)
    if not chromosomes:
        raise ValueError(f'{where}: the index has no chromosome')
    # The header is what the table holds before its first row. bgzip can end a block
    # with it, so the rows may start in a later block than the first.
    start_block, start_text = divmod(chromosomes[0].start, MAX_BLOCK)
    if start_block > HEADER_LIMIT:
        raise ValueError(
            f'{where}: the rows of {table.member} start {start_block} bytes into it: '
            f'is that a header?'
        )
    header = b''
    header_end = min(start_block + MAX_BLOCK, table.size) if start_text else start_block
    header_bytes = stream(source, table.start, table.start + header_end)
    for offset, _, text in blocks(header_bytes, f'{where}, the header'):
        header += text[:start_text] if offset == start_block else text
        if offset == start_block:
            break
    if (
        not header.startswith(b'#CHROM\t')
        or header.count(b'\n') != 1
        or not header.endswith(b'\n')
    ):
        raise ValueError(
            f'{where}: the rows of {table.member} do not start after one header line'
        )
    for before, after in itertools.pairwise(chromosomes):
        if before.end != after.start:
            raise ValueError(
                f'{where}: the rows of {after.name} do not start where those of '
                f"{before.name} end, so the chromosomes' ranges are not the whole "
                f'table'
            )
    # A chromosome's range ends at the start of the block after its last row, and
    # after the last block comes the empty one that ends the file.
    end = (table.size - len(END_OF_FILE)) * MAX_BLOCK
    marker_start = table.start + table.size - len(END_OF_FILE)
    marker = source.read(marker_start, table.start + table.size)
    if chromosomes[-1].end != end or marker != END_OF_FILE:
        raise ValueError(
            f'{where}: the rows of {chromosomes[-1].name} do not end where '
            f'{table.member} ends'
        )
    return {
        'source': dataclasses.asdict(source.described),
        'member': table.member,
        'member_start': table.start,
        'member_bytes': table.size,
        'index_member': table.index_member,
        'header': header.decode(),
        'rows': sum(chromosome.rows for chromosome in chromosomes),
        'chromosomes': [
            {
                'name': chromosome.name,
                'rows': chromosome.rows,
                'start': list(divmod(chromosome.start, MAX_BLOCK)),
                'end': list(divmod(chromosome.end, MAX_BLOCK)),
                'bytes': (chromosome.end >> 16) - (chromosome.start >> 16),
                'sample_blocks': sample_blocks(chromosome, places),
            }
            for chromosome in chromosomes
        ],
    }


def planned_source(plan: Mapping[str, Any], get: Fetch = fetch) -> Source:
    """
    The table's source as the plan describes it, at the plan's generation.

    Raises:
        UnreadableError: if that generation is gone, or the source no longer has the
            plan's size
    """
    described = plan['source']
    source = open_source(described['url'], described['generation'], get)
    if source.described.size != described['size']:
        raise UnreadableError(
            f'{described["url"]} has {source.described.size} bytes, and '
            f'{described["size"]} when the run began'
        )
    return source


def rows(
    source: Source, plan: Mapping[str, Any], chromosome: str, out: BinaryIO
) -> int:
    """
    Writes the table's header and the rows of one chromosome as text.

    Returns:
        the number of rows written

    Raises:
        ValueError: if the plan has no such chromosome, or what the range holds is
            not the chromosome's rows as the index counts them
    """
    spans = {span['name']: span for span in plan['chromosomes']}
    if chromosome not in spans:
        raise ValueError(
            f'the table has no rows on {chromosome}: it has {", ".join(spans)}'
        )
    span = spans[chromosome]
    (start_block, start_text), (end_block, end_text) = span['start'], span['end']
    where = f'{plan["source"]["url"]}, {chromosome}'
    key = f'{chromosome}\t'.encode()
    # The range ends inside a block when the next chromosome's rows start in it.
    # That block is read whole, and a block has at most 64 KiB.
    last = min(end_block + MAX_BLOCK, plan['member_bytes']) if end_text else end_block
    out.write(plan['header'].encode())
    written = 0
    begun = False
    tail = b''
    reached = start_block
    base = plan['member_start']
    span_bytes = stream(source, base + start_block, base + last)
    for offset, size, text in blocks(span_bytes, where):
        block = start_block + offset
        if block > end_block or (block == end_block and not end_text):
            break
        begin = start_text if block == start_block else 0
        stop = end_text if block == end_block else len(text)
        piece = text[begin:stop]
        if piece and not begun:
            begun = True
            if not piece.startswith(key):
                raise ValueError(
                    f'{where}: the range does not start with a row of the chromosome'
                )
        out.write(piece)
        written += piece.count(b'\n')
        tail = piece[-TAIL:] if len(piece) >= TAIL else (tail + piece)[-TAIL:]
        reached = block + size
        if block == end_block:
            reached = end_block
            break
    if reached != end_block:
        raise ValueError(
            f'{where}: the blocks of the range end at byte {reached} of the table, '
            f'not at {end_block}'
        )
    if not tail.endswith(b'\n') or not tail[:-1].rpartition(b'\n')[2].startswith(key):
        raise ValueError(
            f'{where}: the range does not end with a row of the chromosome'
        )
    if written != span['rows']:
        raise ValueError(
            f'{where}: the range holds {written} rows, and the index counts '
            f'{span["rows"]}'
        )
    return written


def sample(source: Source, plan: Mapping[str, Any], out: BinaryIO) -> int:
    """
    Writes the header and the rows of the plan's sample blocks, in the table's order.

    A block's first line is left out, because a row can begin in the block before,
    and so is what follows its last line break.

    Returns:
        the number of rows written
    """
    base, end = plan['member_start'], plan['member_bytes']
    where = plan['source']['url']
    out.write(plan['header'].encode())
    written = 0
    sampled = {
        block for span in plan['chromosomes'] for block in span['sample_blocks']
    }
    for block in sorted(sampled):
        data = source.read(base + block, base + min(block + MAX_BLOCK, end))
        size = block_bytes(data, 0, f'{where}, byte {block} of {plan["member"]}')
        lines = zlib.decompress(data[:size], wbits=31).split(b'\n')[1:-1]
        out.writelines(line + b'\n' for line in lines)
        written += len(lines)
    if not written:
        raise ValueError(f'{where}: the sample blocks hold no whole row')
    return written


# ---- merging per-chromosome echtvar archives ----


def other_extra(extra: bytes, where: str) -> bytes:
    """
    An entry's extra fields without the zip64 one, whose values `zipfile` has already
    read into the entry.

    Raises:
        ValueError: if the bytes are not whole extra fields
    """
    kept = b''
    at = 0
    while at < len(extra):
        if at + 4 > len(extra):
            raise ValueError(f'{where}: an extra field is cut short')
        kind, size = struct.unpack_from('<2H', extra, at)
        if at + 4 + size > len(extra):
            raise ValueError(f'{where}: an extra field is cut short')
        if kind != ZIP64_EXTRA:
            kept += extra[at : at + 4 + size]
        at += 4 + size
    return kept


def copy_entry(
    archive: BinaryIO, info: zipfile.ZipInfo, out: BinaryIO, where: str
) -> None:
    """
    Copies one entry of a zip, its local header and its data, as the bytes they are.

    Raises:
        ValueError: if the entry's local header is not the one that the zip's
            directory describes, or the zip ends inside the entry
    """
    archive.seek(info.header_offset)
    header = archive.read(LOCAL_HEADER.size)
    if len(header) != LOCAL_HEADER.size or header[:4] != b'PK\x03\x04':
        raise ValueError(f'{where}: no entry starts at byte {info.header_offset}')
    flags = LOCAL_HEADER.unpack(header)[2]
    name_length, extra_length = LOCAL_HEADER.unpack(header)[-2:]
    # With this flag the sizes follow the data, in a record that the copy would
    # leave behind. echtvar writes them into the header.
    if flags & 0x08:
        raise ValueError(f'{where}: the sizes of the entry follow its data')
    name = archive.read(name_length)
    if name != info.filename.encode():
        raise ValueError(
            f'{where}: the entry at byte {info.header_offset} is named {name!r}'
        )
    out.write(header)
    out.write(name)
    left = extra_length + info.compress_size
    while left:
        data = archive.read(min(left, COPY_STEP))
        if not data:
            raise ValueError(f'{where}: the zip ends inside the entry')
        out.write(data)
        left -= len(data)


def central_record(
    info: zipfile.ZipInfo, offset: int, where: str, zip64_from: int = ZIP64_BYTES
) -> bytes:
    """
    The record of one entry in a zip's central directory, for the entry's local
    header at `offset`. The record is the entry's own with the new offset.

    Args:
        info: the entry, as the directory of its own zip describes it
        offset: where the entry's local header starts in the new zip
        where: the entry's place, for a message
        zip64_from: the smallest size or offset that goes into the zip64 field

    Raises:
        ValueError: if the entry has an extra field or a comment, which echtvar does
            not write and this function would not carry over
    """
    if other_extra(info.extra, where) or info.comment:
        raise ValueError(f'{where}: the entry has an extra field or a comment')
    name = info.filename.encode()
    # The zip64 field holds, in this order, the values that the record's 32 bits do
    # not hold. The record holds 0xFFFFFFFF in the place of each.
    values = (info.file_size, info.compress_size, offset)
    large = [value for value in values if value >= zip64_from]
    file_size, compress_size, start = (
        ZIP64_BYTES if value >= zip64_from else value for value in values
    )
    extra = b''
    version = info.extract_version
    if large:
        extra = struct.pack(f'<2H{len(large)}Q', ZIP64_EXTRA, 8 * len(large), *large)
        version = max(version, ZIP64_VERSION)
    year, month, day, hour, minute, second = info.date_time
    return (
        CENTRAL_HEADER.pack(
            b'PK\x01\x02',
            info.create_system << 8 | info.create_version,
            version,
            info.flag_bits,
            info.compress_type,
            hour << 11 | minute << 5 | second // 2,
            (year - 1980) << 9 | month << 5 | day,
            info.CRC,
            compress_size,
            file_size,
            len(name),
            len(extra),
            0,
            0,
            info.internal_attr,
            info.external_attr,
            start,
        )
        + name
        + extra
    )


def end_records(
    entries: int, size: int, start: int, zip64_from: int = ZIP64_BYTES
) -> bytes:
    """
    The records that end a zip whose central directory has `entries` records in
    `size` bytes from byte `start`: the zip64 records where a value needs them, laid
    out as the zip library of echtvar writes them, and the end record.

    Args:
        entries: the number of entries of the zip
        size: the bytes of the central directory
        start: where the central directory starts
        zip64_from: the smallest size or offset that goes into the zip64 records
    """
    records = b''
    if entries >= ZIP64_ENTRIES or max(size, start) >= zip64_from:
        records += END_RECORD_64.pack(
            b'PK\x06\x06', 44, 46, 46, 0, 0, entries, entries, size, start
        )
        records += END_LOCATOR_64.pack(b'PK\x06\x07', 0, start + size, 1)
    counted = min(entries, ZIP64_ENTRIES)
    records += END_RECORD.pack(
        b'PK\x05\x06',
        0,
        0,
        counted,
        counted,
        ZIP64_BYTES if size >= zip64_from else size,
        ZIP64_BYTES if start >= zip64_from else start,
        0,
    )
    return records


def merge(
    archives: Sequence[Path], out: Path, zip64_from: int = ZIP64_BYTES
) -> int:
    """
    Writes one echtvar archive with the entries of several, each of which was
    encoded from other chromosomes with the same configuration.

    The entries keep their order: the directory entry and `config.json` of the first
    archive, then the blocks of each archive in the order given. That is the order
    in which one `echtvar encode` over the same inputs writes them. `echtvar anno`
    opens entries by name, so the order does not matter to it.

    Args:
        archives: the archives, in the order of their chromosomes
        out: the archive to write
        zip64_from: the smallest size or offset that is written in the zip64 form.
            A test gives a small one, to write that form for a small archive.

    Returns:
        the number of entries written

    Raises:
        ValueError: if an archive is not as `echtvar encode` writes one, the
            archives' `config.json` differ, two archives hold an entry of the same
            name, or an archive has a string field
    """
    if not archives:
        raise ValueError('no archive to merge')
    config = None
    written: list[tuple[zipfile.ZipInfo, int]] = []
    names: dict[str, Path] = {}
    with out.open('wb') as merged:
        for path in archives:
            with zipfile.ZipFile(path) as archive, path.open('rb') as raw:
                infos = archive.infolist()
                first = [info.filename for info in infos[:2]]
                if first != [ARCHIVE_ROOT, ARCHIVE_CONFIG]:
                    raise ValueError(
                        f'{path}: its first entries are {first}, and not '
                        f'{ARCHIVE_ROOT} and {ARCHIVE_CONFIG} as echtvar writes them'
                    )
                this_config = archive.read(ARCHIVE_CONFIG)
                if config is None:
                    config = this_config
                elif this_config != config:
                    raise ValueError(
                        f'{path}: its {ARCHIVE_CONFIG} differs from that of '
                        f'{archives[0]}, so the two were not encoded alike'
                    )
                for number, info in enumerate(infos):
                    where = f'{path}, {info.filename}'
                    if number < len(first):
                        # The directory entry and the configuration are written
                        # once, from the first archive.
                        if path != archives[0]:
                            continue
                    elif info.filename.startswith(ARCHIVE_STRINGS):
                        raise ValueError(
                            f'{where}: the archive has a string field, whose '
                            f'numbers differ between two encodes'
                        )
                    elif not ARCHIVE_ENTRY.fullmatch(info.filename):
                        raise ValueError(f'{where}: not an entry that echtvar writes')
                    if info.filename in names:
                        raise ValueError(
                            f'{where}: {names[info.filename]} holds an entry of the '
                            f'same name'
                        )
                    names[info.filename] = path
                    offset = merged.tell()
                    copy_entry(raw, info, merged, where)
                    written.append((info, offset))
        start = merged.tell()
        for info, offset in written:
            merged.write(central_record(info, offset, info.filename, zip64_from))
        size = merged.tell() - start
        merged.write(end_records(len(written), size, start, zip64_from))
    return len(written)


def verify_archive(path: Path, entries: int) -> None:
    """
    Inflates every entry of a merged archive and compares the bytes with the size
    and the checksum that echtvar recorded for the entry when it encoded it.

    Raises:
        ValueError: if the archive does not list `entries` entries, or an entry does
            not inflate to the bytes of its checksum
    """
    with zipfile.ZipFile(path) as archive:
        listed = len(archive.infolist())
        if listed != entries:
            raise ValueError(
                f'{path} lists {listed} entries, and {entries} were written'
            )
        try:
            bad = archive.testzip()
        except zlib.error as error:
            raise ValueError(f'{path}: an entry does not inflate: {error}') from error
    if bad is not None:
        raise ValueError(f'{path}: the entry {bad} does not match its checksum')


# ---- the receipt of a finished archive ----


def terms_text(plan: Mapping[str, Any]) -> str:
    """The text of TERMS.txt for an archive of a table that is not for permissive
    use."""
    table = plan['table']
    return (
        f'{table["name"]}.zip holds the values of "{table["title"]}", a download of '
        f"Google DeepMind's AlphaGenome that its downloads page lists under "
        f'"Downloadable artifacts for non-commercial use only". The values are '
        f'subject to the AlphaGenome Output Terms of Use, and so is every file that '
        f'holds values read from this archive.\n\n{table["notice"]}\n'
    )


def tool_versions() -> dict[str, str]:
    """The versions of the tools in the image that this runs in."""

    def said(*command: str) -> list[str]:
        done = subprocess.run(command, check=True, capture_output=True, text=True)
        return done.stdout.splitlines()

    bcftools = said('bcftools', '--version')
    return {
        'case-annotation-tools': importlib.metadata.version('case-annotation-tools'),
        'echtvar': said('echtvar', '--version')[0],
        'bcftools': bcftools[0],
        'htslib': bcftools[1],
        'python': sys.version.split()[0],
    }


def receipt(
    plan: Mapping[str, Any],
    archive_bytes: int,
    sha256: str,
    check: str,
    versions: Mapping[str, str],
) -> dict[str, Any]:
    """What an archive was built from, with what, and by which run."""
    table = plan['table']
    now = datetime.datetime.now(datetime.timezone.utc)
    return {
        'archive': table['name'],
        'file': {
            'name': f"{table['name']}.zip",
            'bytes': archive_bytes,
            'sha256': sha256,
        },
        'written': now.isoformat(timespec='seconds'),
        'source': plan['source'],
        'reference_index': plan['reference_index'],
        'table': {
            'title': table['title'],
            'member': plan['member'],
            'bytes': plan['member_bytes'],
            'rows': plan['rows'],
            'chromosomes': [
                {key: span[key] for key in ('name', 'rows', 'bytes')}
                for span in plan['chromosomes']
            ],
        },
        'tools': {**plan['tools'], 'versions': dict(versions)},
        'script': plan['script'],
        'check': check,
        'terms': table['notice'],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description='Steps of alphagenome_echtvar.py')
    commands = parser.add_subparsers(dest='command', required=True)
    for name, text in (
        ('rows', "write the header and one chromosome's rows"),
        ('sample', 'write the rows that the check compares'),
    ):
        reading = commands.add_parser(name, help=f'{text} on standard output')
        reading.add_argument('--plan', type=Path, required=True)
        if name == 'rows':
            reading.add_argument('--chromosome', required=True)
    merging = commands.add_parser('merge', help='merge per-chromosome archives')
    merging.add_argument('--out', type=Path, required=True)
    merging.add_argument('archives', type=Path, nargs='+')
    recording = commands.add_parser(
        'records', help="write a checked archive's configuration, receipt and terms"
    )
    recording.add_argument('--plan', type=Path, required=True)
    recording.add_argument('--archive', type=Path, required=True)
    recording.add_argument('--sha256', required=True, help="the archive's SHA-256")
    recording.add_argument('--check', type=Path, required=True, help='the report')
    recording.add_argument('--config', type=Path, required=True, help='to write')
    recording.add_argument('--receipt', type=Path, required=True, help='to write')
    recording.add_argument('--terms', type=Path, help='to write, for a notice')
    args = parser.parse_args(argv)
    try:
        if args.command == 'merge':
            entries = merge(args.archives, args.out)
            verify_archive(args.out, entries)
            print(
                f'{args.out}: {entries} entries of {len(args.archives)} archives, '
                f'each equal to its checksum',
                file=sys.stderr,
            )
            return
        plan = json.loads(args.plan.read_text('utf-8'))
        if args.command == 'records':
            if (args.terms is None) != (plan['table']['notice'] is None):
                raise ValueError(
                    '--terms is given for a table with a notice, and only for one'
                )
            with zipfile.ZipFile(args.archive) as archive:
                args.config.write_bytes(archive.read(ARCHIVE_CONFIG))
            written = receipt(
                plan,
                args.archive.stat().st_size,
                args.sha256,
                args.check.read_text('utf-8').strip(),
                tool_versions(),
            )
            args.receipt.write_text(json.dumps(written, indent=2) + '\n', 'utf-8')
            if args.terms is not None:
                args.terms.write_text(terms_text(plan), 'utf-8')
            return
        source = planned_source(plan)
        if args.command == 'rows':
            count = rows(source, plan, args.chromosome, sys.stdout.buffer)
        else:
            count = sample(source, plan, sys.stdout.buffer)
        sys.stdout.buffer.flush()
        print(f'{count} rows of {plan["member"]}', file=sys.stderr)
    except BrokenPipeError:
        # The command that read the rows has ended, and its own message says why.
        # Python would report the closed pipe once more when it exits, unless
        # standard output leads elsewhere by then.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        raise SystemExit(
            'error: the command that read the rows ended before the last row'
        ) from None
    except (UnreadableError, ValueError, zipfile.BadZipFile) as error:
        raise SystemExit(f'error: {error}') from error


if __name__ == '__main__':
    main()
