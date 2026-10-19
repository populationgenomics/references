"""
Tests for extract.py: the P value filter, the per-file extraction, the tabix index
reader and the command-line guards.

From the repo root:

    uv run --no-project --with pytest --with pandas==2.2.3 --with pyarrow==18.1.0 \
        --with polars==1.34.0 --with pysam==0.23.3 pytest gwas_sumstats/test_extract.py
"""

import gzip
import json
import os
import sys
from pathlib import Path

import download
import extract
import pandas as pd
import pyarrow.parquet as pq
import pytest

import format as fmt

HEADER = '\t'.join(fmt.OUTPUT_COLUMNS)


def row(**fields) -> str:
    """One formatted row, as format.py writes it: every column, NA where empty."""
    values = {
        'chromosome': '1',
        'base_pair_location': '100',
        'effect_allele': 'A',
        'other_allele': 'G',
        'beta': '0.01',
        'standard_error': '0.01',
        'effect_allele_frequency': '0.3',
        'p_value': '0.9',
        'neg_log_10_p_value': '0.045757',
        'rsid': 'rs1',
        'n': '1000',
        'z': 'NA',
        'effect_allele_is_alt': 'true',
    }
    values.update({key: str(value) for key, value in fields.items()})
    return '\t'.join(values[column] for column in fmt.OUTPUT_COLUMNS)


def hit(rsid: str = 'rs1', **fields) -> str:
    """A row well past the threshold."""
    return row(
        rsid=rsid, beta='0.5', p_value='1e-20', neg_log_10_p_value='20', **fields
    )


def write_file(path, rows: list[str]) -> None:
    """A minimal formatted file: format.py's header plus the given rows."""
    with gzip.open(path, 'wt') as handle:
        handle.write(HEADER + '\n')
        for line in rows:
            handle.write(line + '\n')


@pytest.mark.parametrize(
    'p_text, neg_log_10, z, beta, se, expected',
    [
        ('1e-9', '9', None, None, None, True),  # well below the cutoff
        ('1e-5', '5', None, None, None, True),  # exactly at it
        ('0.01', '2', None, None, None, False),  # above it
        ('8.51e-3553', '3552.07', None, None, None, True),  # below the float range
        ('8.51e-3553', 'NA', None, None, None, True),  # same, from the text: reads as 0
        ('NA', 'NA', None, None, None, False),  # no P value and no effect size
        ('1e-9', None, None, None, None, True),  # no neg_log column: uses the text
        ('1e-9', 'NA', None, None, None, True),  # neg_log missing for this row
        ('1e-9', 'nan', None, None, None, True),  # neg_log unreadable for this row
        ('0.01', '9', None, None, None, True),  # neg_log wins where both are readable
        ('not a number', 'NA', None, None, None, False),  # unparseable, no effect
        # format.py writes a published P of 0 as NA in both columns: kept on z
        ('NA', 'NA', None, '-0.4774', '0.01', True),  # |beta/se| = 47.7
        ('NA', 'NA', None, '0.04', '0.01', False),  # |beta/se| = 4.0 < 4.42
        ('NA', 'NA', '-5.0', 'NA', 'NA', True),  # z where there is no beta
        ('NA', 'NA', '4.0', '1.0', '0.01', False),  # z preferred to beta/se
        ('0.5', '0.30103', None, '1.0', '0.01', False),  # a readable P wins over z
    ],
)
def test_threshold(p_text, neg_log_10, z, beta, se, expected):
    def series(value):
        return None if value is None else pd.Series([value])

    kept = extract.passing(
        pd.Series([p_text]),
        series(neg_log_10),
        1e-5,
        z=series(z),
        beta=series(beta),
        se=series(se),
    )
    assert kept.tolist() == [expected]


@pytest.fixture
def fasta(tmp_path):
    """GRCh38 stand-in: contigs 1 and 2, all 'C' except 1:10 = A."""
    pysam = pytest.importorskip('pysam')
    bases = ['C'] * 1000
    bases[10 - 1] = 'A'
    path = tmp_path / 'genome.fa'
    path.write_text(f'>1\n{"".join(bases)}\n>2\n{"C" * 1000}\n')
    pysam.faidx(str(path))
    return path


def test_underflowed_p_from_format_py_is_kept(fasta, tmp_path):
    """
    End to end through format.py: a published P of exactly 0 (APOE in Timsina 2026)
    comes out NA in both P columns, and must still be among the hits.
    """
    pytest.importorskip('polars')
    source = tmp_path / 'source.tsv'
    source.write_text(
        'chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t'
        'standard_error\tp_value\n'
        '1\t10\tA\tC\t-0.4774\t0.01\t0.0\n'
        '1\t20\tC\tT\t0.01\t0.01\t0.5\n'
    )
    row_csv = {
        'file_id': 'test',
        'columns': '',
        'n_study': '100',
        'source_build': 'GRCh38',
    }
    record = {'md5': 'x', 'md5_checked': False, 'downloaded_at': 'x'}
    out_root = tmp_path / 'formatted'
    fmt.format_one(
        row_csv, source, out_root, tmp_path / 'stats.json', fasta, {}, record
    )
    assert json.loads((tmp_path / 'stats.json').read_text())['p_value_zero'] == 1

    dest = tmp_path / 'hits.parquet'
    counts = extract.extract_one(Path(f'{out_root}.tsv.gz'), dest, 'test', 1e-5)

    hits = pd.read_parquet(dest)
    assert hits['base_pair_location'].tolist() == ['10']
    assert hits['p_value'].isna().all()
    assert counts['rows_kept_without_p'] == 1


def test_extract_keeps_only_hits_and_records_how(tmp_path):
    source = tmp_path / 'in.tsv.gz'
    write_file(
        source,
        [
            hit('rs1'),
            row(rsid='rs2'),
            row(rsid='rs3', chromosome='2', p_value='1e-6', neg_log_10_p_value='6'),
        ],
    )
    dest = tmp_path / 'out.parquet'
    counts = extract.extract_one(source, dest, 'STUDY_1', 1e-5)

    assert counts == {
        'file_id': 'STUDY_1',
        'rows_in': 3,
        'rows_kept': 2,
        'rows_kept_without_p': 0,
    }
    hits = pd.read_parquet(dest)
    assert list(hits['rsid']) == ['rs1', 'rs3']
    assert set(hits['file_id']) == {'STUDY_1'}
    record = json.loads(pq.read_schema(dest).metadata[b'gwas_sumstats'])
    assert record == {**counts, 'threshold': 1e-5}


def test_extract_handles_a_file_with_no_hits(tmp_path):
    """An empty result still gets a parquet with the same all-text columns, so that
    reading the whole folder back never has to special-case it."""
    source = tmp_path / 'in.tsv.gz'
    write_file(source, [row()])
    dest = tmp_path / 'out.parquet'

    counts = extract.extract_one(source, dest, 'STUDY_2', 1e-5)
    assert counts['rows_kept'] == 0
    schema = pq.read_schema(dest)
    assert schema.names == ['file_id', *fmt.OUTPUT_COLUMNS]
    assert {str(field.type) for field in schema} == {'string'}


@pytest.mark.parametrize('empty_name', ['A_first.parquet', 'Z_last.parquet'])
def test_folder_with_an_empty_file_reads_as_one_table(tmp_path, empty_name):
    """
    Files with and without hits share one all-text schema, so the folder reads back
    whichever sorts first. Inferred types made the empty file's columns null, which
    failed the read when that file came first.
    """
    folder = tmp_path / 'hits'
    folder.mkdir()
    with_hits, no_hits = tmp_path / 'with.tsv.gz', tmp_path / 'without.tsv.gz'
    write_file(with_hits, [hit('rs1')])
    write_file(no_hits, [row()])
    full_name = (
        'Z_last.parquet' if empty_name == 'A_first.parquet' else 'A_first.parquet'
    )
    extract.extract_one(with_hits, folder / full_name, 'WITH', 1e-5)
    extract.extract_one(no_hits, folder / empty_name, 'WITHOUT', 1e-5)

    table = pd.read_parquet(folder)
    assert table['rsid'].tolist() == ['rs1']
    assert table['file_id'].tolist() == ['WITH']


def test_extract_reads_across_chunks(tmp_path, monkeypatch):
    """Chunking must not lose or duplicate rows, whatever the chunk boundary."""
    monkeypatch.setattr(extract, 'CHUNK_ROWS', 2)
    source = tmp_path / 'in.tsv.gz'
    write_file(
        source,
        [
            hit(f'rs{i}', base_pair_location=100 + i) if i % 2 else row(rsid=f'rs{i}')
            for i in range(10)
        ],
    )
    dest = tmp_path / 'out.parquet'

    counts = extract.extract_one(source, dest, 'STUDY_3', 1e-5)
    assert (counts['rows_in'], counts['rows_kept']) == (10, 5)


def test_formatted_name_ignores_the_manifest_column():
    """
    The manifest's formatted_file says 'out.tsv.gz' for every row, so the name has
    to be rebuilt from file_id.
    """
    assert extract.formatted_name('2017_Wheeler_PLoSMed_HbA1c_AFR_GCST007952') == (
        '2017_Wheeler_PLoSMed_HbA1c_AFR_GCST007952_GRCh38_formatted.tsv.gz'
    )


def test_read_files_csv_rejects_an_unknown_id(tmp_path):
    files_csv = tmp_path / 'files.csv'
    files_csv.write_text('file_id,biomarker\nA,HbA1c\nB,Urate\n')

    assert len(download.read_files_csv(files_csv)) == 2
    assert [r['file_id'] for r in download.read_files_csv(files_csv, ['B'])] == ['B']
    with pytest.raises(ValueError, match='not in files.csv'):
        download.read_files_csv(files_csv, ['C'])


def test_tabix_index_rows_and_windows_match_pysam(tmp_path):
    """The index parser lists the same contigs, in order, as pysam does, with the
    rows on each and the 16-kb windows that hold them, sparse or dense."""
    pysam = pytest.importorskip('pysam')
    x_positions = [500, 501, 16_384, 16_385, 1_000_000] + list(
        range(2_000_000, 2_060_000, 50)
    )
    plain = tmp_path / 'formatted.tsv'
    lines = [hit(chromosome='1')]
    lines += [hit(chromosome='2', base_pair_location=300 + i) for i in range(2)]
    lines += [hit(chromosome='23', base_pair_location=p) for p in x_positions]
    plain.write_text(HEADER + '\n' + '\n'.join(lines) + '\n')
    compressed = pysam.tabix_index(
        str(plain), seq_col=0, start_col=1, end_col=1, line_skip=1, force=True
    )
    with pysam.TabixFile(compressed) as tabix:
        expected = list(tabix.contigs)

    with Path(f'{compressed}.tbi').open('rb') as handle:
        index = extract.read_tabix_index(handle)
    assert list(index) == expected == ['1', '2', '23']
    assert {name: contig['rows'] for name, contig in index.items()} == {
        '1': 1,
        '2': 2,
        '23': len(x_positions),
    }
    assert index['23']['windows'] == {(p - 1) >> 14 for p in x_positions}


def test_tabix_index_rejects_other_files(tmp_path):
    other = tmp_path / 'not_an_index.gz'
    with gzip.open(other, 'wb') as handle:
        handle.write(b'CSI\x01' + bytes(64))
    with other.open('rb') as handle, pytest.raises(ValueError, match='not a tabix'):
        extract.read_tabix_index(handle)


# hg19 -> hg38 in UCSC chain form. chr7 lands on X at 100,000-100,050 (window 6),
# chrY on X's minus strand, X on X: only the first two are liftover strays.
CHAIN = """\
chain 1000 chr7 159345973 + 5000 5050 chrX 156040895 + 100000 100050 1
50

chain 1000 chrY 59373566 + 0 20 chrX 156040895 - 1000 1010 2
5 5 5
5

chain 1000 chrX 155270560 + 0 3000000 chrX 156040895 + 0 3000000 3
3000000

chain 1000 chr7 159345973 + 0 100 chr7 159345973 + 0 100 4
100
"""


def test_liftover_windows_are_where_other_chromosomes_land_on_x():
    windows = extract.liftover_windows_on_x(CHAIN.splitlines())
    # minus strand: 1000-1005 and 1010-1015 from the end of X
    end = 156_040_895
    assert windows == {100_000 >> 14, (end - 1015) >> 14}


def test_liftover_windows_from_the_real_chain_hold_known_strays():
    """
    hg19 chr7:100614152, chr7:74223913, chr22:18047650 and chr1:249058621 each lift
    to GRCh38 X; all must fall in the windows the chain marks. Skipped unless the
    UCSC chain is given in GWAS_HG19_CHAIN.
    """
    path = os.environ.get('GWAS_HG19_CHAIN')
    if not path:
        pytest.skip('GWAS_HG19_CHAIN not set')
    with gzip.open(path, 'rt') as handle:
        windows = extract.liftover_windows_on_x(handle)
    for grch38 in (45216764, 92475099, 154768770, 122549712):
        assert (grch38 - 1) >> 14 in windows
    assert len(windows) < 50


FILE = {'file_id': 'F', 'biomarker': 'B', 'ancestry': 'EUR', 'source_build': 'GRCh37'}


@pytest.mark.parametrize(
    'x_windows, liftover_windows, has_x',
    [
        ({6}, {6}, False),  # only liftover strays: X never published
        ({6, 500, 501}, {6}, True),  # a published X, some of it in a stray window
        (set(), {6}, False),  # no X at all
        ({6}, set(), True),  # a GRCh38 source is not lifted, so any X is its own
    ],
)
def test_has_x_ignores_rows_liftover_put_on_x(x_windows, liftover_windows, has_x):
    contigs = {'1': {'rows': 9, 'windows': {0}}}
    if x_windows:
        contigs['23'] = {'rows': 3 * len(x_windows), 'windows': x_windows}
    record = extract.coverage_record(FILE, contigs, liftover_windows)
    assert record['has_x'] is has_x
    assert record['x_windows'] == len(x_windows)
    assert record['x_windows_from_source_x'] == len(x_windows - liftover_windows)


def test_coverage_lists_unformatted_files():
    record = extract.coverage_record(FILE, None, set())
    assert record['formatted'] is False
    assert record['has_x'] == ''


def test_passing_accepts_the_nullable_string_dtype():
    """Newer pandas reads text as its nullable string dtype, where NA is pd.NA."""
    kept = extract.passing(
        pd.Series(['1e-9', None, '0.5'], dtype='string'),
        pd.Series([None, None, None], dtype='string'),
        1e-5,
        beta=pd.Series([None, '1.0', None], dtype='string'),
        se=pd.Series([None, '0.01', None], dtype='string'),
    )
    assert kept.tolist() == [True, True, False]


def test_run_local_reads_local_files_in_place(tmp_path):
    """Local formatted files are not copied, so the output folder holds only
    parquets and reads back as one table."""
    formatted = tmp_path / 'v1'
    formatted.mkdir()
    write_file(formatted / extract.formatted_name('STUDY_4'), [hit('rs1')])
    out = tmp_path / 'hits'

    extract.run_local([{'file_id': 'STUDY_4'}], out, str(formatted), 1e-5)

    assert [path.name for path in out.iterdir()] == ['STUDY_4.parquet']
    assert pd.read_parquet(out)['rsid'].tolist() == ['rs1']


def test_run_local_redoes_files_only_with_force(tmp_path):
    formatted = tmp_path / 'v1'
    formatted.mkdir()
    source = formatted / extract.formatted_name('STUDY_5')
    out = tmp_path / 'hits'
    write_file(source, [hit('rs1')])
    extract.run_local([{'file_id': 'STUDY_5'}], out, str(formatted), 1e-5)

    write_file(source, [hit('rs2')])
    extract.run_local([{'file_id': 'STUDY_5'}], out, str(formatted), 1e-5)
    assert pd.read_parquet(out)['rsid'].tolist() == ['rs1']
    extract.run_local([{'file_id': 'STUDY_5'}], out, str(formatted), 1e-5, force=True)
    assert pd.read_parquet(out)['rsid'].tolist() == ['rs2']


def test_hits_older_than_their_source_are_redone(tmp_path):
    """After format.py --only <id> --force, the formatted file is newer than its
    hits, so the hits are no longer current."""
    source, dest = tmp_path / 'formatted.tsv.gz', tmp_path / 'hits.parquet'
    assert not extract.is_current(dest, source)
    source.write_text('x')
    dest.write_text('x')
    os.utime(source, (1000, 1000))
    os.utime(dest, (2000, 2000))
    assert extract.is_current(dest, source)
    os.utime(source, (3000, 3000))
    assert not extract.is_current(dest, source)


@pytest.mark.parametrize(
    'argv, message',
    [
        (['--force'], '--force needs --only'),
        (['--coverage', '--only', 'X'], '--coverage-out'),
        (['--threshold', '5e-8'], 'give another --out'),
    ],
)
def test_main_refuses_unsafe_combinations(monkeypatch, capsys, argv, message):
    monkeypatch.setattr(sys, 'argv', ['extract.py', *argv])
    with pytest.raises(SystemExit):
        extract.main()
    assert message in capsys.readouterr().err
