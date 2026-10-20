"""
Tests for extract.py: the P value filter, the per-file extraction, the chromosome
table and the command-line guards.

From the repo root:

    uv run --no-project --with pytest --with pandas==2.2.3 --with pyarrow==18.1.0 \
        --with polars==1.34.0 --with pysam==0.23.3 pytest gwas_sumstats/test_extract.py
"""

import gzip
import json
import os
import sys
from pathlib import Path

import extract
import pandas as pd
import polars  # noqa: F401  format.py imports it inside functions
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
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


def test_a_standard_error_of_zero_is_not_a_hit(tmp_path):
    """
    With no P value and no z, beta / 0 would be an infinite z and always pass. An SE
    of 0 is an error in the source, so the row is dropped and counted instead.
    """
    na_p = {'p_value': 'NA', 'neg_log_10_p_value': 'NA'}
    source = tmp_path / 'in.tsv.gz'
    write_file(
        source,
        [
            row(rsid='rs_se0', beta='1e-4', standard_error='0', **na_p),
            row(rsid='rs_seneg', beta='1e-4', standard_error='-0.01', **na_p),
            row(rsid='rs_real', beta='1.0', standard_error='0.01', **na_p),
        ],
    )
    counts = extract.extract_one(source, tmp_path / 'out.parquet', 'S', 1e-5)

    assert pd.read_parquet(tmp_path / 'out.parquet')['rsid'].tolist() == ['rs_real']
    assert counts['rows_dropped_se_zero'] == 2
    assert counts['rows_kept_without_p'] == 1


@pytest.mark.parametrize(
    'cutoff, neg_log_10, expected',
    [
        (2e-5, '4.69897', True),  # P = 2e-5 as format.py stores it: at the cutoff
        (2e-5, '4.698969', False),  # just above the cutoff
        (1e-5, '5', True),  # exactly at the default: kept, so the test needs >=
        (1e-5, '4.999999', False),
        (5e-8, '7.30103', True),
        (5e-8, '7.301029', False),
    ],
)
def test_a_p_value_exactly_at_the_cutoff_is_kept(cutoff, neg_log_10, expected):
    na = pd.Series(['NA'])
    kept = extract.passing(na, pd.Series([neg_log_10]), cutoff)
    assert kept.tolist() == [expected]


def test_z_cutoff_stays_exact_far_below_1e_16():
    """At 1e-20, 1 - cutoff / 2 rounds to 1 and inv_cdf fails; the lower tail does
    not. |z| = 9.3 is P of about 1.4e-20, so it passes 1e-19 but not 1e-21."""
    beta, se = pd.Series(['9.3']), pd.Series(['1'])
    na = pd.Series(['NA'])
    assert extract.passing(na, na, 1e-19, beta=beta, se=se).tolist() == [True]
    assert extract.passing(na, na, 1e-21, beta=beta, se=se).tolist() == [False]


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
        ('8.51e-3553', 'NA', None, None, None, True),  # same, from the text: reads as 0
        ('NA', 'NA', None, None, None, False),  # no P value and no effect size
        ('1e-9', None, None, None, None, True),  # no neg_log column: uses the text
        ('1e-9', 'NA', None, None, None, True),  # neg_log missing for this row
        ('0.01', '9', None, None, None, True),  # neg_log wins where both are readable
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
    assert hits['base_pair_location'].tolist() == [10]
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
        'rows_dropped_se_zero': 0,
        'rows_per_chromosome': {'1': 2, '2': 1},  # every row read, not only hits
    }
    hits = pd.read_parquet(dest)
    assert list(hits['rsid']) == ['rs1', 'rs3']
    assert set(hits['file_id']) == {'STUDY_1'}
    record = json.loads(pq.read_metadata(dest).metadata[b'gwas_sumstats'])
    assert record == {**counts, 'threshold': 1e-5}


def test_extract_handles_a_file_with_no_hits(tmp_path):
    """
    An empty result still gets a parquet with the same typed columns, so the folder
    reads back as one table even when the empty file sorts first. Types inferred per
    file made an empty file's columns null, which failed that read.
    """
    folder = tmp_path / 'hits'
    folder.mkdir()
    no_hits, with_hits = tmp_path / 'without.tsv.gz', tmp_path / 'with.tsv.gz'
    write_file(no_hits, [row()])
    write_file(with_hits, [hit('rs1')])

    counts = extract.extract_one(no_hits, folder / 'A_empty.parquet', 'WITHOUT', 1e-5)
    extract.extract_one(with_hits, folder / 'B_full.parquet', 'WITH', 1e-5)

    assert counts['rows_kept'] == 0
    schema = pq.read_schema(folder / 'A_empty.parquet')
    assert schema.names == ['file_id', *fmt.OUTPUT_COLUMNS]
    assert [field.type for field in schema] == [
        pa.type_for_alias(kind) for kind in extract.HIT_TYPES.values()
    ]
    table = pd.read_parquet(folder)
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
    # Each chunk's hits are written as they come, one row group per chunk with hits.
    assert pd.read_parquet(dest)['rsid'].tolist() == ['rs1', 'rs3', 'rs5', 'rs7', 'rs9']
    assert pq.read_metadata(dest).num_row_groups == 5


def test_extract_handles_a_file_with_a_header_and_no_rows(tmp_path):
    source, dest = tmp_path / 'in.tsv.gz', tmp_path / 'out.parquet'
    write_file(source, [])

    counts = extract.extract_one(source, dest, 'EMPTY', 1e-5)
    assert (counts['rows_in'], counts['rows_kept']) == (0, 0)
    empty = pd.read_parquet(dest)
    assert len(empty) == 0 and list(empty.columns) == ['file_id', *HEADER.split('\t')]


def test_pins_match_everywhere():
    # Batch jobs install PIP_PACKAGES; the README and this module's docstring must
    # tell people to test and run locally with the same versions.
    readme = Path(extract.__file__).with_name('README.md').read_text()
    for pin in extract.PIP_PACKAGES.split():
        assert pin in __doc__, f'{pin} missing from the test_extract.py docstring'
        assert readme.count(pin) >= 2, f'{pin} missing from a README command'


def test_formatted_name_ignores_the_manifest_column():
    """
    The manifest's formatted_file says 'out.tsv.gz' for every row, so the name has
    to be rebuilt from file_id.
    """
    assert extract.formatted_name('2017_Wheeler_PLoSMed_HbA1c_AFR_GCST007952') == (
        '2017_Wheeler_PLoSMed_HbA1c_AFR_GCST007952_GRCh38_formatted.tsv.gz'
    )


FILE = {'file_id': 'F', 'biomarker': 'B', 'ancestry': 'EUR', 'source_build': 'GRCh37'}


@pytest.mark.parametrize(
    'per_chromosome, changed, has_x',
    [
        ({'1': 9, '23': 2}, 2, False),  # no more X rows than liftover moved: strays
        ({'1': 9, '23': 3}, 2, True),  # one more than liftover could explain
        ({'1': 9}, 2, False),  # no X at all
        ({'1': 9, '23': 1}, 0, True),  # GRCh38 source, not lifted: any X is its own
    ],
)
def test_has_x_needs_more_x_rows_than_liftover_moved(per_chromosome, changed, has_x):
    record = {'rows_per_chromosome': per_chromosome}
    stats = {'liftover_chromosome_changed': changed} if changed else {}
    result = extract.chromosomes_record(FILE, record, stats)
    assert result['has_x'] is has_x
    assert result['x_rows'] == per_chromosome.get('23', 0)
    assert result['liftover_chromosome_changed'] == changed


def test_chromosome_table_from_hits_and_stats(monkeypatch, tmp_path):
    """
    End to end through main, nothing mocked: hits written by extract_one, stats as
    format.py writes them, and the table read back. One file with a published X,
    one with only liftover strays on X, and one not extracted yet.
    """
    formatted = tmp_path / 'v1'
    hits = formatted / 'hits' / 'p_1e-05'
    stats = formatted / 'stats'
    hits.mkdir(parents=True)
    stats.mkdir()
    published_x = [row(chromosome='23', base_pair_location=i + 1) for i in range(3)]
    strays = [row(chromosome='23', base_pair_location=i + 1) for i in range(2)]
    for file_id, rows, changed in [
        ('WITH_X', [hit('rs1'), *published_x], 2),
        ('STRAYS', [hit('rs1'), *strays], 2),
    ]:
        source = tmp_path / f'{file_id}.tsv.gz'
        write_file(source, rows)
        extract.extract_one(source, hits / f'{file_id}.parquet', file_id, 1e-5)
        (stats / f'{file_id}.json').write_text(
            json.dumps({'liftover_chromosome_changed': changed})
        )
    files_csv = tmp_path / 'files.csv'
    files_csv.write_text(
        'file_id,biomarker,ancestry,source_build\n'
        'WITH_X,HbA1c,AFR,GRCh37\nSTRAYS,HbA1c,EUR,GRCh37\nLATER,HbA1c,EAS,GRCh38\n'
    )
    out = tmp_path / 'table.csv'
    monkeypatch.setattr(sys, 'argv', [
        'extract.py', '--chromosomes', '--files', str(files_csv),
        '--formatted', str(formatted), '--chromosomes-out', str(out),
    ])
    extract.main()

    table = pd.read_csv(out, dtype=str, keep_default_na=False).set_index('file_id')
    assert table.loc['WITH_X', ['extracted', 'has_x', 'x_rows']].tolist() == [
        'True', 'True', '3',
    ]
    assert table.loc['STRAYS', ['extracted', 'has_x', 'x_rows']].tolist() == [
        'True', 'False', '2',
    ]
    assert table.loc['WITH_X', 'chromosomes'] == '1 23'
    assert table.loc['LATER', ['extracted', 'has_x', 'chromosomes']].tolist() == [
        'False', '', '',
    ]


def test_chromosome_table_needs_stats_for_every_extracted_file(monkeypatch, tmp_path):
    formatted = tmp_path / 'v1'
    hits = formatted / 'hits' / 'p_1e-05'
    hits.mkdir(parents=True)
    source = tmp_path / 'A.tsv.gz'
    write_file(source, [hit('rs1')])
    extract.extract_one(source, hits / 'A.parquet', 'A', 1e-5)
    files_csv = tmp_path / 'files.csv'
    files_csv.write_text('file_id,biomarker,ancestry,source_build\nA,HbA1c,AFR,GRCh37\n')
    with pytest.raises(ValueError, match='A: no format.py stats'):
        extract.run_chromosomes(
            [{'file_id': 'A'}], str(tmp_path / 'out.csv'), str(hits), str(tmp_path)
        )


def test_hit_types_cover_format_py_columns_in_order():
    """HIT_TYPES is declared here because each Batch job runs extract.py alone; it
    must still list exactly the columns format.py writes."""
    assert list(extract.HIT_TYPES) == ['file_id', *fmt.OUTPUT_COLUMNS]


def test_hits_are_typed_so_numbers_compare_as_numbers(tmp_path):
    """As text, '10.2' > '7' is False and a strong hit drops out of a filter."""
    source, dest = tmp_path / 'in.tsv.gz', tmp_path / 'out.parquet'
    write_file(
        source,
        [
            row(
                rsid='rs1',
                chromosome='23',
                neg_log_10_p_value='10.2',
                p_value='6.3e-11',
                n='5000',
                z='6.6',
                effect_allele_is_alt='false',
            ),
            row(rsid='rs2', neg_log_10_p_value='7.5', p_value='3.2e-8', n='NA'),
        ],
    )
    extract.extract_one(source, dest, 'A', 1e-5)

    schema = pq.read_schema(dest)
    for name, kind in extract.HIT_TYPES.items():
        assert schema.field(name).type == pa.type_for_alias(kind), name
    hits = pd.read_parquet(dest)
    assert hits.loc[hits['neg_log_10_p_value'] > 7, 'rsid'].tolist() == ['rs1', 'rs2']
    assert hits['chromosome'].tolist() == [23, 1]
    assert hits['p_value'].tolist() == ['6.3e-11', '3.2e-8']
    assert hits['effect_allele_is_alt'].tolist() == [False, True]
    assert hits['n'].iloc[0] == 5000 and pd.isna(hits['n'].iloc[1])


def test_a_value_that_does_not_fit_its_type_fails_loudly(tmp_path):
    source = tmp_path / 'in.tsv.gz'
    write_file(source, [hit('rs1', effect_allele_is_alt='maybe')])
    with pytest.raises(ValueError, match='effect_allele_is_alt: not true/false'):
        extract.extract_one(source, tmp_path / 'out.parquet', 'A', 1e-5)


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


@pytest.mark.parametrize(
    'out, threshold, expected',
    [
        ('gs://b/v1/hits', 1e-5, 'gs://b/v1/hits/p_1e-05'),
        ('gs://b/v1/hits/', 1e-5, 'gs://b/v1/hits/p_1e-05'),  # trailing slash
        ('gs://b/v1/hits', 5e-8, 'gs://b/v1/hits/p_5e-08'),
        ('hits', 0.001, 'hits/p_0.001'),
    ],
)
def test_each_threshold_gets_its_own_folder(out, threshold, expected):
    assert extract.hits_folder(out, threshold) == expected


def test_runs_at_two_thresholds_never_share_a_folder(monkeypatch, tmp_path):
    """The case the old --out guard missed: a custom --out, with or without a
    trailing slash, rerun at another threshold. Each run gets its own subfolder."""
    files_csv = tmp_path / 'files.csv'
    files_csv.write_text('file_id,biomarker\nA,HbA1c\n')
    targets = []
    monkeypatch.setattr(extract, 'run_batch', lambda rows, args, out: targets.append(out))
    for out, threshold in [('gs://b/custom', '5e-8'), ('gs://b/custom/', '1e-5')]:
        monkeypatch.setattr(sys, 'argv', [
            'extract.py', '--files', str(files_csv), '--out', out,
            '--threshold', threshold,
        ])
        extract.main()

    assert targets == ['gs://b/custom/p_5e-08', 'gs://b/custom/p_1e-05']


def test_one_file_runs_on_a_laptop_from_the_command_line(monkeypatch, tmp_path):
    """The single-file command, which every Batch job runs, needs only local paths."""
    source, dest = tmp_path / 'in.tsv.gz', tmp_path / 'A.parquet'
    write_file(source, [hit('rs1'), row(rsid='rs2')])
    monkeypatch.setattr(sys, 'argv', [
        'extract.py', '--file-id', 'A', '--source', str(source), '--dest', str(dest),
    ])
    extract.main()
    assert pd.read_parquet(dest)['rsid'].tolist() == ['rs1']


def test_every_path_defaults_to_a_place_in_formatted(monkeypatch, tmp_path):
    """The chromosome table, the hits and the stats all follow --formatted, so a
    local run reads and writes locally rather than over gs://cpg-common-main."""
    called = {}
    monkeypatch.setattr(
        extract,
        'run_chromosomes',
        lambda rows, out, hits, stats: called.update(out=out, hits=hits, stats=stats),
    )
    monkeypatch.setattr(sys, 'argv', [
        'extract.py', '--chromosomes', '--formatted', f'{tmp_path}/',
    ])
    extract.main()
    assert called == {
        'out': f'{tmp_path}/chromosomes_reported.csv',
        'hits': f'{tmp_path}/hits/p_1e-05',
        'stats': f'{tmp_path}/stats',
    }


@pytest.mark.parametrize(
    'source_bytes, expected',
    [(1, '10Gi'), (2**30, '13Gi'), (int(6.4 * 2**30), '29Gi')],
)
def test_job_disk_follows_the_formatted_file(source_bytes, expected):
    assert extract.job_storage(source_bytes, 'A') == expected


@pytest.mark.parametrize('source_bytes', [0, -1])
def test_job_disk_refuses_an_unreadable_size(source_bytes):
    with pytest.raises(ValueError, match='A: formatted file has size'):
        extract.job_storage(source_bytes, 'A')


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
        (['--chromosomes', '--only', 'X'], '--chromosomes-out'),
        (['--chromosomes', '--files', 'subset.csv'], '--chromosomes-out'),
        (['--threshold', '5'], 'between 0 and 1'),
        (['--threshold', '1.5'], 'between 0 and 1'),
        (['--threshold', '0'], 'between 0 and 1'),
        (['--threshold', 'nan'], 'between 0 and 1'),
        (['--file-id', 'A', '--threshold', '5'], 'between 0 and 1'),
    ],
)
def test_main_refuses_unsafe_combinations(monkeypatch, capsys, argv, message):
    monkeypatch.setattr(sys, 'argv', ['extract.py', *argv])
    with pytest.raises(SystemExit):
        extract.main()
    assert message in capsys.readouterr().err
