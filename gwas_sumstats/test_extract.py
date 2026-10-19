"""
Tests for extract.py: the P value filter and the per-file extraction.

From the repo root:

    uv run --no-project --with pytest --with pandas --with pyarrow \
        pytest gwas_sumstats/test_extract.py
"""

import gzip

import extract
import pandas as pd
import pytest

HEADER = (
    'chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t'
    'standard_error\teffect_allele_frequency\tp_value\tneg_log_10_p_value\trsid\tn'
)


def write_file(path, rows: list[str]) -> None:
    """A minimal formatted file: the standard header plus the given rows."""
    with gzip.open(path, 'wt') as handle:
        handle.write(HEADER + '\n')
        for row in rows:
            handle.write(row + '\n')


@pytest.mark.parametrize(
    'p_text, neg_log_10, expected',
    [
        ('1e-9', '9', True),  # well below the cutoff
        ('1e-5', '5', True),  # exactly at it
        ('0.01', '2', False),  # above it
        ('8.51e-3553', '3552.07', True),  # below the float range: the strongest hits
        ('0', 'NA', True),  # underflowed to zero: stronger than any float, so kept
        ('NA', 'NA', False),  # no P value at all
        ('1e-9', '', True),  # no neg_log column: falls back to the text
        ('1e-9', 'NA', True),  # neg_log missing for this row
        ('not a number', 'NA', False),  # unparseable
    ],
)
def test_threshold(p_text, neg_log_10, expected):
    assert extract.passes(p_text, neg_log_10, 1e-5) is expected


def test_zero_p_value_is_rescued_by_neg_log_10():
    """
    A P value that underflowed to zero in the source is still usable when the file
    carries neg_log_10_p_value, and those rows are the strongest associations there
    are: both UGT1A1 bilirubin variants reach zero this way.
    """
    assert extract.passes('0', '3552.07', 1e-5) is True


def test_extract_keeps_only_hits_and_counts_rows(tmp_path):
    source = tmp_path / 'in.tsv.gz'
    write_file(
        source,
        [
            '1\t100\tA\tG\t0.5\t0.01\t0.3\t1e-20\t20\trs1\t1000',
            '1\t200\tC\tT\t0.1\t0.01\t0.3\t0.5\t0.30\trs2\t1000',
            '2\t300\tG\tA\t0.4\t0.01\t0.3\t1e-6\t6\trs3\t1000',
        ],
    )
    dest = tmp_path / 'out.parquet'
    counts = extract.extract_one(source, dest, 'STUDY_1', 1e-5)

    assert counts == {'file_id': 'STUDY_1', 'rows_in': 3, 'rows_kept': 2}
    hits = pd.read_parquet(dest)
    assert list(hits['rsid']) == ['rs1', 'rs3']
    assert set(hits['file_id']) == {'STUDY_1'}


def test_extract_handles_a_file_with_no_hits(tmp_path):
    """An empty result still gets a parquet with the right columns, so that reading
    the whole folder back never has to special-case it."""
    source = tmp_path / 'in.tsv.gz'
    write_file(source, ['1\t100\tA\tG\t0.01\t0.01\t0.3\t0.9\t0.05\trs1\t1000'])
    dest = tmp_path / 'out.parquet'

    counts = extract.extract_one(source, dest, 'STUDY_2', 1e-5)
    assert counts['rows_kept'] == 0
    empty = pd.read_parquet(dest)
    assert len(empty) == 0
    assert 'p_value' in empty.columns and 'file_id' in empty.columns


def test_extract_reads_across_chunks(tmp_path, monkeypatch):
    """Chunking must not lose or duplicate rows, whatever the chunk boundary."""
    monkeypatch.setattr(extract, 'CHUNK_ROWS', 2)
    source = tmp_path / 'in.tsv.gz'
    write_file(
        source,
        [
            f'1\t{100 + i}\tA\tG\t0.5\t0.01\t0.3\t'
            f'{"1e-20" if i % 2 else "0.5"}\t{20 if i % 2 else 0.3}\trs{i}\t1000'
            for i in range(10)
        ],
    )
    dest = tmp_path / 'out.parquet'

    counts = extract.extract_one(source, dest, 'STUDY_3', 1e-5)
    assert counts == {'file_id': 'STUDY_3', 'rows_in': 10, 'rows_kept': 5}


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

    assert len(extract.read_files_csv(files_csv)) == 2
    assert [r['file_id'] for r in extract.read_files_csv(files_csv, ['B'])] == ['B']
    with pytest.raises(ValueError, match='not in files.csv'):
        extract.read_files_csv(files_csv, ['C'])
