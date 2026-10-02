"""
Tests for the format.py transforms, on small made-up inputs (no downloads, no cloud).

From the repo root:

    uv run --no-project --with pytest --with polars==1.34.0 --with pysam==0.23.3 \
        --with pyliftover==0.4.1 --with numpy pytest gwas_sumstats
"""

import math

import polars as pl
import pysam
import pytest

import format as fmt

GENOME_LENGTH = 1000


def source_table(**columns) -> pl.DataFrame:
    """A source file as read_table returns it: every column text."""
    return pl.DataFrame({name: [str(v) for v in values] for name, values in columns.items()})


def standardise(df: pl.DataFrame, n_study: int = 100):
    counts: dict = {}
    table, notes = fmt.standardise(df, fmt.resolve_columns(df, ''), n_study, counts)
    return table, notes, counts


@pytest.fixture
def fasta(tmp_path):
    """
    GRCh38 stand-in with Ensembl contig names 1 and 2, all 'C' except chosen bases.
    1-based positions: 1:50 = A, 1:251 = A, 2:990 = G.
    """
    sequences = {'1': ['C'] * GENOME_LENGTH, '2': ['C'] * GENOME_LENGTH}
    sequences['1'][50 - 1] = 'A'
    sequences['1'][251 - 1] = 'A'
    sequences['2'][990 - 1] = 'G'
    path = tmp_path / 'genome.fa'
    path.write_text(
        ''.join(f'>{name}\n{"".join(bases)}\n' for name, bases in sequences.items())
    )
    pysam.faidx(str(path))
    return path


@pytest.fixture
def chain(tmp_path):
    """
    Source chr1 0-100 -> target chr1 200-300 (plus strand, +200).
    Source chr2 0-100 -> target chr2 900-1000 (minus strand: 1-based p -> 1001 - p).
    """
    path = tmp_path / 'test.over.chain'
    path.write_text(
        f'chain 1000 chr1 {GENOME_LENGTH} + 0 100 chr1 {GENOME_LENGTH} + 200 300 1\n'
        '100\n\n'
        f'chain 1000 chr2 {GENOME_LENGTH} + 0 100 chr2 {GENOME_LENGTH} - 0 100 2\n'
        '100\n\n'
    )
    return path


def test_odds_ratio_becomes_ln_or_with_se_from_ci():
    df = source_table(
        chromosome=[1],
        base_pair_location=[50],
        effect_allele=['A'],
        other_allele=['G'],
        odds_ratio=[2],
        ci_lower=[1.5],
        ci_upper=[2.6],
        p_value=['0.01'],
    )
    table, notes, _ = standardise(df)
    assert table['beta'][0] == pytest.approx(math.log(2))
    assert table['standard_error'][0] == pytest.approx(
        (math.log(2.6) - math.log(1.5)) / (2 * 1.959963984540054)
    )
    assert notes['effect_source'] == 'ln(odds_ratio), SE from 95% CI'


def test_p_value_below_float_range_keeps_text():
    df = source_table(
        chromosome=[1],
        base_pair_location=[50],
        effect_allele=['A'],
        other_allele=['G'],
        beta=[0.1],
        p_value=['8.51e-3553'],
    )
    table, _, _ = standardise(df)
    assert table['p_value'][0] == '8.51e-3553'
    assert table['neg_log_10_p_value'][0] == pytest.approx(3552.070070, abs=1e-6)


def test_p_value_rebuilt_from_neg_log10():
    df = source_table(
        chromosome=[1],
        base_pair_location=[50],
        effect_allele=['A'],
        other_allele=['G'],
        beta=[0.1],
        neg_log_10_p_value=[5],
    )
    table, notes, _ = standardise(df)
    assert table['p_value'][0] == '1.0e-5'
    assert notes['p_value_source'] == 'from neg_log_10_p_value'


def test_negative_neg_log10_fails():
    df = source_table(
        chromosome=[1],
        base_pair_location=[50],
        effect_allele=['A'],
        other_allele=['G'],
        beta=[0.1],
        neg_log_10_p_value=[-5],
    )
    with pytest.raises(ValueError, match='log10'):
        standardise(df)


def test_sequence_free_indels_are_dropped_and_counted():
    df = source_table(
        chromosome=[1, 1],
        base_pair_location=[50, 60],
        effect_allele=['A', 'D'],
        other_allele=['G', 'I'],
        beta=[0.1, 0.2],
        p_value=['0.01', '0.02'],
    )
    table, _, counts = standardise(df)
    assert table.height == 1
    assert counts['dropped_non_acgt_allele'] == 1


def test_hapmap_frequency_is_not_used_as_eaf():
    df = source_table(
        snp=['rs1'],
        chr=[1],
        pos=[50],
        effect_allele=['A'],
        other_allele=['G'],
        eaf_hapmap_YRI=[0.3],
        beta=[0.1],
        pvalue=['0.01'],
    )
    assert 'effect_allele_frequency' not in fmt.resolve_columns(df, '')
    table, _, _ = standardise(df)
    assert table['effect_allele_frequency'][0] is None


def test_other_strand_snv_is_complemented_and_keeps_beta(fasta):
    # GRCh38 base at 1:50 is A. C/T matches only as its complement G/A.
    df = source_table(
        chromosome=[1, 1],
        base_pair_location=[50, 60],
        effect_allele=['C', 'A'],
        other_allele=['T', 'T'],
        beta=[0.1, 0.2],
        p_value=['0.01', '0.02'],
    )
    table, _, _ = standardise(df)
    counts: dict = {}
    checked = fmt.check_reference(table, fasta, counts)
    flipped = checked.filter(pl.col('base_pair_location') == 50).row(0, named=True)
    assert (flipped['effect_allele'], flipped['other_allele']) == ('G', 'A')
    assert flipped['beta'] == pytest.approx(0.1)
    # 1:60 is C, so A/T matches neither strand and is dropped.
    assert checked.height == 1
    assert counts['reference_strand_flipped'] == 1
    assert counts['reference_mismatch'] == 1


def test_liftover_plus_and_minus_strand(chain, monkeypatch):
    monkeypatch.setitem(fmt.LIFTOVER_KNOWN, 'GRCh37', [(1, 51, 251), (2, 11, 990)])
    df = source_table(
        chromosome=[1, 2, 2],
        base_pair_location=[51, 11, 20],
        effect_allele=['A', 'C', 'AT'],
        other_allele=['G', 'T', 'A'],
        beta=[0.1, 0.2, 0.3],
        p_value=['0.01', '0.02', '0.03'],
    )
    table, _, _ = standardise(df)
    counts: dict = {}
    lifted = fmt.lift_to_grch38(table, chain, 'GRCh37', counts)
    rows = {r['chromosome']: r for r in lifted.iter_rows(named=True)}
    assert rows[1]['base_pair_location'] == 251
    assert (rows[1]['effect_allele'], rows[1]['other_allele']) == ('A', 'G')
    assert rows[2]['base_pair_location'] == 990
    assert (rows[2]['effect_allele'], rows[2]['other_allele']) == ('G', 'A')
    assert rows[2]['beta'] == pytest.approx(0.2)
    assert counts['liftover_minus_strand_complemented'] == 1
    assert counts['liftover_minus_strand_multibase_dropped'] == 1


def test_liftover_offset_check_catches_off_by_one(chain, monkeypatch):
    monkeypatch.setitem(fmt.LIFTOVER_KNOWN, 'GRCh37', [(1, 51, 252)])
    table, _, _ = standardise(
        source_table(
            chromosome=[1],
            base_pair_location=[51],
            effect_allele=['A'],
            other_allele=['G'],
            beta=[0.1],
            p_value=['0.01'],
        )
    )
    with pytest.raises(AssertionError, match='offset check failed'):
        fmt.lift_to_grch38(table, chain, 'GRCh37', {})


def test_every_row_dropped_gives_clear_error(fasta):
    # Quoted CSV values keep their quotes, so no chromosome is recognised.
    df = source_table(
        chromosome=['"1"'],
        base_pair_location=[50],
        effect_allele=['A'],
        other_allele=['G'],
        beta=[0.1],
        p_value=['0.01'],
    )
    table, _, counts = standardise(df)
    with pytest.raises(ValueError, match='dropped_non_canonical_chromosome'):
        fmt.check_reference(table, fasta, counts)


def test_mt_in_lifted_file_gets_its_own_drop_label(chain, monkeypatch):
    monkeypatch.setitem(fmt.LIFTOVER_KNOWN, 'GRCh37', [(1, 51, 251)])
    table, _, _ = standardise(
        source_table(
            chromosome=[1, 'MT'],
            base_pair_location=[51, 100],
            effect_allele=['A', 'A'],
            other_allele=['G', 'G'],
            beta=[0.1, 0.2],
            p_value=['0.01', '0.02'],
        )
    )
    counts: dict = {}
    lifted = fmt.lift_to_grch38(table, chain, 'GRCh37', counts)
    assert lifted.height == 1
    assert counts['liftover_mt_not_lifted'] == 1
    assert counts['liftover_unmapped'] == 0
