"""
Tests for the format.py transforms, on small made-up inputs (no downloads, no cloud).

From the repo root:

    uv run --no-project --with pytest --with polars==1.34.0 --with pysam==0.23.3 \
        --with pyliftover==0.4.1 --with numpy pytest gwas_sumstats
"""

import argparse
import csv
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


def test_tool_specific_allele_names_need_an_override():
    # METAL's Allele1 is the effect allele, SAIGE's is not, so neither is guessed.
    df = source_table(
        chromosome=[1],
        base_pair_location=[50],
        Allele1=['A'],
        Allele2=['G'],
        Freq1=[0.3],
        beta=[0.1],
        p_value=['0.01'],
    )
    with pytest.raises(ValueError, match='no alleles'):
        standardise(df)
    columns = fmt.resolve_columns(
        df, 'effect_allele=Allele2;other_allele=Allele1;effect_allele_frequency=Freq1'
    )
    table, _ = fmt.standardise(df, columns, 100, {})
    assert (table['effect_allele'][0], table['other_allele'][0]) == ('G', 'A')
    assert table['effect_allele_frequency'][0] == pytest.approx(0.3)


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
    with pytest.raises(RuntimeError, match='offset check failed'):
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


def test_manifest_refuses_only(monkeypatch, capsys):
    # A partial manifest would replace the full index.
    monkeypatch.setattr(
        'sys.argv', ['format.py', '--manifest', '--only', 'any_file_id']
    )
    with pytest.raises(SystemExit) as stop:
        fmt.main()
    assert stop.value.code == 2
    assert 'takes no --only' in capsys.readouterr().err


def test_manifest_fails_on_unformatted_files_unless_allowed(tmp_path):
    rows = [{'file_id': 'done'}, {'file_id': 'missing'}]
    (tmp_path / 'stats').mkdir()
    (tmp_path / 'stats' / 'done.json').write_text('{"file_id": "done", "rows_out": 1}')
    (tmp_path / 'done_GRCh38_formatted.tsv.gz').write_text('')
    args = argparse.Namespace(
        out=str(tmp_path),
        stats=str(tmp_path / 'stats'),
        web=str(tmp_path / 'web'),
        allow_missing=False,
    )
    with pytest.raises(ValueError, match='not formatted: missing'):
        fmt.write_manifest(rows, args)
    assert not (tmp_path / 'manifest.csv').exists()

    args.allow_missing = True
    fmt.write_manifest(rows, args)
    with (tmp_path / 'manifest.csv').open() as handle:
        assert [line['file_id'] for line in csv.DictReader(handle)] == ['done']
    page = (tmp_path / 'web' / 'manifest.html').read_text()
    assert 'not formatted:</b> missing' in page


def test_scripts_have_no_assert_statements():
    # python -O strips assert, so a data check written as one would stop gating.
    import ast
    from pathlib import Path

    for name in ('download.py', 'format.py'):
        tree = ast.parse(Path(fmt.__file__).with_name(name).read_text())
        asserts = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)]
        assert not asserts, f'{name} uses assert at lines {asserts}'


@pytest.mark.parametrize(
    'markers',
    [
        ['3:100796605:A:G', '1:174816432:C:T'],  # Moksnes 2021 troponin, real lines
        ['Chr3:100796605', 'CHR1:174816432'],
        ['chrx_100796605', 'chr1_174816432'],
    ],
)
def test_position_parsed_from_marker_in_any_case(markers):
    df = source_table(
        MarkerName=markers,
        effect_allele=['A', 'C'],
        other_allele=['G', 'T'],
        beta=[0.1, 0.2],
        p_value=['0.01', '0.02'],
    )
    table, notes, counts = standardise(df)
    assert notes['position_source'] == 'parsed from MarkerName'
    assert counts['dropped_non_canonical_chromosome'] == 0
    expected_first = 23 if markers[0].lower().startswith('chrx') else 3
    assert table['chromosome'].to_list() == [expected_first, 1]
    assert table['base_pair_location'].to_list() == [100796605, 174816432]


def test_format_one_works_chromosome_by_chromosome(fasta, tmp_path):
    # Interleaved chromosomes, one unrecognised, whitespace-separated and gzipped.
    import gzip
    import json

    source = tmp_path / 'source.txt.gz'
    with gzip.open(source, 'wt') as handle:
        handle.write(
            'chromosome base_pair_location effect_allele other_allele beta p_value\n'
        )
        handle.write('2 990 G A 0.3 0.03\n')
        handle.write('Un 5 C T 0.9 0.09\n')
        handle.write('1 251 A G 0.2 0.02\n')
        handle.write('1 50 A C 0.1 0.01\n')
    row = {
        'file_id': 'test',
        'columns': '',
        'n_study': '100',
        'source_build': 'GRCh38',
    }
    out_root = tmp_path / 'out'
    stats_path = tmp_path / 'stats.json'
    record = {'md5': 'x', 'md5_checked': False, 'downloaded_at': 'x'}
    fmt.format_one(row, source, out_root, stats_path, fasta, {}, record)

    lines = gzip.open(f'{out_root}.tsv.gz', 'rt').read().splitlines()
    assert lines[0].split('\t') == fmt.OUTPUT_COLUMNS
    assert [line.split('\t')[:2] for line in lines[1:]] == [
        ['1', '50'],
        ['1', '251'],
        ['2', '990'],
    ]
    stats = json.loads(stats_path.read_text())
    assert (stats['rows_in'], stats['rows_out']) == (4, 3)
    assert stats['dropped_non_canonical_chromosome'] == 1
    assert stats['reference_build_checked'] is False
    assert pysam.TabixFile(f'{out_root}.tsv.gz').contigs == ['1', '2']


def test_every_row_failing_grch38_writes_nothing(fasta, tmp_path):
    # Rows reach the reference check, all mismatch, and too few are palindromic for
    # the build check to judge: the file must still fail, with no output written.
    source = tmp_path / 'source.tsv'
    source.write_text(
        'chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tp_value\n'
        '1\t10\tA\tT\t0.1\t0.01\n'
        '1\t20\tA\tT\t0.2\t0.02\n'
    )
    row = {'file_id': 'x', 'columns': '', 'n_study': '100', 'source_build': 'GRCh38'}
    out_root = tmp_path / 'out'
    record = {'md5': 'x', 'md5_checked': False, 'downloaded_at': 'x'}
    stats = tmp_path / 'stats.json'
    with pytest.raises(ValueError, match='no rows survived'):
        fmt.format_one(row, source, out_root, stats, fasta, {}, record)
    assert not list(tmp_path.glob('out*'))
    assert not stats.exists()


def test_unparseable_p_values_become_na_and_are_counted():
    df = source_table(
        chromosome=[1, 1, 1, 1],
        base_pair_location=[10, 20, 30, 40],
        effect_allele=['A'] * 4,
        other_allele=['G'] * 4,
        beta=[0.1] * 4,
        p_value=['<1e-300', '1e-5', '-0.01', '1.05e-1046'],
    )
    table, notes, counts = standardise(df)
    assert table['p_value'].to_list() == [None, '1e-5', None, '1.05e-1046']
    assert table['neg_log_10_p_value'][1] == pytest.approx(5.0)
    assert counts['p_value_unparseable'] == 2
    assert counts['p_value_zero'] == 0
    assert sorted(notes['p_value_unparseable_examples']) == ['-0.01', '<1e-300']


def test_many_unparseable_p_values_fail_the_file(fasta, tmp_path):
    source = tmp_path / 'source.tsv'
    source.write_text(
        'chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tp_value\n'
        '1\t50\tA\tC\t0.1\t<0.001\n'
        '1\t251\tA\tG\t0.2\t0.02\n'
    )
    row = {'file_id': 'x', 'columns': '', 'n_study': '100', 'source_build': 'GRCh38'}
    record = {'md5': 'x', 'md5_checked': False, 'downloaded_at': 'x'}
    with pytest.raises(ValueError, match=r'not numbers in \(0, 1\].*<0.001'):
        fmt.format_one(
            row, source, tmp_path / 'out', tmp_path / 'stats.json', fasta, {}, record
        )


def test_single_space_file_keeps_empty_fields_in_place(tmp_path):
    # Teumer 2019 UACR: single spaces, and a missing rsID is an empty field.
    source = tmp_path / 'source.txt'
    source.write_text(
        'Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum\n'
        '1 729679 rs4951859 c g 0.1585 -0.0059704 0.0030381 0.04939 491125\n'
        '1 2556125  t c 0.3260 -0.0035140 0.0021112 0.09603 547349\n'
    )
    assert fmt.detect_separator(source) == ' '
    df = fmt.read_table(source, tmp_path)
    assert df['RSID'].to_list() == ['rs4951859', None]
    assert df['P-value'].to_list() == ['0.04939', '0.09603']
    parts = fmt.split_by_chromosome(
        source, tmp_path, {'chromosome': 'Chr', 'base_pair_location': 'Pos_b37'}
    )
    lines = parts[1].read_text().splitlines()
    assert lines[2].split('\t')[2:5] == ['', 't', 'c']


def test_space_padded_columns_still_split_on_runs(tmp_path):
    source = tmp_path / 'source.txt'
    source.write_text('CHR   POS  A1 A2\n1     100  A  G\n22  20000  C  T\n')
    assert fmt.detect_separator(source) is None
    assert fmt.read_table(source, tmp_path)['POS'].to_list() == ['100', '20000']


def test_zero_p_value_is_na_counted_apart_and_keeps_beta():
    # Software writes an underflowed P (APOE in Timsina 2026) as 0.0.
    df = source_table(
        chromosome=[19, 19],
        base_pair_location=[10, 20],
        effect_allele=['A', 'A'],
        other_allele=['G', 'G'],
        beta=[-0.4774, 0.1],
        p_value=['0.0', '0.5'],
    )
    table, _, counts = standardise(df)
    assert table['p_value'].to_list() == [None, '0.5']
    assert table['beta'][0] == pytest.approx(-0.4774)
    assert (counts['p_value_zero'], counts['p_value_unparseable']) == (1, 0)


def test_build_check_is_recorded_as_run_or_skipped():
    enough = {'reference_ok': 500, 'reference_palindromic': 200}
    fmt.check_build(enough)
    assert enough['reference_build_checked'] is True
    too_few = {'reference_ok': 500, 'reference_palindromic': 50}
    fmt.check_build(too_few)
    assert too_few['reference_build_checked'] is False


def test_infinite_neg_log10_is_an_underflowed_p_not_a_crash():
    df = source_table(
        chromosome=[1, 1, 1],
        base_pair_location=[10, 20, 30],
        effect_allele=['A'] * 3,
        other_allele=['G'] * 3,
        beta=[0.5, 0.1, 0.2],
        neg_log_10_p_value=['Infinity', '5', 'inf'],
    )
    table, _, counts = standardise(df)
    assert table['p_value'].to_list() == [None, '1.0e-5', None]
    assert table['neg_log_10_p_value'].to_list() == [None, 5.0, None]
    assert (counts['p_value_zero'], counts['p_value_unparseable']) == (2, 0)


def test_infinite_beta_becomes_na():
    df = source_table(
        chromosome=[1],
        base_pair_location=[10],
        effect_allele=['A'],
        other_allele=['G'],
        beta=['Infinity'],
        p_value=['0.5'],
    )
    table, _, _ = standardise(df)
    assert table['beta'][0] is None


def test_rsid_column_found_when_its_first_rows_are_empty():
    n = 1500
    df = source_table(
        chromosome=[1] * n,
        base_pair_location=list(range(1, n + 1)),
        effect_allele=['A'] * n,
        other_allele=['G'] * n,
        beta=[0.1] * n,
        p_value=['0.5'] * n,
        SNP=[''] * 600 + ['.'] * 400 + [f'rs{i}' for i in range(500)],
    )
    assert fmt.resolve_columns(df, '')['rsid'] == 'SNP'


def test_pins_match_everywhere():
    # Batch jobs install PIP_PACKAGES; the README and this module's docstring must
    # tell people to test and run locally with the same versions.
    from pathlib import Path

    readme = Path(fmt.__file__).with_name('README.md').read_text()
    for pin in fmt.PIP_PACKAGES.split():
        assert pin in __doc__, f'{pin} missing from the test_format.py docstring'
        assert readme.count(pin) >= 2, f'{pin} missing from a README command'


def test_se_from_bad_ci_is_na_not_inf_or_nan():
    df = source_table(
        chromosome=[1, 1, 1, 1],
        base_pair_location=[10, 20, 30, 40],
        effect_allele=['A'] * 4,
        other_allele=['G'] * 4,
        odds_ratio=[2, 2, 2, 2],
        ci_lower=[0, -1, 'NA', 1.5],
        ci_upper=[2.6, 2.6, 2.6, 2.6],
        p_value=['0.1'] * 4,
    )
    table, _, _ = standardise(df)
    se = table['standard_error'].to_list()
    assert se[:3] == [None, None, None]
    expected = (math.log(2.6) - math.log(1.5)) / (2 * 1.959963984540054)
    assert se[3] == pytest.approx(expected)
    assert table['standard_error'].is_nan().sum() == 0


def test_offset_check_catches_one_based_positions_on_minus_strand(chain, monkeypatch):
    # The harmoniser bug: 1-based positions passed to the 0-based chain, result used
    # as is. It cancels out on plus-strand blocks, so only a minus-strand known
    # variant catches it. Fixture chr2 is minus strand: 1-based p -> 1001 - p.
    monkeypatch.setitem(fmt.LIFTOVER_KNOWN, 'GRCh37', [(1, 51, 251), (2, 10, 991)])

    def one_based_bug(lifter, contig, position):
        hits = lifter.convert_coordinate(contig, position) or []
        return [(hit[0], hit[1], hit[2]) for hit in hits]

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
    fmt.lift_to_grch38(table, chain, 'GRCh37', {})  # correct code passes
    monkeypatch.setattr(fmt, 'lift_position', one_based_bug)
    with pytest.raises(RuntimeError, match=r'chr2:10 .* expected 991'):
        fmt.lift_to_grch38(table, chain, 'GRCh37', {})


def test_beta_from_z_and_standard_error():
    df = source_table(
        chromosome=[1, 1],
        base_pair_location=[10, 20],
        effect_allele=['A', 'A'],
        other_allele=['G', 'G'],
        z=[2.5, -4],
        standard_error=[0.02, 0.01],
        p_value=['0.01', '0.0001'],
    )
    table, notes, _ = standardise(df)
    assert table['beta'].to_list() == pytest.approx([0.05, -0.04])
    assert notes['effect_source'] == 'z x standard_error'


def test_build_check_falls_back_to_indels():
    # Too few A/T and C/G SNPs: judged on indels instead.
    right = {
        'reference_ok': 900,
        'reference_palindromic': 0,
        'reference_multibase': 500,
        'reference_multibase_mismatch': 2,
    }
    fmt.check_build(right)
    assert right['reference_build_check'] == 'indels'
    assert right['reference_build_checked'] is True
    wrong = dict(right, reference_multibase_mismatch=380)
    with pytest.raises(ValueError, match=r'76.0% of indels do not match GRCh38'):
        fmt.check_build(wrong)
    neither = {
        'reference_ok': 900,
        'reference_palindromic': 10,
        'reference_multibase': 20,
    }
    fmt.check_build(neither)
    assert neither['reference_build_checked'] is False
