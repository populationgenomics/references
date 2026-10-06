#!/usr/bin/env python3

"""
Convert the downloaded GWAS summary statistics to one GRCh38 layout.

For each file in files.csv: map its columns onto the GWAS-SSF names, lift GRCh37 or
NCBI36 positions to GRCh38, check both alleles against the GRCh38 reference, sort, and
write a bgzipped TSV with a tabix index. Nothing is filtered on frequency or P value;
rows are dropped only when they cannot be placed on GRCh38 (reasons are counted).

The largest sources hold ~200 million rows, more than a job's memory. So a file is
first split, streaming, into one file per chromosome on disk, and formatted one
chromosome at a time: peak memory follows the largest chromosome (~7 GB for 17
million rows), not the file.

Output columns, in order (missing values are NA):

    chromosome  base_pair_location  effect_allele  other_allele  beta
    standard_error  effect_allele_frequency  p_value  neg_log_10_p_value
    rsid  n  z

chromosome is GWAS-SSF numeric (1-22, X=23, Y=24, MT=25). base_pair_location is
1-based. beta is per effect allele: odds ratios become ln(OR), with the standard error
taken from the 95% CI if no SE is given; with only z and SE, beta = z x SE. A source
with no effect size (P value only) keeps beta as NA. z is filled only where the source has it.
p_value keeps the source's text, so values below the float range (e.g. 8.51e-3553)
survive; neg_log_10_p_value is computed from that text.

Liftover uses pyliftover, whose chain coordinates are 0-based: positions go in as
pos - 1 and come back + 1. check_liftover_offsets() checks this on variants with
known positions (LIFTOVER_KNOWN) before any file is lifted.

The build in files.csv is checked against the data: at the right positions both
alleles of an A/T or C/G SNP match GRCh38 (on one strand or the other), at wrong
positions about half do not. A file with too few of those is judged on its indels,
which must match GRCh38 as written. Either way, a file mismatching more than
MAX_BUILD_MISMATCH fails rather than being written (check_build).

On Hail Batch, one job per file, after download.py has finished:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Format biomarker GWAS summary statistics" \
        python3 gwas_sumstats/format.py

Each job also writes a stats JSON to v1/stats/ (row counts, drop reasons, and the
download record's MD5 and time), so the manifest needs nothing from the tmp bucket.
Once every job has finished, the manifest:

    analysis-runner --dataset common --access-level full \
        --output-dir gwas_sumstats \
        --description "Biomarker GWAS summary statistics manifest" \
        python3 gwas_sumstats/format.py --manifest

Before either, check that every source's columns resolve, from the first 200 kB of
each file at its source URL (no download, no cloud); writes columns_used.tsv, the
mapping per file, and exits non-zero if any file fails:

    python3 gwas_sumstats/format.py --check-columns

Locally (needs polars, pysam, pyliftover, numpy; where to get the fasta and chain
files: README.md, Running locally):

    python3 gwas_sumstats/format.py --local --originals ./original --out ./v1 \
        --stats ./stats --web ./web --fasta GRCh38.fa \
        --chain-grch37 hg19ToHg38.over.chain.gz \
        --chain-ncbi36 hg18ToHg38.over.chain.gz \
        --only 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951
"""

import argparse
import csv
import gzip
import html
import itertools
import json
import re
import shlex
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from download import FILES_CSV, fetch_head, original_name, read_files_csv

DEFAULT_ORIGINALS = 'gs://cpg-common-main-tmp/gwas_sumstats/original'
DEFAULT_STATS = 'gs://cpg-common-main/references/gwas_sumstats/v1/stats'
DEFAULT_OUT = 'gs://cpg-common-main/references/gwas_sumstats/v1'
DEFAULT_WEB = 'gs://cpg-common-main-web/gwas_sumstats'
COLUMNS_USED_TSV = Path(__file__).with_name('columns_used.tsv')
# The same pins appear in README.md (Running locally, Tests) and the test_format.py
# docstring; test_pins_match_everywhere fails if they drift apart.
PIP_PACKAGES = 'polars==1.34.0 pysam==0.23.3 pyliftover==0.4.1 numpy'

FLOAT_COLUMNS = [
    'beta',
    'standard_error',
    'effect_allele_frequency',
    'neg_log_10_p_value',
    'z',
]
OUTPUT_COLUMNS = [
    'chromosome',
    'base_pair_location',
    'effect_allele',
    'other_allele',
    'beta',
    'standard_error',
    'effect_allele_frequency',
    'p_value',
    'neg_log_10_p_value',
    'rsid',
    'n',
    'z',
]

# Source column names accepted for each field, matched case-insensitively, first
# hit wins. files.csv `columns` ('field=SourceName;...') overrides these per file.
ALIASES = {
    'chromosome': ['chromosome', 'chr', 'chrom', '#chrom', 'chr_num'],
    'base_pair_location': ['base_pair_location', 'pos', 'position', 'bp', 'pos_b37'],
    # Only the unambiguous GWAS-SSF names. Allele1/A1/ALT mean the effect allele
    # in some tools and the other allele in others (METAL: Allele1 is the effect
    # allele; SAIGE: Allele2 is), so a file using them needs a `columns` override
    # saying which is which, and its frequency column with it.
    'effect_allele': ['effect_allele'],
    'other_allele': ['other_allele'],
    'beta': ['beta', 'effect', 'b'],
    'standard_error': ['standard_error', 'se', 'stderr', 'sebeta'],
    'effect_allele_frequency': ['effect_allele_frequency', 'eaf'],
    'p_value': ['p_value', 'p', 'pvalue', 'p-value', 'pval', 'p.value'],
    'neg_log_10_p_value': ['neg_log_10_p_value', 'mlog10p', 'log10p'],
    'odds_ratio': ['odds_ratio', 'or'],
    'ci_lower': ['ci_lower'],
    'ci_upper': ['ci_upper'],
    'z': ['z', 'zscore', 'z_score'],
    'n': ['n', 'sample_size', 'n_total_sum', 'n_total'],
}
RSID_CANDIDATES = ['rsid', 'rs_id', 'rs_number', 'variant_id', 'snp', 'markername']
MARKER_CANDIDATES = ['markername', 'variant_id', 'variant', 'snp', 'id']
# (?i) in the pattern itself: looks_like (re) and the polars extract must agree.
MARKER_PATTERN = r'(?i)^(?:chr)?([0-9XYMT]+)[:_](\d+)'
# Infinities are not here: number() turns every non-finite value into NA, and a
# +inf -log10 P is counted as an underflowed P.
NULL_TOKENS = ['', 'NA', '#NA', '.', 'NAN', 'NULL', 'NONE']

# GWAS-SSF chromosome codes, and the contig names each needs in the chain file
# (UCSC hg19ToHg38) and in a GRCh38 fasta named either way (Ensembl '1', UCSC 'chr1').
# 25 = MT follows GWAS-SSF; PLINK numbers 25 = XY and 26 = MT.
CHROMOSOME_CODES = {str(i): i for i in range(1, 23)} | {
    'X': 23,
    'Y': 24,
    'M': 25,
    'MT': 25,
    '23': 23,
    '24': 24,
    '25': 25,
}
CHAIN_CONTIGS = {i: f'chr{i}' for i in range(1, 23)} | {23: 'chrX', 24: 'chrY'}
ENSEMBL_CONTIGS = {i: str(i) for i in range(1, 23)} | {23: 'X', 24: 'Y', 25: 'MT'}
UCSC_CONTIGS = {i: f'chr{i}' for i in range(1, 23)} | {
    23: 'chrX',
    24: 'chrY',
    25: 'chrM',
}
Z_95 = 1.959963984540054
# A P value is kept only as plain numeric text in (0, 1]. A P of exactly 0 is an
# underflow (the strongest hits, e.g. APOE in Timsina 2026): NA, counted as
# p_value_zero, beta and SE kept. Anything else ('<1e-300', '-0.01', '1.5') becomes
# NA as p_value_unparseable; more than MAX_UNPARSEABLE_P of those fails the file.
P_TEXT_PATTERN = r'^(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$'
MAX_UNPARSEABLE_P = 0.01
MAX_BUILD_MISMATCH = 0.10
MIN_SITES_FOR_BUILD_CHECK = 100

# Known positions that check_liftover_offsets lifts before any file, through the same
# lift_position() the data uses, failing the job if one lands elsewhere: a wrong chain
# file, or 1-based positions passed to the 0-based chain (the GWAS Catalog harmoniser
# bug, sumstats-harmoniser#52). That bug cancels out on plus-strand chain blocks and
# shifts positions by 2 on minus-strand ones (~0.45% of positions), so each build has
# a minus-strand variant. (chromosome, source position, GRCh38 position), 1-based;
# GRCh38 positions from Ensembl.
LIFTOVER_KNOWN = {
    'GRCh37': [
        (19, 45412079, 44908822),  # rs7412 (APOE), plus strand
        (19, 45411941, 44908684),  # rs429358 (APOE), plus strand
        (1, 144865454, 149019017),  # rs2798904, minus strand
    ],
    'NCBI36': [
        (2, 27584444, 27508073),  # rs1260326 (GCKR), plus strand
        (1, 711153, 785910),  # rs12565286, plus strand
        (1, 2476391, 2566588),  # rs2495365, minus strand
    ],
}


def formatted_name(row: dict) -> str:
    """File name of the formatted output, without the .tsv.gz / .tbi suffix."""
    return f'{row["file_id"]}_GRCh38_formatted'


def as_path(path: str):
    """A local Path, or a cloud path for gs:// (cpg_utils is only needed then)."""
    if path.startswith('gs://'):
        from cpg_utils import to_path

        return to_path(path)
    return Path(path)


def is_compressed(path: Path) -> bool:
    """True for gzip or bgzip."""
    with path.open('rb') as handle:
        return handle.read(2) == b'\x1f\x8b'


def open_text(path: Path):
    """Open a text file that may or may not be gzip/bgzip compressed."""
    if is_compressed(path):
        return gzip.open(path, 'rt', newline='')
    return path.open(newline='')


def detect_separator(path: Path) -> str | None:
    """
    '\t' or ',' from the header; ' ' if every sampled line splits on single spaces
    into as many fields as the header (an empty field, e.g. a missing rsID, is then
    two spaces and must not shift the row); None for space-padded columns, which
    are split on runs of whitespace.
    """
    with open_text(path) as handle:
        lines = (line.rstrip('\r\n') for line in handle if not line.startswith('##'))
        sample = list(itertools.islice(lines, 1001))
    header = sample[0]
    if '\t' in header:
        return '\t'
    if ',' in header:
        return ','
    names = header.split(' ')
    if '' not in names and all(len(line.split(' ')) == len(names) for line in sample):
        return ' '
    return None


def read_table(
    path: Path,
    workdir: Path,
    columns: list[str] | None = None,
    n_rows: int | None = None,
):
    """
    Read a summary statistics file with every column as text.

    polars would decompress a gzip input wholly in memory, so a compressed input is
    first streamed, decompressed, to `workdir`. Whitespace-separated files (runs of
    spaces) are rewritten as tab-separated on the way. Lines starting with '##' are
    skipped. format_one reads a 10,000-row sample through this, then the
    per-chromosome files from split_by_chromosome with only the columns it uses.

    Args:
        path: the original file
        workdir: scratch folder for the decompressed / rewritten copy
        columns: source columns to read; all if None
        n_rows: read only the first n_rows data rows (for a sample)
    """
    import polars as pl

    separator = detect_separator(path)
    if separator is None or is_compressed(path):
        plain = workdir / ('sample.tsv' if n_rows else 'plain.tsv')
        with open_text(path) as source, plain.open('w') as out:
            lines = (line for line in source if not line.startswith('##'))
            for count, line in enumerate(lines):
                if n_rows is not None and count > n_rows:
                    break
                out.write(line if separator else '\t'.join(line.split()) + '\n')
        path, separator = plain, separator or '\t'
    return pl.read_csv(
        path,
        separator=separator,
        infer_schema=False,
        quote_char=None,
        comment_prefix='##',
        columns=columns,
        n_rows=n_rows,
    )


def looks_like(column, pattern: str) -> bool:
    """
    True if most of the first 1000 non-missing values match the pattern. Missing
    values (null, '', NULL_TOKENS) are skipped first, so a column whose early rows
    are empty is still judged on the values it has.
    """
    present = column.drop_nulls()
    present = present.filter(
        ~present.str.strip_chars().str.to_uppercase().is_in(NULL_TOKENS)
    )
    sample = present.head(1000).to_list()
    hits = sum(bool(re.match(pattern, v, re.IGNORECASE)) for v in sample)
    return bool(sample) and hits / len(sample) > 0.9


def resolve_columns(df, override: str) -> dict[str, str]:
    """
    Map output fields to source column names.

    Args:
        df: the source table (all text)
        override: files.csv `columns`, e.g. 'effect_allele=ALT;other_allele=REF'

    Returns:
        field -> source column, for every field found, plus 'marker' (a column
        holding chr:pos, used when there is no chromosome/position column)
    """
    by_lower = {name.lower(): name for name in df.columns}
    found = {}
    for field, aliases in ALIASES.items():
        name = next((by_lower[a] for a in aliases if a in by_lower), None)
        if field == 'effect_allele_frequency' and name is None:
            # eaf_hapmap_* (Wheeler 2017) is a reference panel frequency, not the
            # study's, so it is not used.
            name = next(
                (
                    c
                    for c in df.columns
                    if c.lower().startswith('eaf')
                    and not c.lower().startswith('eaf_hapmap')
                ),
                None,
            )
        if name is not None:
            found[field] = name
    for name in (by_lower[c] for c in RSID_CANDIDATES if c in by_lower):
        if looks_like(df[name], r'^rs\d+'):
            found['rsid'] = name
            break
    if 'chromosome' not in found or 'base_pair_location' not in found:
        for name in (by_lower[c] for c in MARKER_CANDIDATES if c in by_lower):
            if looks_like(df[name], MARKER_PATTERN):
                found['marker'] = name
                break
    for item in filter(None, override.split(';')):
        field, name = item.split('=', 1)
        if name not in df.columns:
            raise ValueError(f'override {item}: no column {name}')
        found[field] = name
    return found


def standardise(df, columns: dict[str, str], n_study: int, counts: dict):
    """
    Build the output fields from the source columns, before liftover.

    Rows without a usable chromosome, position or pair of A/C/G/T alleles are
    dropped and counted in `counts`.

    Returns:
        (table, notes) where notes says where beta, SE and N came from
    """
    import polars as pl

    def text(field):
        expr = pl.col(columns[field]).str.strip_chars()
        return (
            pl.when(expr.str.to_uppercase().is_in(NULL_TOKENS))
            .then(None)
            .otherwise(expr)
        )

    def number(field):
        # Non-finite values (inf, Infinity, 1e999) become NA like any unparseable one.
        value = text(field).cast(pl.Float64, strict=False)
        return pl.when(value.is_finite()).then(value).otherwise(None)

    def has(field):
        return field in columns

    if not (has('effect_allele') and has('other_allele')):
        raise ValueError(f'no alleles: {columns}')
    if not (has('p_value') or has('neg_log_10_p_value')):
        raise ValueError(f'no P value: {columns}')
    notes = {}
    out = {}

    if has('chromosome') and has('base_pair_location'):
        chromosome, position = text('chromosome'), number('base_pair_location')
    else:
        if not has('marker'):
            raise ValueError(f'no chromosome/position or chr:pos column: {columns}')
        marker = text('marker')
        chromosome = marker.str.extract(MARKER_PATTERN, 1)
        position = marker.str.extract(MARKER_PATTERN, 2).cast(pl.Float64)
        notes['position_source'] = f'parsed from {columns["marker"]}'
    out['chromosome'] = (
        chromosome.str.to_uppercase()
        .str.replace(r'^CHR', '')
        .replace_strict(CHROMOSOME_CODES, default=None, return_dtype=pl.Int64)
    )
    out['base_pair_location'] = (
        pl.when((position >= 1) & (position == position.floor()))
        .then(position)
        .otherwise(None)
        .cast(pl.Int64)
    )
    out['effect_allele'] = text('effect_allele').str.to_uppercase()
    out['other_allele'] = text('other_allele').str.to_uppercase()
    if has('indels_from'):
        # Indels coded D/I (deletion = shorter allele, insertion = longer; checked for
        # Chen 2020 against 1000 Genomes frequencies), with the sequences in an ID
        # such as 10:100009580_C_CCT: take each allele's sequence from there.
        ids = text('indels_from')
        first = ids.str.extract(r'_([ACGTacgt]+)_([ACGTacgt]+)$', 1).str.to_uppercase()
        second = ids.str.extract(r'_([ACGTacgt]+)_([ACGTacgt]+)$', 2).str.to_uppercase()
        longer_first = first.str.len_chars() > second.str.len_chars()
        indel = first.str.len_chars() != second.str.len_chars()
        shorter = pl.when(longer_first).then(second).otherwise(first)
        longer = pl.when(longer_first).then(first).otherwise(second)
        for field in ('effect_allele', 'other_allele'):
            coded = out[field]
            out[field] = (
                pl.when((coded == 'D') & indel)
                .then(shorter)
                .when((coded == 'I') & indel)
                .then(longer)
                .otherwise(coded)
            )
        notes['indel_alleles'] = f'D/I sequences from {columns["indels_from"]}'

    if has('beta'):
        out['beta'] = number('beta')
        notes['effect_source'] = 'beta'
    elif has('odds_ratio'):
        odds = number('odds_ratio')
        out['beta'] = pl.when(odds > 0).then(odds.log()).otherwise(None)
        notes['effect_source'] = 'ln(odds_ratio)'
    elif has('z') and has('standard_error'):
        out['beta'] = number('z') * number('standard_error')
        notes['effect_source'] = 'z x standard_error'
    else:
        out['beta'] = pl.lit(None, pl.Float64)
        notes['effect_source'] = 'z only' if has('z') else 'none (P value only)'
    if has('standard_error'):
        out['standard_error'] = number('standard_error')
    elif has('odds_ratio') and has('ci_lower') and has('ci_upper'):
        lower, upper = number('ci_lower'), number('ci_upper')
        out['standard_error'] = (
            pl.when((lower > 0) & (upper > lower))
            .then((upper.log() - lower.log()) / (2 * Z_95))
            .otherwise(None)
        )
        notes['effect_source'] += ', SE from 95% CI'
    else:
        out['standard_error'] = pl.lit(None, pl.Float64)

    frequency = (
        number('effect_allele_frequency') if has('effect_allele_frequency') else None
    )
    out['effect_allele_frequency'] = (
        pl.when((frequency >= 0) & (frequency <= 1)).then(frequency).otherwise(None)
        if frequency is not None
        else pl.lit(None, pl.Float64)
    )

    if has('p_value'):
        p_text = text('p_value').str.replace(r'^\.', '0.')
        mantissa = p_text.str.extract(r'^([0-9.]+)', 1).cast(pl.Float64, strict=False)
        exponent = p_text.str.extract(r'[eE]([-+]?\d+)$', 1).cast(pl.Int64).fill_null(0)
        minus_log = -mantissa.log10() - exponent
        numeric = p_text.str.contains(P_TEXT_PATTERN).fill_null(False)
        zero = numeric & (mantissa == 0).fill_null(False)
        valid = (numeric & (mantissa > 0) & (minus_log > -1e-6)).fill_null(False)
        out['p_value'] = pl.when(valid).then(p_text).otherwise(None)
        out['neg_log_10_p_value'] = (
            pl.when(valid).then(minus_log).otherwise(None).round(6)
        )
        checked = df.select(p=p_text, ok=valid, zero=zero)
        counts['p_value_zero'] = int(checked['zero'].sum())
        unparseable = checked.filter(
            pl.col('p').is_not_null() & ~pl.col('ok') & ~pl.col('zero')
        )['p']
        counts['p_value_unparseable'] = unparseable.len()
        notes['p_value_unparseable_examples'] = unparseable.unique().head(5).to_list()
    else:
        raw_text = text('neg_log_10_p_value')
        checked = df.select(t=raw_text, v=raw_text.cast(pl.Float64, strict=False))
        lowest = checked['v'].min()
        if lowest is not None and lowest < 0:
            raise ValueError(
                f'{columns["neg_log_10_p_value"]} has negative values: is it '
                'log10(P) rather than -log10(P)?'
            )
        # +inf is an underflowed P, like P = 0 in a P value column.
        counts['p_value_zero'] = checked.filter(pl.col('v').is_infinite()).height
        unparseable = checked.filter(pl.col('t').is_not_null() & pl.col('v').is_null())
        counts['p_value_unparseable'] = unparseable.height
        notes['p_value_unparseable_examples'] = (
            unparseable['t'].unique().head(5).to_list()
        )
        minus_log = number('neg_log_10_p_value')
        exponent = minus_log.ceil()
        mantissa = (10 ** (exponent - minus_log)).round(4)
        out['p_value'] = (
            mantissa.cast(pl.Utf8)
            + pl.lit('e-')
            + exponent.cast(pl.Int64).cast(pl.Utf8)
        )
        out['neg_log_10_p_value'] = minus_log
        notes['p_value_source'] = 'from neg_log_10_p_value'

    out['rsid'] = text('rsid') if has('rsid') else pl.lit(None, pl.Utf8)
    if has('n'):
        out['n'] = number('n').round(0).cast(pl.Int64)
        notes['n_source'] = 'per variant'
    else:
        out['n'] = pl.lit(n_study or None, pl.Int64)
        notes['n_source'] = 'study total (files.csv n_study)' if n_study else 'none'
    out['z'] = number('z') if has('z') else pl.lit(None, pl.Float64)

    # Safety net: no inf or NaN reaches the output, whatever path made it (write_csv
    # would print them as text, not as NA).
    table = df.select(**out).with_columns(
        pl.when(pl.col(column).is_finite()).then(pl.col(column)).otherwise(None)
        for column in FLOAT_COLUMNS
    )
    counts['rows_in'] = table.height
    for reason, condition in [
        ('non_canonical_chromosome', pl.col('chromosome').is_null()),
        ('missing_or_bad_position', pl.col('base_pair_location').is_null()),
        (
            'non_acgt_allele',
            ~pl.col('effect_allele').str.contains(r'^[ACGT]+$').fill_null(False)
            | ~pl.col('other_allele').str.contains(r'^[ACGT]+$').fill_null(False),
        ),
        ('identical_alleles', pl.col('effect_allele') == pl.col('other_allele')),
    ]:
        dropped = table.filter(condition)
        counts[f'dropped_{reason}'] = dropped.height
        table = table.filter(~condition)
    return table, notes


def complement(expr):
    """Complement each base of an allele column (A<->T, C<->G), not reversed."""
    return expr.str.replace_many(['A', 'C', 'G', 'T'], ['T', 'G', 'C', 'A'])


def lift_position(lifter, contig: str, position: int) -> list[tuple[str, int, str]]:
    """
    Lift one 1-based position with pyliftover, whose chain coordinates are 0-based.

    Returns:
        (contig, 1-based position, strand) per hit
    """
    hits = lifter.convert_coordinate(contig, position - 1) or []
    return [(hit[0], hit[1] + 1, hit[2]) for hit in hits]


def check_liftover_offsets(lifter, build: str) -> None:
    """Fail unless every LIFTOVER_KNOWN variant lands on its known GRCh38 position."""
    for chromosome, source, grch38 in LIFTOVER_KNOWN[build]:
        hits = lift_position(lifter, CHAIN_CONTIGS[chromosome], source)
        if not (hits and hits[0][1] == grch38):
            raise RuntimeError(
                f'{build} liftover offset check failed: chr{chromosome}:{source} -> '
                f'{hits}, expected {grch38}'
            )


def lift_to_grch38(table, chain: Path, build: str, counts: dict):
    """
    Lift GRCh37 or NCBI36 positions to GRCh38.

    Kept: positions with exactly one hit on a canonical chromosome. Minus-strand
    hits have their alleles complemented (SNVs only; multi-base alleles on the
    minus strand are dropped, as their left-aligned form would change).
    """
    import polars as pl
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    check_liftover_offsets(lifter, build)
    canonical = {name: code for code, name in CHAIN_CONTIGS.items()}
    new_chromosome, new_position, minus = [], [], []
    reasons = {
        'unmapped': 0,
        'multiple_hits': 0,
        'non_canonical_target': 0,
        # hg19/hg18 chrM is not the GRCh38 MT sequence and the chains omit it.
        'mt_not_lifted': 0,
    }
    for chromosome, position in zip(
        table['chromosome'].to_list(),
        table['base_pair_location'].to_list(),
        strict=True,
    ):
        contig = CHAIN_CONTIGS.get(chromosome)
        hits = lift_position(lifter, contig, position) if contig else None
        reason = None
        if chromosome == 25:
            reason = 'mt_not_lifted'
        elif not hits:
            reason = 'unmapped'
        elif len(hits) > 1:
            reason = 'multiple_hits'
        elif hits[0][0] not in canonical:
            reason = 'non_canonical_target'
        if reason:
            reasons[reason] += 1
            new_chromosome.append(None)
            new_position.append(None)
            minus.append(None)
            continue
        new_chromosome.append(canonical[hits[0][0]])
        new_position.append(hits[0][1])
        minus.append(hits[0][2] == '-')
    counts.update({f'liftover_{k}': v for k, v in reasons.items()})
    lifted = table.with_columns(
        _chromosome=pl.Series(new_chromosome, dtype=pl.Int64),
        _position=pl.Series(new_position, dtype=pl.Int64),
        _minus=pl.Series(minus, dtype=pl.Boolean),
    ).filter(pl.col('_position').is_not_null())
    counts['liftover_chromosome_changed'] = lifted.filter(
        pl.col('_chromosome') != pl.col('chromosome')
    ).height
    snv = (pl.col('effect_allele').str.len_chars() == 1) & (
        pl.col('other_allele').str.len_chars() == 1
    )
    minus_multibase = pl.col('_minus') & ~snv
    counts['liftover_minus_strand_multibase_dropped'] = lifted.filter(
        minus_multibase
    ).height
    counts['liftover_minus_strand_complemented'] = lifted.filter(
        pl.col('_minus') & snv
    ).height
    lifted = lifted.filter(~minus_multibase).with_columns(
        chromosome=pl.col('_chromosome'),
        base_pair_location=pl.col('_position'),
        effect_allele=pl.when('_minus')
        .then(complement(pl.col('effect_allele')))
        .otherwise('effect_allele'),
        other_allele=pl.when('_minus')
        .then(complement(pl.col('other_allele')))
        .otherwise('other_allele'),
    )
    return lifted.drop('_chromosome', '_position', '_minus')


def check_build(counts: dict) -> None:
    """
    Fail if the positions are not on GRCh38, or if no row survived the reference
    check. Any other SNP matches on one strand or the other wherever it is placed,
    so the test uses A/T and C/G SNPs: at the right positions they match, at wrong
    ones about half do not. A file with too few of them (Timsina 2026, whose authors
    removed them) is judged on its indels instead, whose alleles must match the
    reference as written: about 75% do not at wrong positions. Either way, more than
    MAX_BUILD_MISMATCH mismatching fails. Takes the reference_* counts summed over
    every chromosome, and adds the mismatch rates, reference_build_check (which test
    ran) and reference_build_checked.
    """
    tests = {
        'palindromic SNPs': ('palindromic', 'A/T and C/G SNPs'),
        'indels': ('multibase', 'indels'),
    }
    for key, _ in tests.values():
        total = counts.get(f'reference_{key}', 0)
        rate = counts.get(f'reference_{key}_mismatch', 0) / total if total else 0.0
        counts[f'reference_{key}_mismatch_rate'] = round(rate, 5)
    for name, (key, label) in tests.items():
        rate = counts[f'reference_{key}_mismatch_rate']
        if counts.get(f'reference_{key}', 0) >= MIN_SITES_FOR_BUILD_CHECK:
            counts['reference_build_check'] = name
            counts['reference_build_checked'] = True
            if rate > MAX_BUILD_MISMATCH:
                raise ValueError(
                    f'{rate:.1%} of {label} do not match GRCh38: positions are not '
                    'on the expected build. Correct source_build in files.csv, then '
                    'rerun format.py --only <file_id> --force (the download is reused)'
                )
            break
    else:
        counts['reference_build_check'] = 'none: too few A/T, C/G SNPs and indels'
        counts['reference_build_checked'] = False
    if not counts.get('reference_ok', 0) + counts.get('reference_strand_flipped', 0):
        raise ValueError(
            f'no rows survived the GRCh38 reference check, so nothing is written; '
            f'counts: {dict(counts)}'
        )


def check_reference(table, fasta: Path, counts: dict, verdict: bool = True):
    """
    Check that one of the two alleles matches GRCh38 at the position.

    SNVs that match only after complementing are strand-flipped (both alleles
    complemented; beta still refers to the same, now complemented, effect allele).
    Multi-base alleles must match the reference as written. Rows matching
    neither way are dropped.

    With verdict, also runs check_build on this table's counts. format_one checks
    one chromosome at a time, so it passes verdict=False and calls check_build once
    on the totals.
    """
    import numpy as np
    import polars as pl
    import pysam

    reference = pysam.FastaFile(str(fasta))
    contigs = ENSEMBL_CONTIGS if '1' in reference.references else UCSC_CONTIGS
    kept = []
    tally = {
        'ok': 0,
        'strand_flipped': 0,
        'mismatch': 0,
        'palindromic': 0,
        'palindromic_mismatch': 0,
        'multibase': 0,
        'multibase_mismatch': 0,
    }
    for (chromosome,), part in table.group_by(['chromosome'], maintain_order=True):
        sequence = reference.fetch(contigs[chromosome]).upper()
        bases = np.frombuffer(sequence.encode(), dtype='S1')
        positions = part['base_pair_location'].to_numpy()
        inside = positions <= len(bases)
        ref_base = np.where(inside, bases[np.minimum(positions, len(bases)) - 1], b'N')
        part = part.with_columns(_ref=pl.Series(ref_base.astype(str)))
        effect, other = pl.col('effect_allele'), pl.col('other_allele')
        snv = (effect.str.len_chars() == 1) & (other.str.len_chars() == 1)
        flipped_effect, flipped_other = complement(effect), complement(other)
        palindromic = snv & (effect == flipped_other)
        snvs = part.filter(snv).with_columns(
            _status=pl.when((effect == pl.col('_ref')) | (other == pl.col('_ref')))
            .then(pl.lit('ok'))
            .when(
                (flipped_effect == pl.col('_ref')) | (flipped_other == pl.col('_ref'))
            )
            .then(pl.lit('strand_flipped'))
            .otherwise(pl.lit('mismatch')),
            _palindromic=palindromic,
        )
        snvs = snvs.with_columns(
            effect_allele=pl.when(pl.col('_status') == 'strand_flipped')
            .then(flipped_effect)
            .otherwise(effect),
            other_allele=pl.when(pl.col('_status') == 'strand_flipped')
            .then(flipped_other)
            .otherwise(other),
        )
        multibase = part.filter(~snv)
        status = [
            'ok'
            if any(
                sequence[pos - 1 : pos - 1 + len(allele)] == allele
                for allele in alleles
            )
            else 'mismatch'
            for pos, *alleles in zip(
                multibase['base_pair_location'].to_list(),
                multibase['effect_allele'].to_list(),
                multibase['other_allele'].to_list(),
                strict=True,
            )
        ]
        multibase = multibase.with_columns(
            _status=pl.Series(status, dtype=pl.Utf8), _palindromic=pl.lit(False)
        )
        tally['multibase'] += multibase.height
        tally['multibase_mismatch'] += status.count('mismatch')
        part = pl.concat([snvs, multibase])
        for key, value in part['_status'].value_counts().iter_rows():
            tally[key] += value
        tally['palindromic'] += part['_palindromic'].sum()
        tally['palindromic_mismatch'] += part.filter(
            pl.col('_palindromic') & (pl.col('_status') == 'mismatch')
        ).height
        kept.append(part.filter(pl.col('_status') != 'mismatch'))
    counts.update({f'reference_{k}': int(v) for k, v in tally.items()})
    if verdict:
        check_build(counts)
    if not kept:
        return table.head(0)
    return pl.concat(kept).drop('_ref', '_status', '_palindromic')


SORT_KEY = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele']


def chromosome_key(value: str, from_marker: bool) -> int:
    """GWAS-SSF code for a raw chromosome (or chr:pos marker) value; 0 if unknown."""
    if from_marker:
        match = re.match(MARKER_PATTERN, value)
        value = match.group(1) if match else ''
    value = value.strip().upper()
    if value.startswith('CHR'):
        value = value[3:]
    return CHROMOSOME_CODES.get(value, 0)


def split_by_chromosome(
    path: Path, workdir: Path, columns: dict[str, str]
) -> dict[int, Path]:
    """
    Stream a source file into one tab-separated file per chromosome, so each can
    be formatted with only that chromosome in memory. Rows whose chromosome is not
    recognised go to key 0, where standardise drops and counts them.

    Returns:
        chromosome code -> file, in chromosome order
    """
    from_marker = not ('chromosome' in columns and 'base_pair_location' in columns)
    field = columns['marker'] if from_marker else columns['chromosome']
    handles: dict = {}
    try:
        with open_text(path) as source:
            lines = (line for line in source if not line.startswith('##'))
            header = next(lines)
            separator = detect_separator(path)
            names = [name.strip() for name in header.rstrip('\r\n').split(separator)]
            index = names.index(field)
            for line in lines:
                fields = line.rstrip('\r\n').split(separator)
                if not fields or fields == ['']:
                    continue
                value = fields[index] if index < len(fields) else ''
                key = chromosome_key(value, from_marker)
                if key not in handles:
                    handles[key] = (workdir / f'chromosome_{key}.tsv').open('w')
                    handles[key].write('\t'.join(names) + '\n')
                handles[key].write('\t'.join(fields) + '\n')
    finally:
        for handle in handles.values():
            handle.close()
    return {key: workdir / f'chromosome_{key}.tsv' for key in sorted(handles)}


def drop_duplicates(table, counts: dict):
    """
    One row per variant (chromosome, position, effect and other allele). Exact copies
    keep one row; copies that disagree (Verma 2024 has the same variant twice with
    different beta and P) are all dropped, since nothing says which is right.
    Counted as dropped_duplicate_exact and dropped_duplicate_conflicting.
    """
    import polars as pl

    unique = table.unique(maintain_order=True)
    counts['dropped_duplicate_exact'] += table.height - unique.height
    conflicting = pl.struct(SORT_KEY).is_duplicated()
    counts['dropped_duplicate_conflicting'] += unique.filter(conflicting).height
    return unique.filter(~conflicting)


def format_one(
    row: dict,
    original: Path,
    out_root: Path,
    stats_path: Path,
    fasta: Path,
    chains: dict[str, Path],
    download_record: dict,
) -> None:
    """
    Format one file: writes <out_root>.tsv.gz, its .tbi, and the stats JSON.

    Args:
        row: one files.csv row
        original: the downloaded file
        out_root: output path without the .tsv.gz suffix
        stats_path: where to write the per-file stats JSON
        fasta: GRCh38 fasta, Ensembl or UCSC contig names; .fai built if missing
        chains: source build -> chain file to GRCh38
        download_record: the JSON record download.py wrote next to the original
    """
    import pysam

    import polars as pl

    build = row['source_build']
    if build not in LIFTOVER_KNOWN and build != 'GRCh38':
        raise ValueError(f'unsupported source_build {build}')
    if not Path(f'{fasta}.fai').exists():
        pysam.faidx(str(fasta))
    counts: Counter = Counter()
    notes: dict = {}
    bad_p: set[str] = set()
    with tempfile.TemporaryDirectory(dir=out_root.parent) as tmp:
        workdir = Path(tmp)
        sample = read_table(original, workdir, n_rows=10_000)
        columns = resolve_columns(sample, row['columns'])
        used = sorted(set(columns.values()))
        # One source chromosome at a time; results are kept on disk per GRCh38
        # chromosome, since liftover can move a few rows to another chromosome.
        parts: dict[int, list[Path]] = defaultdict(list)
        for key, path in split_by_chromosome(original, workdir, columns).items():
            source = read_table(path, workdir, columns=used)
            path.unlink()
            step: dict = {}
            table, notes = standardise(source, columns, int(row['n_study'] or 0), step)
            del source
            bad_p.update(notes.pop('p_value_unparseable_examples', []))
            if build in LIFTOVER_KNOWN:
                table = lift_to_grch38(table, chains[build], build, step)
            table = check_reference(table, fasta, step, verdict=False)
            counts.update(step)
            for (chromosome,), part in table.group_by(['chromosome']):
                part_path = workdir / f'part_{key}_{chromosome}.parquet'
                part.write_parquet(part_path)
                parts[chromosome].append(part_path)
        check_build(counts)
        if counts['p_value_unparseable'] > MAX_UNPARSEABLE_P * counts['rows_in']:
            raise ValueError(
                f'{counts["p_value_unparseable"]:,} of {counts["rows_in"]:,} P values '
                f'are not numbers in (0, 1], e.g. {sorted(bad_p)[:5]}: wrong column, '
                'or a format to handle?'
            )
        if bad_p:
            notes['p_value_unparseable_examples'] = sorted(bad_p)[:5]
        plain = workdir / 'formatted.tsv'
        filled = Counter()
        with plain.open('wb') as out:
            out.write(('\t'.join(OUTPUT_COLUMNS) + '\n').encode())
            for chromosome in sorted(parts):
                table = pl.concat([pl.read_parquet(p) for p in parts[chromosome]])
                table = drop_duplicates(table.sort(SORT_KEY), counts)
                counts['rows_out'] += table.height
                filled.update(
                    {c: table.height - table[c].null_count() for c in OUTPUT_COLUMNS}
                )
                table.select(OUTPUT_COLUMNS).write_csv(
                    out, separator='\t', null_value='NA', include_header=False
                )
        empty = [c for c in OUTPUT_COLUMNS if not filled[c]]
        compressed = pysam.tabix_index(
            str(plain), seq_col=0, start_col=1, end_col=1, line_skip=1, force=True
        )
        shutil.move(compressed, f'{out_root}.tsv.gz')
        shutil.move(f'{compressed}.tbi', f'{out_root}.tsv.gz.tbi')
    stats = {
        'file_id': row['file_id'],
        'formatted_file': f'{Path(out_root).name}.tsv.gz',
        'source_build': row['source_build'],
        'original_md5': download_record['md5'],
        'md5_checked': download_record['md5_checked'],
        'downloaded_at': download_record['downloaded_at'],
        'columns_used': columns,
        'columns_all_missing': empty,
        **notes,
        **dict(counts),
        'formatted_at': datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }
    stats_path.write_text(json.dumps(stats, indent=1) + '\n')
    print(
        f'{row["file_id"]}: {counts["rows_in"]:,} rows in, {counts["rows_out"]:,} out'
    )


COLUMNS_USED_FIELDS = [
    'chromosome',
    'base_pair_location',
    'marker',
    'effect_allele',
    'other_allele',
    'beta',
    'odds_ratio',
    'standard_error',
    'effect_allele_frequency',
    'p_value',
    'neg_log_10_p_value',
    'z',
    'n',
    'rsid',
]


def check_one_source(row: dict) -> dict:
    """
    Resolve and standardise the first lines of one source file, read from its URL.

    Returns:
        one columns_used.tsv line: the source column behind each output field,
        where beta and N come from, and how many sample rows survived
    """
    line = {'file_id': row['file_id'], 'override': row['columns']}
    try:
        with tempfile.TemporaryDirectory() as tmp:
            head = Path(tmp) / 'head.tsv'
            head.write_text(fetch_head(row['source_url']))
            source = read_table(head, Path(tmp))
            columns = resolve_columns(source, row['columns'])
            counts: dict = {}
            table, notes = standardise(
                source, columns, int(row['n_study'] or 0), counts
            )
        line |= {field: columns.get(field, '') for field in COLUMNS_USED_FIELDS}
        line |= {
            'status': 'ok',
            'effect_source': notes['effect_source'],
            'n_source': notes['n_source'],
            'sample_rows_in': counts['rows_in'],
            'sample_rows_kept': table.height,
        }
    except (OSError, ValueError) as error:
        line |= {'status': f'error: {error}'}
    return line


def check_columns(rows: list[dict], out: Path) -> bool:
    """
    Run check_one_source on every row (8 at a time) and write the result to `out`.

    Returns:
        True if every file resolved
    """
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(8) as pool:
        lines = list(pool.map(check_one_source, rows))
    fields = ['file_id', 'status', 'override', *COLUMNS_USED_FIELDS]
    fields += ['effect_source', 'n_source', 'sample_rows_in', 'sample_rows_kept']
    with out.open('w', newline='') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter='\t', lineterminator='\n'
        )
        writer.writeheader()
        writer.writerows(lines)
    failed = [line for line in lines if line['status'] != 'ok']
    for line in failed:
        print(f'{line["file_id"]}: {line["status"]}')
    print(f'{len(lines) - len(failed)} of {len(lines)} files resolved -> {out}')
    return not failed


def local_chains(args) -> dict[str, str | None]:
    return {'GRCh37': args.chain_grch37, 'NCBI36': args.chain_ncbi36}


def check_local_args(rows: list[dict], args, parser) -> None:
    """Stop early, with the flag to add, if a --local run lacks a local input."""
    remote = [
        f'--{name} {getattr(args, name)}'
        for name in ('originals', 'out', 'stats', 'web')
        if getattr(args, name).startswith('gs://')
    ]
    if remote:
        parser.error(f'--local needs local folders, got: {", ".join(remote)}')
    if not args.fasta:
        parser.error('--local needs --fasta (a GRCh38 fasta)')
    for build, chain in local_chains(args).items():
        if not chain and any(row['source_build'] == build for row in rows):
            flag = f'--chain-{build.lower()}'
            parser.error(f'--local: some selected files are {build}, pass {flag}')


def run_local(rows: list[dict], args) -> None:
    """Format rows one after another, local folders in and out."""
    originals, out, stats = Path(args.originals), Path(args.out), Path(args.stats)
    out.mkdir(parents=True, exist_ok=True)
    stats.mkdir(parents=True, exist_ok=True)
    for row in rows:
        out_root = out / formatted_name(row)
        if Path(f'{out_root}.tsv.gz').exists() and not args.force:
            print(f'{row["file_id"]}: exists, skipped')
            continue
        record = originals / f'{original_name(row)}.json'
        format_one(
            row,
            originals / original_name(row),
            out_root,
            stats / f'{row["file_id"]}.json',
            Path(args.fasta),
            {build: Path(chain) for build, chain in local_chains(args).items() if chain},
            json.loads(record.read_text()),
        )


def run_batch(rows: list[dict], args) -> None:
    """
    One Hail Batch job per downloaded file not yet formatted. Each job gets this
    script and download.py (for the shared helpers) inline, the original, the fasta
    and the chain file; Batch copies its outputs only if the job succeeds.
    """
    from cpg_utils import to_path
    from cpg_utils.config import config_retrieve, reference_path
    from cpg_utils.hail_batch import get_batch

    fasta = args.fasta or str(reference_path('ensembl_113/unmasked_reference'))
    chains = {
        'GRCh37': args.chain_grch37 or str(reference_path('liftover_37_to_38')),
        'NCBI36': args.chain_ncbi36 or str(reference_path('liftover_36_to_38')),
    }
    batch = get_batch(name='Format biomarker GWAS summary statistics')
    scripts = {
        name: (Path(__file__).with_name(name)).read_text()
        for name in ('download.py', 'format.py')
    }
    if to_path(f'{fasta}.fai').exists():
        # pysam looks for <fasta>.fai. hailtop 0.2.139 keeps each input's own file
        # name; older versions name group members {root}.<key>. The keys 'fa' and
        # 'fa.fai' put the index next to the fasta under either scheme.
        fasta_group = batch.read_input_group(**{'fa': fasta, 'fa.fai': f'{fasta}.fai'})
        fasta_input = fasta_group['fa']
    else:
        fasta_input = batch.read_input(fasta)
    chain_inputs = {build: batch.read_input(path) for build, path in chains.items()}
    submitted, not_submitted = 0, []
    for row in rows:
        original = f'{args.originals}/{original_name(row)}'
        out_root = f'{args.out}/{formatted_name(row)}'
        if to_path(f'{out_root}.tsv.gz').exists() and not args.force:
            print(f'{row["file_id"]}: exists, skipped')
            continue
        if not to_path(original).exists():
            not_submitted.append(f'{row["file_id"]} (not downloaded)')
            continue
        record = to_path(f'{original}.json')
        if not record.exists():
            not_submitted.append(f'{row["file_id"]} (no download record)')
            continue
        size_gib = int(row['size_bytes'] or 0) / 2**30
        job = batch.new_bash_job(f'format {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(8 if size_gib > 1 else 4)
        job.memory('highmem')
        # decompressed copy (gzip ratio up to ~6) + output + fasta
        job.storage(f'{int(size_gib * 12) + 20}Gi')
        job.declare_resource_group(
            out={'tsv.gz': '{root}.tsv.gz', 'tsv.gz.tbi': '{root}.tsv.gz.tbi'}
        )
        write_scripts = ''.join(
            f"cat > {name} <<'GWAS_SUMSTATS_SCRIPT'\n{text}\nGWAS_SUMSTATS_SCRIPT\n"
            for name, text in scripts.items()
        )
        job.command(
            f'python3 -m pip install --quiet {PIP_PACKAGES}\n'
            f'{write_scripts}'
            'mkdir -p work\n'
            f'python3 format.py --one {shlex.quote(json.dumps(row))} '
            f'--original {batch.read_input(original)} --out-root work/out '
            f'--stats-file {job.stats} --fasta {fasta_input} '
            f'--download-record {shlex.quote(record.read_text())} '
            f'--chain-grch37 {chain_inputs["GRCh37"]} '
            f'--chain-ncbi36 {chain_inputs["NCBI36"]}\n'
            f'mv work/out.tsv.gz {job.out["tsv.gz"]}\n'
            f'mv work/out.tsv.gz.tbi {job.out["tsv.gz.tbi"]}'
        )
        batch.write_output(job.out, out_root)
        batch.write_output(job.stats, f'{args.stats}/{row["file_id"]}.json')
        submitted += 1
    print(f'{submitted} format jobs submitted')
    if not_submitted:
        print(
            f'{len(not_submitted)} files NOT submitted, so they will be missing from '
            f'{args.out} until downloaded:\n  ' + '\n  '.join(not_submitted)
        )
    if submitted:
        batch.run(wait=False)


MANIFEST_FIELDS = [
    'formatted_file',
    'file_id',
    'first_author',
    'year',
    'journal',
    'biomarker',
    'trait',
    'ancestry',
    'gcst',
    'n_study',
    'source_kind',
    'source_build',
    'build_evidence',
    'source_url',
    'original_md5',
    'md5_checked',
    'downloaded_at',
    'formatted_at',
    'effect_source',
    'n_source',
    'rows_in',
    'rows_out',
    'columns_all_missing',
    'notes',
    'columns_used',
]


def write_manifest(rows: list[dict], args) -> None:
    """
    manifest.csv next to the formatted files, plus an HTML copy for the web bucket.
    One line per formatted file, joining files.csv and the formatting stats (every
    count in the stats JSON becomes a column). Fails if a formatted file has no stats,
    rather than writing a manifest that leaves it out. Fails too if any files.csv row
    is not formatted, unless args.allow_missing; then the missing file_ids are listed
    in the log and in manifest.html.
    """
    out_dir, stats_dir = as_path(args.out), as_path(args.stats)
    lines, missing_stats, not_formatted = [], [], []
    for row in rows:
        stats_file = stats_dir / f'{row["file_id"]}.json'
        if not stats_file.exists():
            if (out_dir / f'{formatted_name(row)}.tsv.gz').exists():
                missing_stats.append(row['file_id'])
            else:
                not_formatted.append(row['file_id'])
            continue
        stats = json.loads(stats_file.read_text())
        line = {k: row.get(k, '') for k in MANIFEST_FIELDS}
        for key, value in stats.items():
            if key == 'columns_used':
                line[key] = ';'.join(f'{k}={v}' for k, v in value.items())
            elif key not in ('file_id', 'source_build'):
                line[key] = ';'.join(value) if isinstance(value, list) else value
        lines.append(line)
    if missing_stats:
        raise ValueError(
            f'{len(missing_stats)} formatted files have no stats in {stats_dir}: '
            + ', '.join(missing_stats)
        )
    if not lines:
        raise ValueError(f'no stats found in {stats_dir}')
    if not_formatted:
        listing = f'{len(not_formatted)} files in files.csv are not formatted: ' + (
            ', '.join(not_formatted)
        )
        if not args.allow_missing:
            raise ValueError(f'{listing}. Format them, or pass --allow-missing.')
        print(f'{listing}. Left out of the manifest (--allow-missing).')
    fields = MANIFEST_FIELDS + sorted(
        {k for line in lines for k in line} - set(MANIFEST_FIELDS)
    )
    out = as_path(args.out) / 'manifest.csv'
    with out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(lines)
    cells = ''.join(f'<th>{html.escape(f)}</th>' for f in fields)
    body = ''.join(
        '<tr>'
        + ''.join(f'<td>{html.escape(str(line.get(f, "")))}</td>' for f in fields)
        + '</tr>'
        for line in lines
    )
    page = (
        '<!doctype html><meta charset="utf-8"><title>GWAS summary statistics</title>'
        '<style>body{font:13px sans-serif}table{border-collapse:collapse}'
        'td,th{border:1px solid #ccc;padding:2px 4px;white-space:nowrap}'
        'th{position:sticky;top:0;background:#eef}</style>'
        f'<h1>GWAS summary statistics</h1><p>{len(lines)} files in '
        f'{html.escape(args.out)}, written {datetime.now(timezone.utc):%Y-%m-%d}.</p>'
        + (
            f'<p><b>{len(not_formatted)} files in files.csv are not formatted:</b> '
            f'{html.escape(", ".join(not_formatted))}</p>'
            if not_formatted
            else ''
        )
        + f'<table><tr>{cells}</tr>{body}</table>'
    )
    web = as_path(args.web)
    if isinstance(web, Path):
        web.mkdir(parents=True, exist_ok=True)
    (web / 'manifest.html').write_text(page)
    print(f'manifest: {len(lines)} files -> {out} and {web}/manifest.html')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--files', type=Path, default=FILES_CSV)
    parser.add_argument('--originals', default=DEFAULT_ORIGINALS)
    parser.add_argument('--out', default=DEFAULT_OUT)
    parser.add_argument('--stats', default=DEFAULT_STATS, help='per-file stats folder')
    parser.add_argument('--web', default=DEFAULT_WEB, help='folder for manifest.html')
    parser.add_argument('--fasta', help='GRCh38 fasta; default: references config')
    parser.add_argument('--chain-grch37', help='default: references liftover_37_to_38')
    parser.add_argument('--chain-ncbi36', help='default: references liftover_36_to_38')
    parser.add_argument('--only', nargs='+', help='file_ids to format')
    parser.add_argument(
        '--force',
        action='store_true',
        help='with --only: redo these files even if formatted (after a files.csv fix)',
    )
    parser.add_argument('--local', action='store_true', help='run here, not on Batch')
    parser.add_argument('--manifest', action='store_true', help='write the manifest')
    parser.add_argument(
        '--allow-missing',
        action='store_true',
        help='with --manifest: write it even if some files.csv rows are not formatted',
    )
    parser.add_argument(
        '--check-columns',
        type=Path,
        nargs='?',
        const=COLUMNS_USED_TSV,
        help='resolve every source header from its URL; writes columns_used.tsv',
    )
    parser.add_argument('--one', help='a single files.csv row as JSON (job mode)')
    parser.add_argument('--original', type=Path, help='job mode: input file')
    parser.add_argument('--out-root', type=Path, help='job mode: output path root')
    parser.add_argument('--stats-file', type=Path, help='job mode: stats JSON path')
    parser.add_argument('--download-record', help='job mode: download record JSON')
    args = parser.parse_args()
    if args.force and not args.only:
        parser.error('--force needs --only, so every file is not redone by mistake')
    if args.allow_missing and not args.manifest:
        parser.error('--allow-missing only applies to --manifest')
    if args.check_columns == COLUMNS_USED_TSV and args.only:
        parser.error(
            '--check-columns with --only would replace columns_used.tsv with a '
            'subset; give an output path: --check-columns subset.tsv'
        )
    if args.manifest and args.only:
        parser.error(
            '--manifest indexes every formatted file and replaces the existing '
            'manifest, so it takes no --only'
        )
    for name in ('originals', 'out', 'stats', 'web'):
        setattr(args, name, getattr(args, name).rstrip('/'))

    if args.one:
        format_one(
            json.loads(args.one),
            args.original,
            args.out_root,
            args.stats_file,
            Path(args.fasta),
            {'GRCh37': Path(args.chain_grch37), 'NCBI36': Path(args.chain_ncbi36)},
            json.loads(args.download_record),
        )
        return
    rows = read_files_csv(args.files, args.only)
    if args.check_columns:
        if not check_columns(rows, args.check_columns):
            raise SystemExit(1)
    elif args.manifest:
        write_manifest(rows, args)
    elif args.local:
        check_local_args(rows, args, parser)
        run_local(rows, args)
    else:
        run_batch(rows, args)


if __name__ == '__main__':
    main()
