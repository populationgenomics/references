#!/usr/bin/env python3

"""
Convert the downloaded GWAS summary statistics to one GRCh38 layout.

For each file in files.csv: map its columns onto the GWAS-SSF names, lift GRCh37 or
NCBI36 positions to GRCh38, check both alleles against the GRCh38 reference, sort, and
write a bgzipped TSV with a tabix index. Nothing is filtered on frequency or P value;
rows are dropped only when they cannot be placed on GRCh38 (reasons are counted).

Output columns, in order (missing values are NA):

    chromosome  base_pair_location  effect_allele  other_allele  beta
    standard_error  effect_allele_frequency  p_value  neg_log_10_p_value
    rsid  n  z

chromosome is GWAS-SSF numeric (1-22, X=23, Y=24, MT=25). base_pair_location is
1-based. beta is per effect allele: odds ratios become ln(OR), with the standard error
taken from the 95% CI if no SE is given. A source with no effect size (P value only)
keeps beta as NA. z is filled only where the source has it.
p_value keeps the source's text, so values below the float range (e.g. 8.51e-3553)
survive; neg_log_10_p_value is computed from that text.

Liftover uses pyliftover, whose chain coordinates are 0-based: positions go in as
pos - 1 and come back + 1. check_liftover_offsets() asserts this on two variants with
known positions before any file is lifted.

The build in files.csv is checked against the data: at the right positions both
alleles of an A/T or C/G SNP match GRCh38 (on one strand or the other), at wrong
positions about half do not. A file whose palindromic SNPs mismatch more than
MAX_PALINDROMIC_MISMATCH fails rather than being written.

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

Locally (needs polars, pysam, pyliftover, numpy):

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
import json
import re
import shlex
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from download import FILES_CSV, original_name, read_files_csv

DEFAULT_ORIGINALS = 'gs://cpg-common-main-tmp/gwas_sumstats/original'
DEFAULT_STATS = 'gs://cpg-common-main/references/gwas_sumstats/v1/stats'
DEFAULT_OUT = 'gs://cpg-common-main/references/gwas_sumstats/v1'
DEFAULT_WEB = 'gs://cpg-common-main-web/gwas_sumstats'
PIP_PACKAGES = 'polars==1.34.0 pysam==0.23.3 pyliftover==0.4.1 numpy'

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
    # allele1/a1 = effect allele follows METAL and GCTA. SAIGE and REGENIE use
    # Allele2 as the effect allele: give those files a `columns` override.
    'effect_allele': ['effect_allele', 'allele1', 'a1', 'ea', 'tested_allele'],
    'other_allele': ['other_allele', 'allele2', 'a2', 'nea', 'non_effect_allele'],
    'beta': ['beta', 'effect', 'b'],
    'standard_error': ['standard_error', 'se', 'stderr', 'sebeta'],
    'effect_allele_frequency': ['effect_allele_frequency', 'eaf', 'freq1', 'af1'],
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
MARKER_PATTERN = r'^(?:chr)?([0-9XYMT]+)[:_](\d+)'
NULL_TOKENS = ['', 'NA', '#NA', '.', 'NAN', 'NULL', 'INF', '-INF', 'NONE']

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
MAX_PALINDROMIC_MISMATCH = 0.10
MIN_PALINDROMIC_FOR_CHECK = 100

# (chromosome, source position, GRCh38 position), 1-based. GRCh37: rs7412 and
# rs429358 (APOE). NCBI36: rs1260326 (GCKR) and rs12565286. GRCh38 from Ensembl.
LIFTOVER_KNOWN = {
    'GRCh37': [(19, 45412079, 44908822), (19, 45411941, 44908684)],
    'NCBI36': [(2, 27584444, 27508073), (1, 711153, 785910)],
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


def open_text(path: Path):
    """Open a text file that may or may not be gzip/bgzip compressed."""
    with path.open('rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    if compressed:
        return gzip.open(path, 'rt', newline='')
    return path.open(newline='')


def read_table(path: Path, workdir: Path):
    """
    Read a summary statistics file with every column as text.

    Tab- and comma-separated files are read directly; whitespace-separated ones
    (runs of spaces) are first rewritten as tab-separated. Lines starting with '##'
    are skipped.

    Args:
        path: the original file
        workdir: scratch folder for the rewritten copy
    """
    import polars as pl

    with open_text(path) as handle:
        header = next(line for line in handle if not line.startswith('##'))
    if '\t' in header:
        separator = '\t'
    elif ',' in header:
        separator = ','
    else:
        rewritten = workdir / 'tab_separated.tsv'
        with open_text(path) as source, rewritten.open('w') as out:
            for line in source:
                if not line.startswith('##'):
                    out.write('\t'.join(line.split()) + '\n')
        path, separator = rewritten, '\t'
    return pl.read_csv(
        path,
        separator=separator,
        infer_schema=False,
        quote_char=None,
        comment_prefix='##',
    )


def looks_like(values, pattern: str) -> bool:
    """True if most non-null values in a sample match the pattern."""
    sample = [v for v in values[:1000] if v]
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
        if looks_like(df[name].to_list(), r'^rs\d+'):
            found['rsid'] = name
            break
    if 'chromosome' not in found or 'base_pair_location' not in found:
        for name in (by_lower[c] for c in MARKER_CANDIDATES if c in by_lower):
            if looks_like(df[name].to_list(), MARKER_PATTERN):
                found['marker'] = name
                break
    for item in filter(None, override.split(';')):
        field, name = item.split('=', 1)
        assert name in df.columns, f'override {item}: no column {name}'
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
        return text(field).cast(pl.Float64, strict=False)

    def has(field):
        return field in columns

    assert has('effect_allele') and has('other_allele'), f'no alleles: {columns}'
    assert has('p_value') or has('neg_log_10_p_value'), f'no P value: {columns}'
    notes = {}
    out = {}

    if has('chromosome') and has('base_pair_location'):
        chromosome, position = text('chromosome'), number('base_pair_location')
    else:
        assert has('marker'), f'no chromosome/position or chr:pos column: {columns}'
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

    if has('beta'):
        out['beta'] = number('beta')
        notes['effect_source'] = 'beta'
    elif has('odds_ratio'):
        odds = number('odds_ratio')
        out['beta'] = pl.when(odds > 0).then(odds.log()).otherwise(None)
        notes['effect_source'] = 'ln(odds_ratio)'
    else:
        out['beta'] = pl.lit(None, pl.Float64)
        notes['effect_source'] = 'z only' if has('z') else 'none (P value only)'
    if has('standard_error'):
        out['standard_error'] = number('standard_error')
    elif has('odds_ratio') and has('ci_lower') and has('ci_upper'):
        out['standard_error'] = (
            number('ci_upper').log() - number('ci_lower').log()
        ) / (2 * Z_95)
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
        out['p_value'] = p_text
        out['neg_log_10_p_value'] = (
            pl.when(mantissa > 0)
            .then(-mantissa.log10() - exponent)
            .otherwise(None)
            .round(6)
        )
    else:
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
        notes['n_source'] = 'study total from GWAS Catalog' if n_study else 'none'
    out['z'] = number('z') if has('z') else pl.lit(None, pl.Float64)

    table = df.select(**out)
    if notes.get('p_value_source') == 'from neg_log_10_p_value':
        lowest = table['neg_log_10_p_value'].min()
        if lowest is not None and lowest < 0:
            raise ValueError(
                f'{columns["neg_log_10_p_value"]} has negative values: is it '
                'log10(P) rather than -log10(P)?'
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


def check_liftover_offsets(lifter, build: str) -> None:
    """Assert 1-based in, 1-based out on variants with known positions."""
    for chromosome, source, grch38 in LIFTOVER_KNOWN[build]:
        hits = lifter.convert_coordinate(CHAIN_CONTIGS[chromosome], source - 1)
        assert hits and hits[0][1] + 1 == grch38, (
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
        hits = lifter.convert_coordinate(contig, position - 1) if contig else None
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
        new_position.append(hits[0][1] + 1)
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


def check_reference(table, fasta: Path, counts: dict):
    """
    Check that one of the two alleles matches GRCh38 at the position.

    SNVs that match only after complementing are strand-flipped (both alleles
    complemented; beta still refers to the same, now complemented, effect allele).
    Multi-base alleles must match the reference as written. Rows matching
    neither way are dropped.

    Fails if more than MAX_PALINDROMIC_MISMATCH of A/T and C/G SNPs mismatch.
    Any other SNP matches on one strand or the other wherever it is placed, so
    only palindromic SNPs show whether positions are on the right build.
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
        part = pl.concat([snvs, multibase])
        for key, value in part['_status'].value_counts().iter_rows():
            tally[key] += value
        tally['palindromic'] += part['_palindromic'].sum()
        tally['palindromic_mismatch'] += part.filter(
            pl.col('_palindromic') & (pl.col('_status') == 'mismatch')
        ).height
        kept.append(part.filter(pl.col('_status') != 'mismatch'))
    counts.update({f'reference_{k}': int(v) for k, v in tally.items()})
    palindromic = tally['palindromic']
    rate = tally['palindromic_mismatch'] / palindromic if palindromic else 0.0
    counts['reference_palindromic_mismatch_rate'] = round(rate, 5)
    assert (
        palindromic < MIN_PALINDROMIC_FOR_CHECK or rate <= MAX_PALINDROMIC_MISMATCH
    ), (
        f'{rate:.1%} of A/T and C/G SNPs do not match GRCh38: positions are not on '
        'the expected build (check source_build in files.csv)'
    )
    if not kept:
        raise ValueError(
            f'no rows left to check against GRCh38; drop counts so far: {counts}'
        )
    return pl.concat(kept).drop('_ref', '_status', '_palindromic')


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

    counts: dict = {}
    with tempfile.TemporaryDirectory(dir=out_root.parent) as tmp:
        workdir = Path(tmp)
        source = read_table(original, workdir)
        columns = resolve_columns(source, row['columns'])
        table, notes = standardise(source, columns, int(row['n_study'] or 0), counts)
        del source
        build = row['source_build']
        if build in LIFTOVER_KNOWN:
            table = lift_to_grch38(table, chains[build], build, counts)
        else:
            assert build == 'GRCh38', f'unsupported source_build {build}'
        if not Path(f'{fasta}.fai').exists():
            pysam.faidx(str(fasta))
        table = check_reference(table, fasta, counts)
        table = table.sort(
            ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele']
        )
        counts['duplicate_variants'] = (
            table.height
            - table.unique(
                ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele']
            ).height
        )
        counts['rows_out'] = table.height
        empty = [c for c in OUTPUT_COLUMNS if table[c].null_count() == table.height]
        plain = workdir / 'formatted.tsv'
        table.select(OUTPUT_COLUMNS).write_csv(plain, separator='\t', null_value='NA')
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
        **counts,
        'formatted_at': datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }
    stats_path.write_text(json.dumps(stats, indent=1) + '\n')
    print(
        f'{row["file_id"]}: {counts["rows_in"]:,} rows in, {counts["rows_out"]:,} out'
    )


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
        fasta_input = batch.read_input_group(base=fasta, fai=f'{fasta}.fai').base
    else:
        fasta_input = batch.read_input(fasta)
    chain_inputs = {build: batch.read_input(path) for build, path in chains.items()}
    submitted = 0
    for row in rows:
        original = f'{args.originals}/{original_name(row)}'
        out_root = f'{args.out}/{formatted_name(row)}'
        if to_path(f'{out_root}.tsv.gz').exists() and not args.force:
            print(f'{row["file_id"]}: exists, skipped')
            continue
        if not to_path(original).exists():
            print(f'{row["file_id"]}: not downloaded, skipped')
            continue
        record = to_path(f'{original}.json')
        if not record.exists():
            print(f'{row["file_id"]}: no download record, skipped')
            continue
        size_gib = int(row['size_bytes'] or 0) / 2**30
        job = batch.new_bash_job(f'format {row["file_id"]}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(8 if size_gib > 1 else 4)
        job.memory('highmem')
        job.storage(f'{int(size_gib * 8) + 20}Gi')
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
    rather than writing a manifest that leaves it out.
    """
    out_dir, stats_dir = as_path(args.out), as_path(args.stats)
    lines, missing_stats, not_formatted = [], [], 0
    for row in rows:
        stats_file = stats_dir / f'{row["file_id"]}.json'
        if not stats_file.exists():
            if (out_dir / f'{formatted_name(row)}.tsv.gz').exists():
                missing_stats.append(row['file_id'])
            else:
                not_formatted += 1
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
    assert lines, f'no stats found in {stats_dir}'
    if not_formatted:
        print(f'{not_formatted} files in files.csv are not formatted yet, left out')
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
        f'<table><tr>{cells}</tr>{body}</table>'
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
    parser.add_argument('--one', help='a single files.csv row as JSON (job mode)')
    parser.add_argument('--original', type=Path, help='job mode: input file')
    parser.add_argument('--out-root', type=Path, help='job mode: output path root')
    parser.add_argument('--stats-file', type=Path, help='job mode: stats JSON path')
    parser.add_argument('--download-record', help='job mode: download record JSON')
    args = parser.parse_args()
    if args.force and not args.only:
        parser.error('--force needs --only, so every file is not redone by mistake')
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
    if args.manifest:
        write_manifest(rows, args)
    elif args.local:
        check_local_args(rows, args, parser)
        run_local(rows, args)
    else:
        run_batch(rows, args)


if __name__ == '__main__':
    main()
