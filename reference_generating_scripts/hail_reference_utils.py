"""
Pieces shared by the Hail Table building scripts in this directory.

Imported as a sibling module: running `python3 reference_generating_scripts/<script>.py`
puts this directory on sys.path.
"""

import gzip

import hail as hl
import hailtop.fs as hfs

# CADD and the SpliceAI hg38 files name contigs without the chr prefix; MT is chrM in GRCh38.
CONTIG_RECODING = {
    **{str(c): f'chr{c}' for c in [*range(1, 23), 'X', 'Y']},
    'MT': 'chrM',
    'M': 'chrM',
}


def read_intervals(bed_path: str) -> list[hl.Interval]:
    """
    The variant-balanced intervals as Python Interval objects, for a partitioned read.

    Parsed here rather than with hl.import_bed, which cannot open a .gz path; the same
    parser ourdna_genomic_atlas uses for this file. BED is 0-based half-open, Hail loci
    are 1-based, so start + 1 with both ends included.
    """
    with hfs.open(bed_path, 'rb') as raw, gzip.open(raw, 'rt') as bed:
        rows = (
            line.split('\t')
            for line in bed
            if line.strip() and not line.startswith(('#', 'track'))
        )
        return [
            hl.Interval(
                hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38'),
                hl.Locus(chrom, int(end), reference_genome='GRCh38'),
                includes_start=True,
                includes_end=True,
            )
            for chrom, start, end, *_ in rows
        ]


def refuse_existing(out: str, overwrite: bool) -> None:
    """
    Fail before any work if `out` exists and replacing it was not asked for.

    Hail's write(overwrite=True) deletes the target before writing, so a re-run that
    dies mid-write leaves every consumer of a published table reading a corrupt one.
    """
    if hfs.exists(out) and not overwrite:
        raise FileExistsError(f'{out} exists; pass --overwrite to replace it')


def write_on_intervals(
    ht: hl.Table,
    intervals: list[hl.Interval],
    out: str,
    overwrite: bool,
    source: dict[str, str],
) -> None:
    """
    Checkpoint `ht`, re-read it partitioned on `intervals`, record `source` in the
    globals and write to `out`, then fail if the re-read dropped any row.

    The checkpoint lets the sort or merge pick its own layout; the re-read decides the
    layout that is written. `_intervals` is a private Hail argument, the standard idiom
    for read-time partitioning, and it also filters to the intervals, so the row count
    is compared before and after (both from table metadata, no scan).
    """
    tmp = hl.utils.new_temp_file(out.rstrip('/').rsplit('/', 1)[-1], 'ht')
    ht.checkpoint(tmp)
    n_in = hl.read_table(tmp).count()
    ht = hl.read_table(tmp, _intervals=intervals)
    ht = ht.annotate_globals(source=hl.struct(**source))
    ht.write(out, overwrite=overwrite)
    ht = hl.read_table(out)
    n_out = ht.count()
    if n_out != n_in:
        raise ValueError(
            f'{n_in - n_out:,} rows fell outside the intervals: {n_in:,} in, {n_out:,} written to {out}'
        )
    ht.describe()
    print(f'{n_out:,} rows in {ht.n_partitions()} partitions at {out}')
