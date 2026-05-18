#!/usr/bin/env python3

"""
Parse an Illumina BPM (Bead Pool Manifest) and emit a BED6+3 of biallelic SNV
sites, polarised REF/ALT against a reference assembly fasta.

Output columns (tab-separated):
    chrom  chromStart  chromEnd  name  score(=0)  strand(=RefStrand)
    REF  ALT  SourceStrand

Drops:
    - unmapped probes (Chr == '0' or MapInfo == 0)
    - PAR probes (Chr == 'XY')
    - non-SNV probes (SNP alleles outside {A,C,G,T}; excludes [I/D]/[D/I])
    - probes with Unknown RefStrand
    - sites where neither manifest allele matches the reference base

Chromosome mapping (manifest -> assembly): 'MT' -> 'chrM'; otherwise 'chr' is
prepended.

Dependencies:
    pip install pyfaidx
    pip install git+https://github.com/Illumina/BeadArrayFiles.git
"""

import sys
from argparse import ArgumentParser

from IlluminaBeadArrayFiles import BeadPoolManifest, RefStrand, SourceStrand
from pyfaidx import Fasta


COMP = str.maketrans('ACGT', 'TGCA')
ACGT = set('ACGT')
OUTPUT_SUFFIX = '-biallelic-snps.bed'


def chrom_to_assembly(chrom: str) -> str | None:
    if chrom in ('0', 'XY', ''):
        return None
    if chrom == 'MT':
        return 'chrM'
    return f'chr{chrom}'


def ref_strand_symbol(rs: int) -> str | None:
    if rs == RefStrand.Plus:
        return '+'
    if rs == RefStrand.Minus:
        return '-'
    return None


def source_strand_symbol(ss: int) -> str:
    if ss == SourceStrand.Forward:
        return 'F'
    if ss == SourceStrand.Reverse:
        return 'R'
    return '.'


def parse_snp(snp: str) -> tuple[str, str] | None:
    """Parse SNP designation '[X/Y]' into (X, Y); None if malformed or non-SNV."""
    if len(snp) != 5 or snp[0] != '[' or snp[2] != '/' or snp[4] != ']':
        return None
    a, b = snp[1], snp[3]
    if a not in ACGT or b not in ACGT:
        return None
    return a, b


def main(bpm_path: str, fasta_path: str, output_path: str) -> None:
    manifest = BeadPoolManifest(bpm_path)
    fasta = Fasta(fasta_path, sequence_always_upper=True, as_raw=True)

    counts = {
        'total': manifest.num_loci,
        'unmapped': 0,
        'par': 0,
        'non_snv': 0,
        'unknown_strand': 0,
        'missing_contig': 0,
        'ref_mismatch': 0,
        'written': 0,
    }

    with open(output_path, 'w') as out:
        for i in range(manifest.num_loci):
            chrom_in = manifest.chroms[i]
            pos = manifest.map_infos[i]

            if chrom_in == '0' or pos == 0:
                counts['unmapped'] += 1
                continue
            if chrom_in == 'XY':
                counts['par'] += 1
                continue

            alleles = parse_snp(manifest.snps[i])
            if alleles is None:
                counts['non_snv'] += 1
                continue

            strand = ref_strand_symbol(manifest.ref_strands[i])
            if strand is None:
                counts['unknown_strand'] += 1
                continue

            a1, a2 = alleles
            if strand == '-':
                a1 = a1.translate(COMP)
                a2 = a2.translate(COMP)

            chrom = chrom_to_assembly(chrom_in)
            if chrom is None or chrom not in fasta:
                counts['missing_contig'] += 1
                continue

            ref_base = fasta[chrom][pos - 1]
            if ref_base == a1:
                ref, alt = a1, a2
            elif ref_base == a2:
                ref, alt = a2, a1
            else:
                counts['ref_mismatch'] += 1
                continue

            out.write(
                f'{chrom}\t{pos - 1}\t{pos}\t{manifest.names[i]}\t0\t{strand}\t'
                f'{ref}\t{alt}\t{source_strand_symbol(manifest.source_strands[i])}\n'
            )
            counts['written'] += 1

    width = max(len(k) for k in counts)
    for k, v in counts.items():
        print(f'{k:<{width}}  {v}', file=sys.stderr)


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--bpm', required=True, help='Illumina BPM manifest file')
    parser.add_argument(
        '--fasta',
        required=True,
        help='Reference fasta (must be indexed; .fai alongside)',
    )
    parser.add_argument(
        '--output',
        required=True,
        help=f"Output BED path; must end with '{OUTPUT_SUFFIX}'",
    )
    args = parser.parse_args()
    if not args.output.endswith(OUTPUT_SUFFIX):
        parser.error(f"--output must end with '{OUTPUT_SUFFIX}'")
    main(args.bpm, args.fasta, args.output)


if __name__ == '__main__':
    cli_main()
