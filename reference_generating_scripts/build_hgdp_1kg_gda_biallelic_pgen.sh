#!/usr/bin/env bash
#
# Provenance for the HGDP+1KG v1 PLINK2 trio
# (hgdp-1kg-v1-GDA_8v1_0_D2_biallelic_snps.{pgen,pvar,psam}) staged under
# gs://cpg-common-*/references/genotype_array_reference_data/ via the
# `genotype_array_reference_data` Source in references.py.
#
# Input: HGDP+1KG VCF already subset to GDA-8v1-0_D2 biallelic SNV sites.
# Auxiliary inputs:
#   - hgdp-1kg-cpg.psam        — sex assignments used by --update-sex
#   - hgdp-1kg-cpg-remove.txt  — sample IDs to drop via --remove
#
# Output: PLINK2 .pgen/.pvar/.psam trio with:
#   - SNVs only (--snps-only)
#   - chrM-style contigs (--output-chr chrM)
#   - variant IDs set to CHROM:POS:REF:ALT (--set-all-var-ids)
#   - per-variant call-rate == 100% (--geno 0)
#   - HWE filter, mid-p, keep-fewhet (--hwe 1e-5 0.001 midp keep-fewhet)
#
# Requirements: plink2 on PATH.

set -euo pipefail

VCF=''
PSAM_SEX=''
REMOVE=''
OUT_PREFIX='hgdp-1kg-v1-GDA_8v1_0_D2_biallelic_snps'

usage() {
    cat <<EOF
Usage:
  $0 --vcf VCF --update-sex PSAM --remove REMOVE [--out PREFIX]

  --vcf         Input VCF (.vcf.bgz) subset to GDA biallelic SNV sites
  --update-sex  PLINK2 .psam supplying sex assignments
  --remove      Newline-delimited sample IDs to drop
  --out         Output prefix (default: ${OUT_PREFIX})
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        VCF="$2"; shift 2 ;;
        --update-sex) PSAM_SEX="$2"; shift 2 ;;
        --remove)     REMOVE="$2"; shift 2 ;;
        --out)        OUT_PREFIX="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$VCF" || -z "$PSAM_SEX" || -z "$REMOVE" ]]; then
    echo "ERROR: --vcf, --update-sex, and --remove are required." >&2
    usage >&2
    exit 2
fi

plink2 \
    --vcf "$VCF" \
    --make-pgen \
    --out "$OUT_PREFIX" \
    --update-sex "$PSAM_SEX" \
    --snps-only \
    --output-chr chrM \
    --set-all-var-ids '@:#:$r:$a' \
    --remove "$REMOVE" \
    --geno 0 \
    --hwe 1e-5 0.001 midp keep-fewhet
