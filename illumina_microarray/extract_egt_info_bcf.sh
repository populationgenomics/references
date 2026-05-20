#!/usr/bin/env bash
#
# Re-transcribe the per-variant INFO payload of an Illumina EGT cluster file
# into a sample-less BCF using `bcftools +gtc2vcf` with no `--gtcs`. The output
# is the "expected" cluster geometry (GenTrain_Score, Cluster_Sep, N_AA/AB/BB,
# meanTHETA_*, devTHETA_*, meanR_*, devR_*) keyed on (CHROM, POS, REF, ALT) —
# the reference half of any observed-vs-EGT comparison.
#
# Run locally (or inside the bcftools:1.23-1 image) against the canonical
# BPM / EGT / fasta references staged in cpg-common-main. Staging the output
# to cpg-common-test and promoting to cpg-common-main are handled out-of-band
# (see illumina_microarray/copy_egt_info_bcf_to_main.py for the test->main
# promote step).
#
# Requirements: bcftools (with the gtc2vcf plugin) on PATH.
#
# Example:
#   ./extract_egt_info_bcf.sh \
#     --bpm   /path/to/GDA-8v1-0_D2.bpm \
#     --egt   /path/to/GDA-8v1-0_D1_ClusterFile.egt \
#     --fasta /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#     --out-dir ./out

set -euo pipefail

OUT_NAME='GDA-8v1-0_D1_ClusterFile_info.bcf'

BPM=''
EGT=''
FASTA=''
OUT_DIR='./out'

usage() {
    cat <<EOF
Usage:
  $0 --bpm BPM --egt EGT --fasta FASTA [--out-dir DIR]

  --bpm      Path to the Illumina BPM manifest (.bpm)
  --egt      Path to the Illumina EGT cluster file (.egt)
  --fasta    Path to the reference fasta (.fna); .fai must sit alongside
  --out-dir  Local output directory (default: ./out)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bpm)     BPM="$2"; shift 2 ;;
        --egt)     EGT="$2"; shift 2 ;;
        --fasta)   FASTA="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

mkdir -p "$OUT_DIR"
OUT_BCF="${OUT_DIR}/${OUT_NAME}"
OUT_CSI="${OUT_BCF}.csi"
METADATA_TSV="${OUT_DIR}/${OUT_NAME%.bcf}.metadata.tsv"

if [[ -z "$BPM" || -z "$EGT" || -z "$FASTA" ]]; then
    echo "ERROR: --bpm, --egt, and --fasta are required." >&2
    usage >&2
    exit 2
fi

for f in "$BPM" "$EGT" "$FASTA" "${FASTA}.fai"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing input: $f" >&2
        exit 1
    fi
done

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools not on PATH." >&2
    exit 1
fi

if ! bcftools plugin -l 2>/dev/null | grep -q gtc2vcf; then
    echo "ERROR: bcftools is missing the gtc2vcf plugin. Use the bcftools:1.23-1 image." >&2
    exit 1
fi

# Mirrors src/popgen_genotyping/jobs/gtc_to_bcfs_job.py — same flags + norm
# + sort — but with no `--gtcs`, so the resulting record set has only INFO
# (and an empty FORMAT column, which `view -G` strips defensively).
echo "Running bcftools +gtc2vcf (no GTCs) -> norm -> sort -> view -G | bcf ..."
bcftools +gtc2vcf \
    --no-version \
    --do-not-check-bpm \
    --bpm "$BPM" \
    --egt "$EGT" \
    --fasta-ref "$FASTA" \
    --extra "$METADATA_TSV" \
  | bcftools norm -m -both --no-version -c x -f "$FASTA" \
  | bcftools sort -T "$(mktemp -d)" \
  | bcftools view -G -O b -o "$OUT_BCF" --write-index=csi

echo "Wrote ${OUT_BCF}"
echo "Wrote ${OUT_CSI}"
echo "Wrote ${METADATA_TSV}"

if ! bcftools view -h "$OUT_BCF" | grep -q 'ID=GenTrain_Score'; then
    echo "WARNING: INFO/GenTrain_Score not declared in header — verify gtc2vcf output." >&2
fi
