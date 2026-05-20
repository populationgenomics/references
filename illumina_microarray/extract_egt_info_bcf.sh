#!/usr/bin/env bash
#
# Re-transcribe the per-variant INFO payload of an Illumina EGT cluster file
# into a sample-less BCF using `bcftools +gtc2vcf` with no `--gtcs`. The output
# is the "expected" cluster geometry (GenTrain_Score, Cluster_Sep, N_AA/AB/BB,
# meanTHETA_*, devTHETA_*, meanR_*, devR_*) keyed on (CHROM, POS, REF, ALT) —
# the reference half of any observed-vs-EGT comparison.
#
# Run locally (or inside the bcftools:1.23-1 image) against the canonical
# BPM / EGT / fasta references staged in cpg-common-main. Optionally uploads
# the resulting BCF + index to the cpg-common-test bucket; the references-repo
# CI then rsyncs from there into cpg-common-main on merge.
#
# Requirements: bcftools (with the gtc2vcf plugin) on PATH; gcloud SDK on PATH
# if --upload is passed.
#
# Example — local run, no upload:
#   ./extract_egt_info_bcf.sh \
#     --bpm  /path/to/GDA-8v1-0_D2.bpm \
#     --egt  /path/to/GDA-8v1-0_D1_ClusterFile.egt \
#     --fasta /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#     --out-dir ./out
#
# Example — fetch inputs from cpg-common-main, run, upload to cpg-common-test:
#   gcloud storage cp \
#     gs://cpg-common-main/references/illumina_microarray/{GDA-8v1-0_D2.bpm,GDA-8v1-0_D1_ClusterFile.egt,GCA_000001405.15_GRCh38_no_alt_analysis_set.fna,GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai} \
#     ./inputs/
#   ./extract_egt_info_bcf.sh \
#     --bpm   ./inputs/GDA-8v1-0_D2.bpm \
#     --egt   ./inputs/GDA-8v1-0_D1_ClusterFile.egt \
#     --fasta ./inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#     --out-dir ./out \
#     --upload

set -euo pipefail

OUT_NAME='GDA-8v1-0_D1_egt_info.bcf'
TEST_BUCKET_PREFIX='gs://cpg-common-test/references/illumina_microarray'
BILLING_PROJECT='cpg-common'

BPM=''
EGT=''
FASTA=''
OUT_DIR='./out'
UPLOAD=0
UPLOAD_ONLY=0

usage() {
    cat <<EOF
Usage:
  build:        $0 --bpm BPM --egt EGT --fasta FASTA [--out-dir DIR]
  build+upload: $0 --bpm BPM --egt EGT --fasta FASTA [--out-dir DIR] --upload
  upload only:  $0 --upload-only [--out-dir DIR]

  --bpm          Path to the Illumina BPM manifest (.bpm)
  --egt          Path to the Illumina EGT cluster file (.egt)
  --fasta        Path to the reference fasta (.fna); .fai must sit alongside
  --out-dir      Local output / source directory (default: ./out)
  --upload       After build, copy ${OUT_NAME}{,.csi} to
                 ${TEST_BUCKET_PREFIX}/
  --upload-only  Skip the bcftools build; just upload existing
                 \${OUT_DIR}/${OUT_NAME}{,.csi} to the test bucket. Use this
                 when the build was run inside a docker container and the
                 upload runs from the host (where gcloud is on PATH).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bpm)         BPM="$2"; shift 2 ;;
        --egt)         EGT="$2"; shift 2 ;;
        --fasta)       FASTA="$2"; shift 2 ;;
        --out-dir)     OUT_DIR="$2"; shift 2 ;;
        --upload)      UPLOAD=1; shift ;;
        --upload-only) UPLOAD_ONLY=1; UPLOAD=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

mkdir -p "$OUT_DIR"
OUT_BCF="${OUT_DIR}/${OUT_NAME}"
OUT_CSI="${OUT_BCF}.csi"
METADATA_TSV="${OUT_DIR}/${OUT_NAME%.bcf}.metadata.tsv"

if [[ "$UPLOAD_ONLY" -eq 0 ]]; then
    if [[ -z "$BPM" || -z "$EGT" || -z "$FASTA" ]]; then
        echo "ERROR: --bpm, --egt, and --fasta are required (or pass --upload-only)." >&2
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
fi

if [[ "$UPLOAD" -eq 1 ]]; then
    if ! command -v gcloud >/dev/null 2>&1; then
        echo "ERROR: --upload requested but gcloud not on PATH." >&2
        exit 1
    fi
    for f in "$OUT_BCF" "$OUT_CSI"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: expected artifact missing for upload: $f" >&2
            exit 1
        fi
    done
    echo "Uploading to ${TEST_BUCKET_PREFIX}/ ..."
    gcloud --billing-project "$BILLING_PROJECT" storage cp \
        "$OUT_BCF" "$OUT_CSI" \
        "${TEST_BUCKET_PREFIX}/"
    echo "Upload complete."
fi
