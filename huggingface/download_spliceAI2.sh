#!/usr/bin/env bash
# Download the SpliceAI2 model checkpoints (gated, non-commercial licence) from
# Hugging Face and upload them to GCS.
#
# Requires: hf CLI, gcloud. The token must belong to an account that has
# accepted the SpliceAI2 terms on huggingface.co.
#
# Usage: download_spliceAI2.sh HF_TOKEN gs://bucket/references/spliceai2/<revision>

set -euo pipefail


HF_TOKEN="${1:?Usage: $0 HF_TOKEN gs://destination/prefix}"
DEST="${2:?Usage: $0 HF_TOKEN gs://destination/prefix}"
export HF_TOKEN

REPO='illumina-ai/SpliceAI2-data'
REVISION=${3:-"cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"}

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

hf download "$REPO" \
  --revision "$REVISION" \
  --repo-type dataset \
  --include "precomputed_scores_v2.0/*" \
  --local-dir "$WORKDIR/SpliceAI2"

gcloud storage rsync -r --exclude='^\.cache/' "$WORKDIR/SpliceAI2" "${DEST%/}"
