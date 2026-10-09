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


REPO='illumina-ai/SpliceAI2'
REVISION=${3:-"55948bae9b3638009178aadbf9d9a51c1017814f"}

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

hf download "$REPO" --revision "$REVISION" --local-dir "$WORKDIR/SpliceAI2"

gcloud storage rsync -r --exclude='^\.cache/' "$WORKDIR/SpliceAI2" "${DEST%/}"
