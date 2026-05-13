#!/usr/bin/env bash

# One-shot cleanup for the liftover_37_to_38 chain filename typo. PR #96 landed
# the chain at the typo'd path (missing a `.` before `over`):
#     gs://cpg-common-main/references/liftover/grch37_to_grch38over.chain.gz
# The companion references-repo PR corrects the Source `dst` to
# `grch37_to_grch38.over.chain.gz`; CI auto-curls the 228 KB chain to the
# corrected path on merge. This script removes the typo'd blob.
#
# analysis-runner \
#     --dataset common \
#     --access-level full \
#     --description 'Rename liftover_37_to_38 chain typo' \
#     --output-dir references/liftover-rename \
#     rename_liftover_chain_typo.sh
#
# Idempotent: handles old-only (rename), both-present (delete old), and
# new-only (no-op) states. Errors if neither blob is present.

set -euo pipefail

OLD='gs://cpg-common-main/references/liftover/grch37_to_grch38over.chain.gz'
NEW='gs://cpg-common-main/references/liftover/grch37_to_grch38.over.chain.gz'

if gcloud storage objects describe "$OLD" >/dev/null 2>&1; then
    old_exists=1
else
    old_exists=0
fi

if gcloud storage objects describe "$NEW" >/dev/null 2>&1; then
    new_exists=1
else
    new_exists=0
fi

if (( old_exists == 0 && new_exists == 1 )); then
    echo "Typo blob already absent; corrected blob at $NEW. Nothing to do."
    exit 0
fi

if (( old_exists == 0 && new_exists == 0 )); then
    echo "Neither $OLD nor $NEW exists. Has the references-repo PR correcting liftover_37_to_38 dst merged and CI run?" >&2
    exit 1
fi

if (( old_exists == 1 && new_exists == 0 )); then
    gcloud storage mv "$OLD" "$NEW"
else
    gcloud storage rm "$OLD"
fi
