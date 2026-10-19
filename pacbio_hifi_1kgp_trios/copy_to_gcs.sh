#!/usr/bin/env bash
# Copy the 1kGP PacBio HiFi trio uBAMs (+ .pbi) from the public 1000g-ont S3 bucket
# into GCS, mirroring the pacbio_hifi_1kgp_trios block in popgen:references references.py.
# Run on a GCE VM: each file streams through this machine.
#
# Usage: copy_to_gcs.sh [gs://destination/prefix]

set -euo pipefail

SRC='https://s3.amazonaws.com/1000g-ont/1KGP_PacBio_WGS/Raw_Data'
DEST="${1:-gs://cpg-common-main/references/1kgp/pacbio_hifi_trios}"
DEST="${DEST%/}"

SAMPLES=(
  # m039 (MXL): child NA19796, father NA19795, mother NA19794
  GM19796-W-0 GM19795-W-0 GM19794-W-0
  # PEL005 (PEL): child HG01579, father HG01577, mother HG01578
  HG01579-W-0 HG01577-W-0 HG01578-C-0
)

copy_sample() {
  local sample=$1 ext path
  for ext in bam bam.pbi; do
    path="${sample}/${sample}.hifi_reads.${ext}"
    if gcloud storage objects describe "${DEST}/${path}" &>/dev/null; then
      echo "exists, skipping: ${DEST}/${path}"
      continue
    fi
    echo "copying ${path}"
    # a failed stream still finalises a partial object, so remove it (as _curl_pipe does)
    curl -A "Mozilla/5.0" -fsSL "${SRC}/${path}" | gcloud storage cp - "${DEST}/${path}" \
      || { gcloud storage rm "${DEST}/${path}" || true; echo "FAILED: ${path}"; return 1; }
  done
}

# one background job per sample; wait on each so one failure doesn't orphan the rest
pids=()
for sample in "${SAMPLES[@]}"; do
  copy_sample "$sample" &
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  wait "$pid" || failed=1
done
exit "$failed"
