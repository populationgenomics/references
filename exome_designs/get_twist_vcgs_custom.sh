#!/usr/bin/env bash
# Fetch the Twist Comprehensive Exome + VCGS custom content target bed (hg38).
# Output goes alongside this script; upload to cpg-common-test afterwards.
set -euo pipefail
URL='https://www.twistbioscience.com/content/dam/twistbioscience/resources/2021-10/VCGS_Exome_Covered_Targets_hg38%2040.1MB%20.bed'
curl -L "$URL" -o Twist_VCGS_Exome_Covered_Targets_hg38.bed
