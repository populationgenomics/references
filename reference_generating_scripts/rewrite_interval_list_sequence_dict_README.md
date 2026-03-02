# Purpose

This script is intended to make a new `.interval_list` file by updating an existing interval list with the hashes from the sequence dict of a fasta `.dict` file.

Specifically, it is used to ensure compatibility between CRAM files and interval lists by updating the interval list with the correct sequence dictionary hashes. This is necessary for CRAMs aligned with the masked reference fasta, which has different hashes at the alt contigs than the un-masked fasta used to create the Broad's WGS coverage regions interval list.

## Context

The `Align` pipeline used at the CPG to create aligned CRAM files uses the reference fasta provided by the Broad, sourced from:

```text
gs://gcp-public-data--broad-references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta
```

This masked reference fasta, saved to `hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta`, differs from the original un-masked fasta in `hg38/v0/Homo_sapiens_assembly38.fasta` at the alt contigs (e.g. "chr1_KI270762v1_alt").

Because of this, interval list files derived from the un-masked reference fasta also differ at the alt contigs when comparing the hashes present in the sequence dictionary. E.g. the alt contig hashes in `wgs_coverage_regions.hg38.interval_list` do not match the hashes in the masked fasta.

Some tools (e.g. Picard v3.4.0 CollectWgsMetrics) require an exact match between the sequence dictionary hashes in the CRAM header and the input interval list, even if the interval lists cover the same regions and are the same lengths. For this reason, it is appropriate to re-header an interval list file with the hashes from the reference fasta used by the CRAMs. This enables compatibility between the CRAMs and the interval lists, without modifying the actual regions in the interval file.

## Usage

The default arguments are sufficient for reheadering the Broad's wgs coverage interval list with the hashes from the masked reference fasta. To use the default arguments, simply run:

```bash
python rewrite_interval_list_sequence_dict.py
```

To specify custom arguments, first ensure your inputs and outputs are in the `references.py` file. Pass the appropriate reference strings to the command line arguments as follows:

```bash
python rewrite_interval_list_sequence_dict.py \
    --interval-list-ref <input_interval_list_ref> \
    --fasta-ref <reference_fasta_ref> \
    --out-ref <output_interval_list_ref>
```
