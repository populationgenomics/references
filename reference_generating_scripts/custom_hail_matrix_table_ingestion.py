from cpg_utils import hail_batch
import hail as hl
import argparse

hail_batch.init_batch()

# Parse command line arguments
parser = argparse.ArgumentParser(description='Import TSV file into Hail Table format')
parser.add_argument('input_file', type=str, help='Path to the input TSV file (can be gzipped)')
parser.add_argument('--output', '-o', type=str, default=None,
                    help='Output path for Hail table (default: input_file.ht)')
parser.add_argument('--reference-genome', '-r', type=str, default='GRCh38',
                    choices=['GRCh37', 'GRCh38'],
                    help='Reference genome to use (default: GRCh38)')
args = parser.parse_args()

# Define paths
tsv_path = args.input_file
output_path = args.output if args.output else tsv_path.replace('.tsv.gz', '.ht').replace('.tsv', '.ht')

# 1. Define the input types for the initial table import
# We import chrom/pos/ref/alt as strings/ints first to transform them
input_types = {
    '#CHROM': hl.tstr,
    'POS': hl.tint32,
    'REF': hl.tstr,
    'ALT': hl.tstr,
    'avis': hl.tfloat64
}

# 2. Load the TSV
ht = hl.import_table(
    tsv_path,
    types=input_types,
    delimiter='\t',
    force_bgz=True
)

# 3. Transform to standard Hail genomic format
# We create a 'locus' object and an 'alleles' array
ht = ht.transmute(
    locus=hl.locus(ht['#CHROM'], ht.POS, reference_genome=args.reference_genome),
    alleles=[ht.REF, ht.ALT]
)

# 4. Key the table by locus and alleles
ht = ht.key_by('locus', 'alleles')

ht.describe()
ht.show(5)

# Write the table to disk in Hail format for later use
ht.write(output_path, overwrite=True)
print(f'Table successfully written to: {output_path}')
