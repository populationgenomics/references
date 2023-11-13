"""
Builds a STAR reference genome index through Hail batch.
"""
import hailtop.batch as hb
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, image_path, reference_path
from textwrap import dedent

# Get config and parameters
config = get_config()
BILLING_PROJECT = config['hail']['billing_project']
DEFAULT_IMAGE = config['workflow']['driver_image']
TMP_BUCKET = config['storage']['common']['tmp']
TEST_BUCKET = config['storage']['common']['test']
CPU = int(config['workflow'].get('n_cpu', 8))
MEMORY = config['workflow'].get('memory', 'standard')
STORAGE = config['workflow'].get('storage', '150Gi')
STAR_VERSION = str(config['star']['version'])
SJDB_OVERHANG = int(config['star'].get('sjdb_overhang', 100))
GENCODE_VERSION = str(config['star'].get('gencode_version', 44))
GENCODE_BASE_URL = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_VERSION}'
GENCODE_GTF_BASENAME = f'gencode.v{GENCODE_VERSION}.primary_assembly.annotation.gtf'
GENCODE_GTF_URL = f'{GENCODE_BASE_URL}/{GENCODE_GTF_BASENAME}.gz'
# GENCODE_FASTA_BASENAME = 'GRCh38.primary_assembly.genome.fa'
# GENCODE_FASTA_URL = f'{GENCODE_BASE_URL}/{GENCODE_FASTA_BASENAME}.gz'

sb = hb.ServiceBackend(billing_project=BILLING_PROJECT, remote_tmpdir=remote_tmpdir())
b = hb.Batch(backend=sb, default_image=DEFAULT_IMAGE)

j = b.new_job('build-star-reference')
j.image(image_path('star'))
j.cpu(CPU)
j.memory(MEMORY)
j.storage(STORAGE)

TMPDIR = Path("$BATCH_TMPDIR")
TMP_DL_DIR = TMPDIR / 'dl'
TMP_FASTA_DIR = TMPDIR / 'fasta'
TMP_MKREF_DIR = TMPDIR / 'mkref'
TMP_GENOME_DIR = TMP_MKREF_DIR / 'hg38'
OUT_GENOME_DIR = Path(TEST_BUCKET) / 'references' / 'star' / 'hg38'

star_ref_file_basenames = {
    'chr_len': 'chrLength.txt',
    'chr_name_len': 'chrNameLength.txt',
    'chr_name': 'chrName.txt',
    'chr_start': 'chrStart.txt',
    'exon_ge_tr_info': 'exonGeTrInfo.tab',
    'exon_info': 'exonInfo.tab',
    'gene_info': 'geneInfo.tab',
    'genome': 'Genome',
    'genome_params': 'genomeParameters.txt',
    'sa': 'SA',
    'sa_idx': 'SAindex',
    'sjdb_info': 'sjdbInfo.txt',
    'sjdb_list_gtf': 'sjdbList.fromGTF.out.tab',
    'sjdb_list': 'sjdbList.out.tab',
    'transcript_info': 'transcriptInfo.tab',
}
star_ref_files = {
    key: TMP_MKREF_DIR / file
    for key, file in star_ref_file_basenames.items()
}
j.declare_resource_group(star_ref=star_ref_files)

major_chromosomes = ' '.join([f'chr{str(i)}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

# Build command
cmd = f"""\
    # Download GTF and strip down to just the major chromosomes
    mkdir -p {TMP_DL_DIR}
    cd {TMP_DL_DIR}
    wget {GENCODE_GTF_URL}
    gunzip {GENCODE_GTF_BASENAME}.gz
    awk -v FS="\\t" '$0~/^#/ || $1~/^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)/' {GENCODE_GTF_BASENAME} > hg38.gtf
    # Strip FASTA down to just the major chromosomes
    mkdir -p {TMP_FASTA_DIR}
    cd {TMP_FASTA_DIR}
    samtools faidx {reference_path('broad/fasta')} {major_chromosomes} > hg38.fa
    samtools faidx hg38.fa
    # Build reference
    mkdir -p {TMP_MKREF_DIR}
    cd {TMPDIR}
    STAR
        --runThreadN {str(CPU)}
        --runMode genomeGenerate
        --genomeDir hg38
        --genomeFastaFiles {TMP_FASTA_DIR / 'hg38.fa'}
        --sjdbGTFfile {TMP_DL_DIR / 'hg38.gtf'}
        --sjdbOverhang {str(SJDB_OVERHANG)}
    """
cmd = dedent(cmd)

j.command(cmd)

# Write outputs
# Have to do this one-by-one since the output files don't have a simple naming convention
for key, file in star_ref_files.items():
    j.write_output(file, str(OUT_GENOME_DIR / star_ref_file_basenames[key]))

b.run(wait=False)
