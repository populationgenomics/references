"""
Builds a STAR reference genome index through Hail batch.
"""
from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, reference_path, get_batch
from textwrap import dedent

# Get config and parameters
config = get_config()
# IS_TEST = config['workflow']['access_level'] == 'test'
TEST_BUCKET = config['storage']['common']['tmp']
CPU = int(config['workflow'].get('n_cpu', 8))
MEMORY = config['workflow'].get('memory', 'highmem')
STORAGE = config['workflow'].get('storage', '150Gi')
STAR_VERSION = str(config['star']['version'])
SJDB_OVERHANG = int(config['star'].get('sjdb_overhang', 100))
GENCODE_VERSION = str(config['star'].get('gencode_version', 44))
GENCODE_BASE_URL = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_VERSION}'
GENCODE_GTF_BASENAME = f'gencode.v{GENCODE_VERSION}.primary_assembly.annotation.gtf'
GENCODE_GTF_URL = f'{GENCODE_BASE_URL}/{GENCODE_GTF_BASENAME}.gz'
# GENCODE_FASTA_BASENAME = 'GRCh38.primary_assembly.genome.fa'
# GENCODE_FASTA_URL = f'{GENCODE_BASE_URL}/{GENCODE_FASTA_BASENAME}.gz'

b = get_batch('Build Star Reference', default_image=config['workflow']['driver_image'])

# Job to get and subset reference files
get_ref_j = b.new_job('get-ref-files')
get_ref_j.image(image_path('samtools'))
get_ref_j.cpu(CPU)
get_ref_j.memory(MEMORY)
get_ref_j.storage(STORAGE)

TMPDIR = '$BATCH_TMPDIR'
TMP_DL_DIR = f'{TMPDIR}/dl'
TMP_FASTA_DIR = f'{TMPDIR}/fasta'

major_chromosomes = ' '.join([f'chr{str(i)}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

get_ref_j.declare_resource_group(
    ref_files={
        'fa': 'hg38.fa',
        'fa_idx': 'hg38.fa.fai',
        'gtf': 'hg38.gtf',
    },
)

ref_fasta = b.read_input_group(
    fa=str(reference_path('broad/ref_fasta')),
)

cmd = f"""\
    # Create directories
    mkdir -p {TMP_DL_DIR}
    mkdir -p {TMP_FASTA_DIR}

    # Download GTF and strip down to just the major chromosomes
    cd {TMP_DL_DIR}
    wget {GENCODE_GTF_URL}
    gunzip {GENCODE_GTF_BASENAME}.gz
    awk -v FS="\\t" '$0~/^#/ || $1~/^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)/' {GENCODE_GTF_BASENAME} > {get_ref_j.ref_files.gtf}

    # Strip FASTA down to just the major chromosomes
    cd {TMP_FASTA_DIR}
    samtools faidx {ref_fasta.fa}
    samtools faidx {ref_fasta.fa} {major_chromosomes} > hg38.fa
    samtools faidx hg38.fa
    mv hg38.fa {get_ref_j.ref_files.fa}
    mv hg38.fa.fai {get_ref_j.ref_files.fa_idx}
    """
cmd = dedent(cmd)

get_ref_j.command(cmd)

# Write outputs
OUT_REF_DIR = to_path(TEST_BUCKET) / 'references' / 'star' / 'hg38'
b.write_output(get_ref_j.ref_files.fa, str(OUT_REF_DIR / 'hg38.fa'))
b.write_output(get_ref_j.ref_files.fa_idx, str(OUT_REF_DIR / 'hg38.fa.fai'))
b.write_output(get_ref_j.ref_files.gtf, str(OUT_REF_DIR / 'hg38.gtf'))

# Job to build reference
j = b.new_job('build-star-reference')
j.image(image_path('star'))
j.cpu(CPU)
j.memory(MEMORY)
j.storage(STORAGE)

TMP_MKREF_DIR = f'{TMPDIR}/mkref'
TMP_GENOME_DIR = f'{TMP_MKREF_DIR}/hg38'
OUT_GENOME_DIR = to_path(TEST_BUCKET) / 'references' / 'star' / str(STAR_VERSION) / 'hg38'

star_ref_files = {
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
j.declare_resource_group(star_ref=star_ref_files)


# Build command
cmd = f"""\
    # Create directories
    mkdir -p {TMP_GENOME_DIR}

    # Build reference
    cd {TMP_MKREF_DIR}
    STAR \
        --runThreadN {str(CPU)} \
        --runMode genomeGenerate \
        --genomeDir {TMP_GENOME_DIR} \
        --genomeFastaFiles {get_ref_j.ref_files.fa} \
        --sjdbGTFfile {get_ref_j.ref_files.gtf} \
        --sjdbOverhang {str(SJDB_OVERHANG)}
    """
cmd = dedent(cmd)

# Move reference files to resource group
for key, file in star_ref_files.items():
    cmd += f'mv {TMP_GENOME_DIR}/{file} {j.star_ref[key]}\n'

j.command(cmd)

# Write outputs
# Have to do this one-by-one since the output files don't have a simple naming convention
for key, file in star_ref_files.items():
    b.write_output(j.star_ref[key], str(OUT_GENOME_DIR / file))

b.run(wait=False)
