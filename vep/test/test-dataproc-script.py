import sys
import hail as hl

mt_path = sys.argv[1]
out_mt_path = sys.argv[2]

hl.init(default_reference='GRCh38')

mt = hl.read_matrix_table(mt_path)

# Filter to biallelic loci only
mt = mt.filter_rows(hl.len(mt.alleles) == 2)
mt = mt.filter_rows(mt.alleles[1] != '*')

mt = hl.vep(
    mt, 
    # Because we are not starting the cluster with `--vep` and instead passing custom 
    # startup script with `--init`, Hail would not set VEP_CONFIG_URI for us. Thus, 
    # we need to provide config as a function parameter explicitly:
    config='file:///vep_data/vep-gcloud.json'
)
mt.write(out_mt_path, overwrite=True)
