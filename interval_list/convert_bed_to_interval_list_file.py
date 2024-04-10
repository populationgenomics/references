import click

from cpg_utils.hail_batch import command, get_batch, image_path, reference_path


def get_ref_files(bed_ref: str, sd_ref: str):
    """
    Converts a BED file to a .interval_list file using Picard BedToIntervalList.
    Saves the .interval_list file to the reference directory with the specified outfile name.
    """  
    bed_file = reference_path(bed_ref)
    sd_file = reference_path(sd_ref)
    
    b = get_batch()
    j = b.new_job(f'Convert {bed_file} to .interval_list file')
    j.image(image_path('picard'))
    
    cmd = f"""
    mkdir $BATCH_TMPDIR/out
    
    picard BedToIntervalList \
    -I {b.read_input(str(bed_file))} \
    -O $BATCH_TMPDIR/out/outfile.interval_list \
    -SD {b.read_input(str(sd_file))}
    
    ln $BATCH_TMPDIR/out/outfile.interval_list {j['interval_list']}
    """
    
    j.command(command(cmd))
    b.write_output(j['interval_list'], reference_path(f'{bed_ref}_interval_list'))
    b.run()
    
    
@click.command()
@click.option('--bed-ref', required=True, help='String identifier for the BED file from the references.', default='hg38_telomeres_and_centromeres')
@click.option('--sd-ref', required=True, help='String identifier for the sequence dictionary file from the references.', default='broad/genome_calling_interval_lists')
def main(bed_ref, sd_ref):
    """
    Converts a BED file to a .interval_list file using Picard BedToIntervalList.
    """
    get_ref_files(bed_ref, sd_ref)
