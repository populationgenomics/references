import click

from cpg_utils.hail_batch import command, get_batch, image_path, reference_path


def get_ref_files(bed_ref: str, sd_ref: str, outfile_path: str) -> None:
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
    b.write_output(j['interval_list'], outfile_path)
    b.run()
    
    
@click.command()
@click.option('--bed-ref', required=True, help='String identifier for the BED file from the references.', default='hg38_telomeres_and_centromeres')
@click.option('--sd-ref', required=True, help='String identifier for the sequence dictionary file from the references.', default='broad/genome_calling_interval_lists')
@click.option('--out-ref', required=True, help='Reference path to save the output .interval_list file.', default='hg38_telomeres_and_centromeres_intervals/interval_list')
def main(bed_ref, sd_ref, out_ref):
    """
    Converts a BED file to a .interval_list file using Picard BedToIntervalList.
    The output file must have been added to references.py.
    """
    out_ref = reference_path(out_ref)
    if not out_ref.endswith('.interval_list'):
        raise ValueError('Output reference file must have a .interval_list extension.')
    get_ref_files(bed_ref, sd_ref, out_ref)


if __name__ == '__main__':
    main()
