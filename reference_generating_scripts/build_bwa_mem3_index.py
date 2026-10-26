from argparse import ArgumentParser

from cpg_utils.hail_batch import get_batch


DEFAULT_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bwa_mem3_ora_dupblaster:v0.15.0-3'


def main(fasta: str, outpath: str, image: str):
    """
    Read a FASTA file and build bwa mem3 index.
    Write a new collection of files - the original Fasta, and all adjacent index files.
    """

    assert outpath.endswith('.fasta'), 'this process requires the output to end with .fasta'

    batch = get_batch('Generate bwa-mem3 fasta index files')

    job = batch.new_bash_job(name='Generate bwa-mem3 fasta index files')

    job.image(image=image)

    # chunky resourcing
    job.memory('100GiB')
    job.storage('100GiB')
    job.cpu(16)

    # index elements to be extracted
    # contains fasta as well - the intention is to write the whole collection of files, not just the new elements
    # this should allow for a more intentional output path (e.g. gs://references/bwamem3/masked/FILE.fasta & companions)
    job.declare_resource_group(
        output={
            'fasta': '{root}.fasta',
            'fasta.bwt.2bit.64': '{root}.fasta.bwt.2bit.64',
            'fasta.amb': '{root}.fasta.amb',
            'fasta.ann': '{root}.fasta.ann',
            'fasta.pac': '{root}.fasta.pac',
        }
    )

    # read in the fasta file
    local_fa = batch.read_input(fasta)

    # move the input fasta file to an extractable output location
    job.command(f'mv {local_fa} {job.output.fasta}')

    # generate the new index - memory and threads left as default to auto-assign
    job.command(f'bwa-mem3 index --tmp-dir $BATCH_TMPDIR {job.output.fasta}')

    # extract the output
    batch.write_output(job.output, outpath.removesuffix('.fasta'))

    batch.run(wait=False)


if __name__ == '__main__':
    parser = ArgumentParser(description='Build bwa-mem3 index')
    parser.add_argument('-i', help='fasta input path, in GCS')
    parser.add_argument('-o', help='output fasta to use as a base, in GCS. Must end with ".fasta"')
    parser.add_argument('-b', help='bwa image', default=DEFAULT_IMAGE)
    args = parser.parse_args()

    main(fasta=args.i, outpath=args.o, image=args.b)
