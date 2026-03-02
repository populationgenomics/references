"""
IMPORTANT: Read the README document before using this script.
"""

import click
import re

from cpg_utils import to_path
from cpg_utils.hail_batch import reference_path


def read_fasta_dict(fasta_dict_path):
    """Read the fasta.dict file and return a dictionary mapping sequence names to their hashes."""
    fasta_dict = {}
    with open(to_path(fasta_dict_path)) as f:
        for line in f:
            if not line.startswith('@SQ'):
                continue
            # Lines look like: "@SQ	SN:chr1	LN:XXXX	M5:<hash>	UR:file:/path/to/ref.fasta"
            seq_name, seq_hash = re.search(r'SN:(\S+)', line).group(1), re.search(r'M5:(\S+)', line).group(1)
            fasta_dict[seq_name] = seq_hash
    return fasta_dict


def rewrite_interval_list(interval_list_path, fasta_dict, output_path):
    with open(to_path(interval_list_path)) as f:
        lines = f.readlines()
    with to_path(output_path).open('w') as f:
        for line in lines[:]:
            if line.startswith('@SQ'):
                seq_name = re.search(r'SN:(\S+)', line).group(1)
                current_hash = re.search(r'M5:(\S+)', line).group(1)
                if seq_name in fasta_dict:
                    new_hash = fasta_dict[seq_name]
                    if new_hash != current_hash:
                        new_line = line.replace(f'M5:{current_hash}', f'M5:{new_hash}')
                        print(f'Updated hash for {seq_name}: {current_hash} -> {new_hash}')
                        f.write(new_line)
                        continue

            f.write(line)


def get_ref_files(interval_list_ref: str, fasta_ref: str, outfile_path: str) -> None:
    """
    Rewrites the sequence dictionary in an interval list file using the hashes from a fasta.dict file.
    The output file must have been added to references.py.
    """
    interval_list_file = reference_path(interval_list_ref)
    fasta_dict_file = reference_path(fasta_ref).with_suffix('.dict')
    
    fasta_dict = read_fasta_dict(fasta_dict_file)

    rewrite_interval_list(interval_list_file, fasta_dict, outfile_path)
    print('Done.')


@click.command()
@click.option(
    '--interval-list-ref',
    required=True,
    help='String identifier for the input interval_list file from the references.',
    default='broad/genome_coverage_interval_list',
)
@click.option(
    '--fasta-ref',
    required=True,
    help='String identifier for the fasta sequence file from the references.',
    default='broad/ref_fasta',
)
@click.option(
    '--out-ref',
    required=True,
    help='Reference path to save the output .interval_list file.',
    default='genome_coverage_interval_list_masked/interval_list',
)
def main(interval_list_ref, fasta_ref, out_ref):
    """
    Rewrites the sequence dictionary in an interval list file using the hashes from a fasta.dict file.
    The output file must have been added to references.py.
    """
    out_ref = str(reference_path(out_ref))
    if not out_ref.endswith('.interval_list'):
        raise ValueError('Output reference file must have a .interval_list extension.')
    get_ref_files(interval_list_ref, fasta_ref, out_ref)


if __name__ == '__main__':
    main()
