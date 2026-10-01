"""
Prepare test matrix (to transfer references in parallel)
"""

import argparse
import sys
from os.path import join

from google.cloud import storage

from references import SOURCES

GCS_CLIENT: storage.Client = storage.Client()


def gcs_file_exists(path: str) -> bool:
    """Check if file exists in GCS

    Args:
        path (str): A path to a file in GCS

    Returns:
        bool: True if the file exists, else False
    """
    assert path.startswith('gs://'), f'Invalid path: {path}, must start with gs://'
    bucket_name, blob_name = path.removeprefix('gs://').split('/', maxsplit=1)
    bucket = GCS_CLIENT.get_bucket(bucket_name)
    blob = bucket.get_blob(blob_name)
    if not blob:
        # fallback to see if it's a directory
        return gcs_directory_exists(path)
    return blob.exists()


def gcs_directory_exists(path: str) -> bool:
    """Check if directory exists in GCS

    Args:
        path (str): A path to a directory in GCS

    Returns:
        bool: True if the directory exists, else False
    """
    assert path.startswith('gs://'), f'Invalid path: {path}, must start with gs://'
    bucket_name, blob_name = path.removeprefix('gs://').split('/', maxsplit=1)
    bucket = GCS_CLIENT.get_bucket(bucket_name)
    if not blob_name.endswith('/'):
        # this is surprisingly important
        blob_name = blob_name + '/'

    query = bucket.list_blobs(prefix=blob_name, delimiter='/')
    if next(query, False):
        # has at least one blob
        return True
    return False


def generate_matrix(references_prefix: str) -> dict:
    """Generate matrix for transferring references in parallel

    Args:
        references_prefix (str): References prefix

    Returns:
        dict: {"include": [<list of transfers>]}
    """
    # A source is scheduled only for entries missing from the bucket; transfer.py
    # copies just those. Reference data is never overwritten in place: a different
    # upstream belongs under a new dst, the same way vep/105, 110 and 115 sit side by side.
    transfers = {}
    for source in SOURCES:
        if source.src and source.transfer_cmd:
            missing = [
                dst
                for _, _, dst in source.transfers(references_prefix)
                if not gcs_file_exists(dst)
            ]
            if missing:
                print(f'{missing} do not exist, will transfer', file=sys.stderr)
                transfers[source.name] = {
                    'src': source.src,
                    'dst': join(references_prefix, source.dst),
                }
            else:
                print(f'{source.name} is complete', file=sys.stderr)

    if not transfers:
        return {}
    return {
        'include': [
            {
                'name': name,
                'src': data['src'],
                'dst': data['dst'],
            }
            for name, data in transfers.items()
        ]
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references-prefix', help='Prefix for references')
    return parser.parse_args()


def print_matrix(matrix: dict):
    print(str(matrix).replace(' ', ''), end='', file=sys.stderr)
    print(str(matrix).replace(' ', ''), end='')


if __name__ == '__main__':
    args = parse_args()
    matrix = generate_matrix(references_prefix=args.references_prefix)
    print_matrix(matrix)
