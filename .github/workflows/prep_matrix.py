"""
Prepare test matrix (to transfer references in parallel)
"""

import argparse
import sys
from os.path import join

from google.cloud import storage

from references import SOURCES as NEW_SOURCES

try:
    # copied into place by the github action
    from references_before import SOURCES as OLD_SOURCES
except ImportError:
    OLD_SOURCES = []

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
    transfers = {}
    for source in NEW_SOURCES:
        old_sources_d = {s.name: s for s in OLD_SOURCES}
        dst_path = join(references_prefix, source.dst)

        if source.src and source.transfer_cmd:
            if not gcs_file_exists(dst_path):
                print(f'{dst_path} does not exist, will transfer', file=sys.stderr)
                transfers[source.name] = {'src': source.src, 'dst': dst_path}
                continue
            elif (
                source.name not in old_sources_d
                or source.src != old_sources_d[source.name].src
                or source.dst != old_sources_d[source.name].dst
            ):
                print(f'{source.name} has changed, will transfer', file=sys.stderr)
                transfers[source.name] = {'src': source.src, 'dst': dst_path}
            else:
                print(
                    f'{source.name} has not changed since previous revision',
                    file=sys.stderr,
                )

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
