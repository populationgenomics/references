"""
Prepare test matrix (to transfer references in parallel)
"""

import os
import sys
from os.path import join
from cloudpathlib import AnyPath
from references import SOURCES as NEW_SOURCES

try:
    from references_before import SOURCES as OLD_SOURCES
except ImportError:
    OLD_SOURCES = []

REFERENCES_PREFIX = os.environ['REFERENCES_PREFIX']

transfers = {}
for source in NEW_SOURCES:
    old_sources_d = {s.name: s for s in OLD_SOURCES}
    dst_path = join(REFERENCES_PREFIX, source.dst)
    is_changed = (
        source.name not in old_sources_d
        or source.src != old_sources_d[source.name].src
        or source.dst != old_sources_d[source.name].dst
        or not AnyPath(dst_path).exists()
    )
    if is_changed and source.src and source.transfer_cmd:
        print(f'{source.name} has changed, will transfer', file=sys.stderr)
        transfers[source.name] = {'src': source.src, 'dst': dst_path}
    else:
        print(f'{source.name} has not changed since previous revision', file=sys.stderr)

if transfers:
    matrix = {
        'include': [
            {
                'name': name,
                'src': data['src'],
                'dst': data['dst'],
            }
            for name, data in transfers.items()
        ]
    }
else:
    matrix = {}
print(str(matrix).replace(' ', ''), file=sys.stderr)
print(' '.join(str(matrix).replace(' ', '')), file=sys.stderr)
print(str(matrix).replace(' ', ''), end='')
