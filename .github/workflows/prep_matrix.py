import os
from os.path import join
from references import SOURCES as NEW_SOURCES
try:
    from references_before import SOURCES as OLD_SOURCES
except ImportError:
    OLD_SOURCES = []

REFERENCES_PREFIX = os.environ['REFERENCES_PREFIX']

transfers = {}
for source in NEW_SOURCES:
    old_sources_d = {s.name: s for s in OLD_SOURCES}
    is_changed = (
        source.name not in old_sources_d 
        or source != old_sources_d[source.name]
    )
    dst_path = join(REFERENCES_PREFIX, source.dst)
    if is_changed and source.src:
        transfers[source.name] = {'src': source.src, 'dst': dst_path}

d = {'include': [{
      'name': name, 
      'src': data['src'],
      'dst': data['dst'],
    } for name, data in transfers.items()
]}
print(str(d).replace(' ', ''), end='')
