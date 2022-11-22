import os
from os.path import join
from references import SOURCES as NEW_SOURCES, GENOME_BUILD
try:
    from references_before import SOURCES as OLD_SOURCES
except ImportError:
    OLD_SOURCES = []
    
PREFIX = os.environ['REFERENCES_PREFIX']

toml = {'genome_build': GENOME_BUILD}
transfer_commands = {}
for source in NEW_SOURCES:
    old_sources_d = {s.name: s for s in OLD_SOURCES}
    is_changed = (
        source.name not in old_sources_d 
        or source != old_sources_d[source.name]
    )
    dst_path = join(PREFIX, source.dst)
    if is_changed and source.src:
        type_ = None
        if source.src.startswith('gs://'):
            type_ = 'gcs'
        if source.src.startswith('https://'):
            type_ = 'https'
        if type_:
            transfer_commands[source.name] = {
                'src': source.src, 'dst': dst_path, 'type': type_
            }
    if not source.files:
        toml[source.name] = str(dst_path)
    else:
        toml[source.name] = {
            k: join(dst_path, suffix)
            for k, suffix in source.files.items()
        }
            
d = {"include": [{
      "name": name, 
      "src": data['src'],
      "dst": data['dst'],
      "type": data['type'],
    } for name, data in transfer_commands.items()
]}
print(str(d).replace(" ", ""), end='')
