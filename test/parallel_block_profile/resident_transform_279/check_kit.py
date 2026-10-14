"""Read-only seal/source verification, also usable before submission."""
import json,sys,tarfile,tempfile,importlib
from pathlib import Path
from common import require_identity,verify,source_scope,verify_packaged_grid,sha

def check(root):
    root=Path(root).resolve(strict=True)
    manifest=verify(root);delta=source_scope(root);verify_packaged_grid(root)
    require_identity(json.loads((root/'IDENTITY.json').read_text()))
    gate=json.loads((root/'ENTRY_GATE.json').read_text())
    tests=json.loads((root/'TEST_VALIDATION.json').read_text())
    if not gate.get('sealing_ready') or tests.get('status')!='passed':raise ValueError('Qualification incomplete')
    for name,h in tests['tested_helper_sha256'].items():
        if sha(root/name)!=h:raise ValueError('Helper differs from tests: '+name)
    for module in ('build','run','review','diagnostic','attribution','control','contract','capability_223','placement_278'):
        importlib.import_module(module)
    families=json.loads((root/'source-files.json').read_text())
    with tempfile.TemporaryDirectory(prefix='rk279-source-check-') as tmp:
        with tarfile.open(root/'sources.tar.gz') as t:t.extractall(tmp,filter='data')
        for family,files in families.items():
            for name,h in files.items():
                if sha(Path(tmp)/family/name)!=h:raise ValueError('Extracted source differs: '+family+'/'+name)
    return dict(status='passed',sealed_files=len(manifest),source_families=list(families),qualified_overlay_files=len(delta),source_files={k:len(v) for k,v in families.items()},submitted=False)
if __name__=='__main__':print(json.dumps(check(Path(__file__).resolve().parent),indent=2))
