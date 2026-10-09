"""Seal Stage266/main and the two-file fused primitive prototype."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import BASELINE, MAIN, save, sha, verify, source_scope
HERE=Path(__file__).resolve().parent

def inventory(root):
    return {str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*'))
            if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts and p.name!='files-sha256.json'}

def prepare(repo,anchor,destination):
    repo,anchor=map(lambda p:Path(p).resolve(strict=True),(repo,anchor))
    verify(anchor)
    identity=json.loads((anchor/'IDENTITY.json').read_text())
    if identity['stage']!=266 or identity['revision']!='deferred-installation' or identity['main_commit']!=MAIN:
        raise ValueError('Require sealed Stage266 campaign anchor')
    destination=Path(destination).resolve()
    if destination.exists():raise ValueError('Fresh kit directory required')
    destination.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk276-sources-') as tmp:
        root=Path(tmp)
        with tarfile.open(anchor/'sources.tar.gz') as stream:stream.extractall(root,filter='data')
        prior=json.loads((anchor/'source-files.json').read_text())
        for kind in ('candidate','main'):
            for name,expected in prior[kind].items():
                if sha(root/kind/name)!=expected:raise ValueError('Anchor source differs: '+kind+'/'+name)
        for name in identity['committed_inputs_checked']:
            actual=subprocess.check_output(['git','show',BASELINE+':'+name],cwd=repo)
            if hashlib.sha256(actual).hexdigest()!=prior['candidate'][name]:raise ValueError('Baseline commit differs: '+name)
        shutil.rmtree(root/'stage265t')
        (root/'candidate').rename(root/'stage266')
        shutil.copytree(root/'stage266',root/'candidate',symlinks=True)
        from instrument import instrument
        changed=instrument(root/'candidate')
        manifests={kind:inventory(root/kind) for kind in ('stage266','main','candidate')}
        save(destination/'source-files.json',manifests)
        with tarfile.open(destination/'sources.tar.gz','w:gz',format=tarfile.PAX_FORMAT) as stream:
            for kind in manifests:stream.add(root/kind,arcname=kind)
    for p in anchor.iterdir():
        if p.is_file() and p.suffix in ('.py','.sh') and p.name!='repeat.sh':shutil.copy2(p,destination/p.name)
    for p in HERE.iterdir():
        if p.is_file() and p.suffix in ('.py','.sh','.sbatch','.md','.f90','.json'):shutil.copy2(p,destination/p.name)
    for name in ('grids','fixture-J4'):shutil.copytree(anchor/name,destination/name)
    identity.update(stage=276,revision='primitive-prototype',analysis_schema_version=3,stage266_commit=BASELINE,changed_source=changed,
                    anchor_manifest_sha256=sha(anchor/'files-sha256.json'),
                    scope='Executable positive-interior fused primitive prototype with resident block halo checks. Reference packets and physics remain Domain-fed validation adapters; no production replacement or speedup claim.')
    identity.pop('stage265t_commit',None)
    identity.update(candidate_source_sha256={name:manifests['candidate'][name] for name in changed},
                    baseline_source_sha256={name:manifests['stage266'][name] for name in changed})
    save(destination/'IDENTITY.json',identity)
    source_scope(destination)
    save(destination/'files-sha256.json',inventory(destination));verify(destination)
    print('Prepared: '+str(destination))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','anchor','destination'):p.add_argument('--'+name,type=Path,required=True)
    prepare(**vars(p.parse_args()))
