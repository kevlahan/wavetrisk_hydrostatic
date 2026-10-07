"""Build a source-bound installation candidate from the accepted sealed Stage265t kit."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import BASELINE, save, sha, verify, source_scope

HERE = Path(__file__).resolve().parent
CHANGED = ('src/parallel_block_inverse.f90', 'src/parallel_block_mpi.f90')


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*'))
            if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts and p.name != 'files-sha256.json'}


def prepare(repo, anchor, destination):
    repo, anchor = [Path(p).resolve(strict=True) for p in (repo, anchor)]
    verify(anchor)
    identity = json.loads((anchor / 'IDENTITY.json').read_text())
    if identity['stage'] != 265 or identity.get('revision') != 'tiled':
        raise ValueError('Require the accepted sealed Stage265 tiled-publication kit')
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError('Fresh kit destination required')
    destination.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk266-source-') as tmp:
        root = Path(tmp)
        with tarfile.open(anchor / 'sources.tar.gz') as stream:
            stream.extractall(root, filter='data')
        original = json.loads((anchor / 'source-files.json').read_text())
        for name in identity['committed_inputs_checked']:
            committed = subprocess.check_output(['git', 'show', BASELINE + ':' + name], cwd=repo)
            if __import__('hashlib').sha256(committed).hexdigest() != original['candidate'][name]:
                raise ValueError('Accepted baseline differs from committed source: ' + name)
        for kind in ('candidate', 'main'):
            for name, expected in original[kind].items():
                if sha(root / kind / name) != expected:
                    raise ValueError('Anchor source differs: ' + kind + '/' + name)
        shutil.rmtree(root / 'stage230')
        (root / 'candidate').rename(root / 'stage265t')
        shutil.copytree(root / 'stage265t', root / 'candidate', symlinks=True)
        for name in CHANGED:
            shutil.copy2(repo / name, root / 'candidate' / name)
        manifests = {kind: inventory(root / kind) for kind in ('stage265t', 'main', 'candidate')}
        save(destination / 'source-files.json', manifests)
        with tarfile.open(destination / 'sources.tar.gz', 'w:gz', format=tarfile.PAX_FORMAT) as stream:
            for kind in manifests:
                stream.add(root / kind, arcname=kind)
    for path in anchor.iterdir():
        if path.is_file() and path.suffix in ('.py', '.sh'):
            shutil.copy2(path, destination / path.name)
    for path in HERE.iterdir():
        if path.is_file() and path.suffix in ('.py', '.sh', '.sbatch', '.md'):
            shutil.copy2(path, destination / path.name)
    for name in ('grids', 'fixture-J4'):
        shutil.copytree(anchor / name, destination / name)
    identity.update(stage=266, revision='deferred-installation', stage265t_commit=BASELINE,
                    baseline_source_sha256={name: manifests['stage265t'][name] for name in CHANGED},
                    changed_source=list(CHANGED), anchor_manifest_sha256=sha(anchor / 'files-sha256.json'),
                    candidate_source_sha256={name: sha(repo / name) for name in CHANGED},
                    scope='Omit first full installation only for provisional Stage203 production; retain final boundary installation, accepted/oracle/fallback first installation and all publication checks.')
    identity.pop('prior_candidate_source_sha256', None)
    identity.pop('prior_candidate_manifest_sha256', None)
    save(destination / 'IDENTITY.json', identity)
    source_scope(destination)
    save(destination / 'files-sha256.json', inventory(destination))
    verify(destination)
    print('Prepared: ' + str(destination))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'anchor', 'destination'):
        parser.add_argument('--' + name, type=Path, required=True)
    prepare(**vars(parser.parse_args()))
