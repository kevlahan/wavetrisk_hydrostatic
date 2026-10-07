"""Build a self-contained publication candidate from the sealed Stage264 source anchor."""
import argparse
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
from common import save, sha, verify, source_scope

HERE = Path(__file__).resolve().parent
CHANGED = ('src/parallel_block_inverse.f90', 'src/parallel_block_mpi.f90')


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*'))
            if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts and p.name != 'files-sha256.json'}


def prepare(repo, anchor, prior, destination):
    repo, anchor, prior = [Path(p).resolve(strict=True) for p in (repo, anchor, prior)]
    verify(prior)
    prior_identity = json.loads((prior / 'IDENTITY.json').read_text())
    if prior_identity['stage'] != 265 or prior_identity.get('revision'):
        raise ValueError('Require the original sealed Stage265 candidate')
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError('Fresh kit destination required')
    verify(anchor)
    identity = json.loads((anchor / 'IDENTITY.json').read_text())
    if identity['stage'] != 264:
        raise ValueError('Require the verified Stage264 reset kit')
    destination.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk265t-source-') as tmp:
        root = Path(tmp)
        with tarfile.open(anchor / 'sources.tar.gz') as stream:
            stream.extractall(root, filter='data')
        original = json.loads((anchor / 'source-files.json').read_text())
        for kind in ('stage230', 'main'):
            for name, expected in original[kind].items():
                if sha(root / kind / name) != expected:
                    raise ValueError('Anchor source differs: ' + kind + '/' + name)
        shutil.copytree(root / 'stage230', root / 'candidate', symlinks=True)
        for name in CHANGED:
            shutil.copy2(repo / name, root / 'candidate' / name)
        manifests = {kind: inventory(root / kind) for kind in ('stage230', 'main', 'candidate')}
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
    identity.update(stage=265, revision='tiled', prior_candidate_source_sha256=prior_identity['candidate_source_sha256'],
                    prior_candidate_manifest_sha256=sha(prior / 'files-sha256.json'), changed_source=list(CHANGED), anchor_manifest_sha256=sha(anchor / 'files-sha256.json'),
                    candidate_source_sha256={name: sha(repo / name) for name in CHANGED},
                    scope='Patch-tiled refinement of direct canonical inverse SOL publication; both installations and preflight retained.')
    save(destination / 'IDENTITY.json', identity)
    source_scope(destination)
    save(destination / 'files-sha256.json', inventory(destination))
    verify(destination)
    print('Prepared: ' + str(destination))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'anchor', 'prior', 'destination'):
        parser.add_argument('--' + name, type=Path, required=True)
    prepare(**vars(parser.parse_args()))
