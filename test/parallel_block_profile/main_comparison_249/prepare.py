"""Export a read-only main snapshot with the established checkpoint I/O fix."""
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inventory(root):
    return {str(path.relative_to(root)): sha(path) for path in sorted(root.rglob('*'))
            if path.is_file() and path.name != 'files-sha256.json'}


def prepare(repo, out):
    repo, out = Path(repo).resolve(), Path(out).resolve()
    archive = out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Use a fresh main snapshot destination')
    main = subprocess.check_output(['git', 'rev-parse', 'main'], cwd=repo, text=True).strip()
    patch = repo/'test/parallel_block_profile/legacy_checkpoint_io_190.patch'
    out.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk249-main-') as temporary:
        export = Path(temporary)/'main.tar'
        subprocess.run(['git', 'archive', '--format=tar', '--output='+str(export), main],
                       cwd=repo, check=True)
        with tarfile.open(export) as stream:
            stream.extractall(out, filter='data')
    before = inventory(out)
    subprocess.run(['git', 'apply', '--check', str(patch)], cwd=out, check=True)
    subprocess.run(['git', 'apply', str(patch)], cwd=out, check=True)
    after = inventory(out)
    changed = sorted(name for name in set(before)|set(after) if before.get(name) != after.get(name))
    if changed != ['src/checkpoint.f90']:
        raise ValueError('Unexpected change to read-only main export: '+str(changed))
    identity = dict(main_commit=main, checkpoint_io_patch_sha256=sha(patch),
                    changed_source=changed,
                    scope='Read-only main snapshot with established Stage190 checkpoint I/O hardening. '
                          'External Simple-physics objects are shared with the candidate at build time.')
    (out/'MAIN_IDENTITY.json').write_text(json.dumps(identity, indent=2)+'\n')
    files = inventory(out)
    (out/'files-sha256.json').write_text(json.dumps(files, indent=2)+'\n')
    with tarfile.open(archive, 'w:gz') as stream:
        stream.add(out, arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk249-main-check-') as temporary:
        with tarfile.open(archive) as stream:
            stream.extractall(temporary, filter='data')
        extracted = Path(temporary)/out.name
        if inventory(extracted) != files:
            raise ValueError('Main archive unpack verification failed')
    return dict(archive=str(archive), sha256=sha(archive), files=len(files),
                identity=identity, unpack_verified=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        raise SystemExit('Usage: prepare.py REPO FRESH_MAIN_DIR')
    print(json.dumps(prepare(*sys.argv[1:]), indent=2))
