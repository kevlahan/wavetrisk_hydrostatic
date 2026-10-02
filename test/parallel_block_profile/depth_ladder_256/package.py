"""Seal and unpack-verify the Stage256 bb13 depth-index trial kit."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile


FILES = ('README.md', 'preflight.py', 'run.py', 'launch.sh',
         'performance.sbatch', 'collect.py', 'collect.sh', 'transform.py')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def package(source, out):
    source, out = map(Path, (source, out))
    archive = out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Use a fresh kit destination')
    out.mkdir(parents=True)
    for name in FILES:
        shutil.copy2(source/name, out/name)
    shutil.copy2(source/'optimized-parallel_block_program.f90',
                 out/'optimized-parallel_block_program.f90')
    files = {path.name: sha(path) for path in sorted(out.iterdir()) if path.is_file()}
    (out/'files-sha256.json').write_text(json.dumps(files, indent=2)+'\n')
    with tarfile.open(archive, 'w:gz') as stream:
        stream.add(out, arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk256-unpack-') as temporary:
        with tarfile.open(archive) as stream:
            stream.extractall(temporary, filter='data')
        unpacked = Path(temporary)/out.name
        actual = {path.name: sha(path) for path in unpacked.iterdir()
                  if path.is_file() and path.name != 'files-sha256.json'}
        if actual != files:
            raise ValueError('Kit archive unpack differs')
    return dict(archive=str(archive), sha256=sha(archive), sealed_files=len(files),
                unpack_verified=True)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: package.py SOURCE FRESH_KIT_DIR')
    print(json.dumps(package(*sys.argv[1:]), indent=2))
