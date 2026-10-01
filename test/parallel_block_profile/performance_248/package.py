"""Seal performance scripts separately from the already qualified source archive."""
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def package(source, out, candidate_archive):
    source,out,candidate_archive=(Path(p).resolve() for p in (source,out,candidate_archive))
    if out.exists() or out.with_suffix('.tar.gz').exists():
        raise ValueError('Fresh package destination required')
    if sha(candidate_archive)!='dc6daba22953f5ab54398f6085b9b39aa04a2161b233fade692fedaa9d5d7d63':
        raise ValueError('Qualified source archive differs')
    out.mkdir(parents=True)
    for path in sorted(source.iterdir()):
        if path.is_file() and path.name!='files-sha256.json':
            shutil.copy2(path,out/path.name)
    files={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'files-sha256.json').write_text(json.dumps(files,indent=2)+'\n')
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as tar:
        tar.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk248-verify-') as directory:
        with tarfile.open(archive) as tar:
            tar.extractall(directory,filter='data')
        unpacked=Path(directory)/out.name
        for name,expected in files.items():
            if sha(unpacked/name)!=expected:raise ValueError('Unpacked script differs: '+name)
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),script_files=len(files),
                          source_archive_sha256=sha(candidate_archive),unpack_verified=True),indent=2))


if __name__=='__main__':
    import sys
    if len(sys.argv)!=4:raise SystemExit('Usage: package.py SOURCE_DIR OUT_DIR QUALIFIED_SOURCE_ARCHIVE')
    package(*sys.argv[1:])
