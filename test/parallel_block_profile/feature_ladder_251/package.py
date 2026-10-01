"""Seal the early-feature ladder scripts against both source archives."""
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

FILES=('README.md','build.py','collect.py','collect.sh','launch.sh','performance.sbatch',
       'placement.py','preflight.py','run.py','timing.py')
CANDIDATE_SHA='7189c938950f277378c30b8f546d7154f818e5262cc9acebfc79d7260a9cbaac'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def package(source,out,candidate_archive,main_archive):
    source,out,candidate_archive,main_archive=map(Path,(source,out,candidate_archive,main_archive))
    archive=out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Fresh performance kit destination required')
    if sha(candidate_archive)!=CANDIDATE_SHA:
        raise ValueError('Qualified candidate source archive differs')
    with tarfile.open(main_archive) as stream:
        member=stream.extractfile('rk249-main-legacy/MAIN_IDENTITY.json')
        if member is None:
            raise ValueError('Main identity absent')
        identity=json.load(member)
    out.mkdir(parents=True)
    for name in FILES:
        shutil.copy2(source/name,out/name)
    archives=dict(candidate_archive_sha256=CANDIDATE_SHA,main_archive_sha256=sha(main_archive),
                  main_identity=identity)
    (out/'ARCHIVES.json').write_text(json.dumps(archives,indent=2)+'\n')
    files={path.name:sha(path) for path in sorted(out.iterdir()) if path.is_file()}
    (out/'files-sha256.json').write_text(json.dumps(files,indent=2)+'\n')
    with tarfile.open(archive,'w:gz') as stream:
        stream.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk251-unpack-') as temporary:
        with tarfile.open(archive) as stream:
            stream.extractall(temporary,filter='data')
        unpacked=Path(temporary)/out.name
        if {path.name:sha(path) for path in unpacked.iterdir() if path.is_file() and path.name!='files-sha256.json'}!=files:
            raise ValueError('Performance kit unpack verification failed')
    return dict(archive=str(archive),sha256=sha(archive),script_files=len(files),
                candidate_archive_sha256=CANDIDATE_SHA,main_archive_sha256=sha(main_archive),
                main_commit=identity['main_commit'],unpack_verified=True)


if __name__=='__main__':
    import sys
    if len(sys.argv)!=5:
        raise SystemExit('Usage: package.py SOURCE_DIR FRESH_KIT_DIR CANDIDATE_ARCHIVE MAIN_ARCHIVE')
    print(json.dumps(package(*sys.argv[1:]),indent=2))
