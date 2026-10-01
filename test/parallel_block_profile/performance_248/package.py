"""Seal performance scripts separately from the corrected source archive."""
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
    if sha(candidate_archive)!='7189c938950f277378c30b8f546d7154f818e5262cc9acebfc79d7260a9cbaac':
        raise ValueError('Corrected source archive differs')
    out.mkdir(parents=True)
    for path in sorted(source.iterdir()):
        if path.is_file() and path.name not in ('files-sha256.json','test_scripts.py','package_nullfix.py','package_visitor_fix.py','package_pentagon_fix.py'):
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
