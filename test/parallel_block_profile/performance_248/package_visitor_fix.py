"""Reseal the bbserv source after removing executable-stack callbacks."""
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

from package_nullfix import sha, verify

ROOT = Path(__file__).resolve().parents[3]
BASE_SHA = '7469cea885059eb38ccb075a67b9751fd11de14b8e21003aee07cbfb72729955'
GDB_SHA = '48ffe334861d7a6cfc52f003c1f750bbeecbc01f2bd4ee304cbeffdcf12d4b4a'
CHANGED = ('src/parallel_block_primitive_walk.f90',
           'src/parallel_block_inverse.f90',
           'src/parallel_block_primitives.f90')


def package(out):
    base = ROOT/'.stage248-validation/bbserv-candidate3'
    archive_base = base.with_suffix('.tar.gz')
    out = Path(out).resolve()
    archive = out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Fresh package destination required')
    if sha(archive_base) != BASE_SHA:
        raise ValueError('Corrected v4 source archive differs')
    previous = verify(base)
    shutil.copytree(base, out, symlinks=True)
    changes = {}
    for name in CHANGED:
        fixed = ROOT/name
        member = 'source/'+name
        if sha(fixed) == previous[member]:
            raise ValueError('Visitor correction is absent: '+name)
        shutil.copy2(fixed, out/member)
        changes[name] = dict(before=previous[member], after=sha(fixed))
    provenance = dict(
        status='prepared-for-bbserv-validation',
        base_archive_sha256=BASE_SHA,
        gdb_job=55199,
        gdb_review_sha256=GDB_SHA,
        finding='All 16 ranks attempted to execute an internal-procedure visitor '
                'trampoline from a non-executable stack at primitive walk emit.',
        correction='Collect primitive visits and invoke the internal handlers '
                   'directly; use a module procedure for pressure registration.',
        changed_sources=changes,
        local_check='Serial J4 feature-off baseline and candidate builds plus '
                    'local four-rank exact trajectory checks passed.',
        cluster_check='Pending nk/bb13 rerun',
    )
    (out/'VISITOR_CORRECTION.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'README.md').write_text(
        '# RK candidate for bbserv after callback correction\n\n'
        'This kit builds on bbserv-candidate3. VISITOR_CORRECTION.json records '
        'the GDB finding and the three source changes that remove internal '
        'procedure trampolines from primitive planning while retaining the '
        'non-executable stack. CORRECTION.json records the earlier optional '
        'argument change. LOCAL_VALIDATION.json is the original Stage 247 '
        'qualification record; this cumulative source still requires bb13 '
        'validation. Use the rk248-performance-kit-v5 launcher.\n')
    files = {str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*'))
             if p.is_file() and p.name != 'files-sha256.json'}
    (out/'files-sha256.json').write_text(json.dumps(files, indent=2)+'\n')
    with tarfile.open(archive, 'w:gz') as stream:
        stream.add(out, arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk248-visitors-verify-') as temp:
        with tarfile.open(archive) as stream:
            stream.extractall(temp, filter='data')
        if verify(Path(temp)/out.name) != files:
            raise ValueError('Unpacked source manifest differs')
    return dict(archive=str(archive), sha256=sha(archive), files=len(files),
                changes=changes, unpack_verified=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        raise SystemExit('Usage: package_visitor_fix.py FRESH_SOURCE_KIT_DIR')
    print(json.dumps(package(sys.argv[1]), indent=2))
