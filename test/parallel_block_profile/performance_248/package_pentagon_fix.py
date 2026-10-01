"""Reseal the bbserv source after removing remaining migrated stack callbacks."""
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

from package_nullfix import sha, verify

ROOT = Path(__file__).resolve().parents[3]
BASE_SHA = 'edc647fb8aff057444d0d4d9e93e2756eba1749c5fdcacff332e0404735b8eb0'
GDB_SHA = '509aa5f975a1bee2c36f307231a7fd56442ad0c0cab0622d91b1e30670649919'
CHANGED = ('src/domain_ops.f90',
           'src/parallel_block_pentagons.f90',
           'src/parallel_block_inverse.f90',
           'src/parallel_block_mpi.f90')


def package(out):
    base = ROOT/'.stage248-validation/bbserv-candidate4'
    archive_base = base.with_suffix('.tar.gz')
    out = Path(out).resolve()
    archive = out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Fresh package destination required')
    if sha(archive_base) != BASE_SHA:
        raise ValueError('V5 source archive differs')
    previous = verify(base)
    shutil.copytree(base, out, symlinks=True)
    changes = {}
    for name in CHANGED:
        fixed = ROOT/name
        member = 'source/'+name
        if sha(fixed) == previous[member]:
            raise ValueError('Callback correction is absent: '+name)
        shutil.copy2(fixed, out/member)
        changes[name] = dict(before=previous[member], after=sha(fixed))
    provenance = dict(
        status='prepared-for-bbserv-validation',
        base_archive_sha256=BASE_SHA,
        failed_performance_job=55200,
        gdb_job=55233,
        gdb_review_sha256=GDB_SHA,
        finding='All 16 ranks attempted to execute an internal-procedure '
                'callback from a non-executable stack in apply_to_penta_d '
                'during pentagon plan construction.',
        correction='Return ordered pentagon visits as data and invoke the '
                   'local handler directly. Lift the remaining retained-RK '
                   'callbacks to module procedures.',
        changed_sources=changes,
        local_check='Serial J4 feature-off baseline and candidate builds, '
                    'a baseline first step with 31 exact VTK members, and '
                    'a four-rank exact candidate trajectory passed.',
        cluster_check='Pending nk/bb13 rerun',
    )
    (out/'PENTAGON_CORRECTION.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'README.md').write_text(
        '# RK candidate for bbserv after pentagon callback correction\n\n'
        'This kit builds on bbserv-candidate4. PENTAGON_CORRECTION.json '
        'records the GDB finding and the additional source changes. '
        'VISITOR_CORRECTION.json and CORRECTION.json record the earlier '
        'callback and optional-argument changes. LOCAL_VALIDATION.json is '
        'the original Stage 247 qualification record; this cumulative '
        'source still requires bb13 validation. Use the '
        'rk248-performance-kit-v6 launcher.\n')
    files = {str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*'))
             if p.is_file() and p.name != 'files-sha256.json'}
    (out/'files-sha256.json').write_text(json.dumps(files, indent=2)+'\n')
    with tarfile.open(archive, 'w:gz') as stream:
        stream.add(out, arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk248-pentagons-verify-') as temp:
        with tarfile.open(archive) as stream:
            stream.extractall(temp, filter='data')
        if verify(Path(temp)/out.name) != files:
            raise ValueError('Unpacked source manifest differs')
    return dict(archive=str(archive), sha256=sha(archive), files=len(files),
                changes=changes, unpack_verified=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        raise SystemExit('Usage: package_pentagon_fix.py FRESH_SOURCE_KIT_DIR')
    print(json.dumps(package(sys.argv[1]), indent=2))
