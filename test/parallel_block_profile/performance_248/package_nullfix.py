"""Reseal the Stage 247 source with the bb13 boundary-import correction."""
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[3]
ORIGINAL_SHA = 'dc6daba22953f5ab54398f6085b9b39aa04a2161b233fade692fedaa9d5d7d63'
SANITIZER_SHA = '4ce3e09aebfbf82639ea0b7a1b5c5b1fef3a886991ea1d5eeb165b4074134c55'
SOURCE = 'source/src/parallel_block_mpi.f90'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify(directory):
    files = json.loads((directory/'files-sha256.json').read_text())
    for name, expected in files.items():
        if sha(directory/name) != expected:
            raise ValueError('Source kit file differs: '+name)
    return files


def package(out):
    original = ROOT/'.stage247-validation/bbserv-candidate2'
    original_archive = original.with_suffix('.tar.gz')
    fixed_source = ROOT/'src/parallel_block_mpi.f90'
    out = Path(out).resolve()
    archive = out.with_suffix('.tar.gz')
    if out.exists() or archive.exists():
        raise ValueError('Fresh package destination required')
    if sha(original_archive) != ORIGINAL_SHA:
        raise ValueError('Original Stage 247 archive differs')
    prior = verify(original)
    if sha(fixed_source) == prior[SOURCE]:
        raise ValueError('Correction is absent')
    shutil.copytree(original, out, symlinks=True)
    shutil.copy2(fixed_source, out/SOURCE)
    provenance = dict(
        status='prepared-for-bbserv-validation',
        original_archive_sha256=ORIGINAL_SHA,
        original_local_validation_sha256=prior['LOCAL_VALIDATION.json'],
        sanitizer_archive_sha256=SANITIZER_SHA,
        sanitizer_job=55186,
        finding='Absent optional domain_sol forwarded at two boundary packing calls; '
                'GCC 13 UBSan reported a null Float_Field load at startup.',
        corrected_file=SOURCE,
        original_file_sha256=prior[SOURCE],
        corrected_file_sha256=sha(fixed_source),
        local_check='Serial feature-off baseline and candidate builds passed; '
                    'four-rank J4 feature-off baseline one-step run exited 0.',
        cluster_check='Pending nk/bb13 rerun',
    )
    (out/'CORRECTION.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'README.md').write_text(
        '# Corrected RK candidate for bbserv\n\n'
        'This kit changes only source/src/parallel_block_mpi.f90 from the sealed '
        'Stage 247 candidate. It omits an absent optional domain_sol argument in '
        'two boundary packing calls. CORRECTION.json records the sanitizer finding '
        'and local checks. LOCAL_VALIDATION.json is the original Stage 247 '
        'qualification record; the corrected source still requires bb13 validation.\n\n'
        'Use the accompanying rk248-performance-kit-v4 launcher. It verifies '
        'this kit, builds matched binaries serially, and runs the J4/J5 checks '
        'and timing on partition nk (bb13).\n')
    files = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'files-sha256.json'}
    (out/'files-sha256.json').write_text(json.dumps(files, indent=2)+'\n')
    with tarfile.open(archive, 'w:gz') as stream:
        stream.add(out, arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk248-nullfix-verify-') as temp:
        with tarfile.open(archive) as stream:
            stream.extractall(temp, filter='data')
        unpacked = Path(temp)/out.name
        if verify(unpacked) != files:
            raise ValueError('Unpacked source manifest differs')
    return dict(archive=str(archive), sha256=sha(archive), files=len(files),
                corrected_source_sha256=sha(fixed_source), unpack_verified=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        raise SystemExit('Usage: package_nullfix.py FRESH_SOURCE_KIT_DIR')
    print(json.dumps(package(sys.argv[1]), indent=2))
