"""Reject missing or changed archives, patch and J5 inputs before submission."""
import hashlib
import json
from pathlib import Path
import sys

from run import J5_HASHES


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def preflight(perf,candidate,main,fixture):
    perf,candidate,main,fixture = map(Path,(perf,candidate,main,fixture))
    for name,expected in json.loads((perf/'files-sha256.json').read_text()).items():
        if sha(perf/name)!=expected:
            raise ValueError('Performance kit file differs: '+name)
    patch=json.loads((perf/'PATCH.json').read_text())
    if sha(perf/'quiet_storage.patch')!=patch['patch_sha256']:
        raise ValueError('Storage diagnostic patch differs')
    archives=json.loads((perf/'ARCHIVES.json').read_text())
    if sha(candidate)!=archives['candidate_archive_sha256']:
        raise ValueError('Qualified candidate archive differs')
    if sha(main)!=archives['main_archive_sha256']:
        raise ValueError('Read-only main archive differs')
    for name,expected in (('simple.in',J5_HASHES['simple']),
                          ('test_checkpoint_0003.bin.zst',J5_HASHES['seed'])):
        if sha(fixture/name)!=expected:
            raise ValueError('Established J5 fixture differs: '+name)
    print('PASS: direct-main kit, both source archives and J5 fixture seals',flush=True)


if __name__=='__main__':
    if len(sys.argv)!=5:
        raise SystemExit('Usage: preflight.py PERF_KIT CANDIDATE_ARCHIVE MAIN_ARCHIVE J5_FIXTURE')
    preflight(*sys.argv[1:])
