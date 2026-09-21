#!/usr/bin/env python3
"""Apply the Stage223 placement repair without touching source/build identities."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile


def digest(path):
    if path.is_symlink():
        return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def apply(repair, kit):
    old = json.loads((repair/'original-helpers-sha256.json').read_text())
    new = json.loads((repair/'helpers-sha256.json').read_text())
    current = json.loads((kit/'helpers-sha256.json').read_text())
    if current not in (old, new):
        raise ValueError('Unexpected installed helper manifest; no files changed')
    changed = [name for name in new if old.get(name) != new[name]]
    # Permit resuming a partially copied repair, but reject unrelated edits.
    for name, expected in new.items():
        path = kit/name
        if name not in old and not path.exists():
            continue
        if digest(path) not in {old.get(name), expected}:
            raise ValueError('Unexpected installed file: '+name)
    for name in changed:
        if digest(repair/name) != new[name]:
            raise ValueError('Damaged repair payload: '+name)
    if current == new and all(digest(kit/n) == new[n] for n in changed):
        print('Placement repair already installed.')
        return
    backup = Path(tempfile.mkdtemp(prefix='stage223-launcher-backup-', dir=kit.parent))
    for name in [*changed, 'helpers-sha256.json']:
        if (kit/name).exists():
            shutil.copy2(kit/name, backup/name)
    for name in [*changed, 'helpers-sha256.json']:
        temporary = kit/(name+'.repair-tmp')
        shutil.copy2(repair/name, temporary)
        temporary.replace(kit/name)
    for name, expected in new.items():
        if digest(kit/name) != expected:
            raise ValueError('Repaired file verification failed: '+name)
    print('Placement repair verified. Original helpers saved in '+str(backup))


if __name__ == '__main__':
    apply(Path(__file__).resolve().parent, Path(sys.argv[1]).resolve(strict=True))
