"""Package only after the real four-rank correctness replay has passed."""
import json
from pathlib import Path
import shutil
import sys
import tarfile
from common import save, sha, verify


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*'))
            if (p.is_file() and '__pycache__' not in p.parts and p.name != 'files-sha256.json')}


def seal(kit, validation):
    kit, validation = Path(kit).resolve(strict=True), Path(validation).resolve(strict=True)
    result = json.loads(validation.read_text())
    if result['status'] != 'passed' or result['ranks'] != 4 or set(result['cases']) != {'stage230', 'checked', 'stage230-oracle', 'main'}:
        raise ValueError('Require all four real local correctness cases')
    if sha(kit / 'source-files.json') != result['source_manifest_sha256']:
        raise ValueError('Validated and packaged sources differ')
    build_notes = json.loads((validation.parent / 'builds.json').read_text())
    if sha(validation.parent / 'builds.json') != result['build_manifest_sha256'] or build_notes['status'] != 'passed':
        raise ValueError('Validated builds differ')
    here = Path(__file__).resolve().parent
    for path in here.iterdir():
        if path.is_file() and path.suffix in ('.py', '.sh', '.sbatch', '.md'):
            shutil.copy2(path, kit / path.name)
    save(kit / 'LOCAL_VALIDATION.json', dict(result, build_notes=build_notes,
         validation_sha256=sha(validation), orchestration_tests='Four workflow tests, full import graph and shell syntax passed.'))
    save(kit / 'files-sha256.json', inventory(kit))
    verify(kit)
    archive = kit.with_suffix('.tar.gz')
    if archive.exists():
        raise ValueError('Use a fresh archive path')
    temporary = archive.with_suffix('.tmp')
    def keep(info):
        return None if '__pycache__' in Path(info.name).parts else info
    with tarfile.open(temporary, 'w:gz', format=tarfile.PAX_FORMAT) as stream:
        stream.add(kit, arcname=kit.name, filter=keep)
    temporary.rename(archive)
    return dict(archive=str(archive), sha256=sha(archive), bytes=archive.stat().st_size,
                files=len(inventory(kit)), source_manifest_sha256=sha(kit / 'source-files.json'),
                local_validation_sha256=sha(kit / 'LOCAL_VALIDATION.json'))


if __name__ == '__main__':
    print(json.dumps(seal(*sys.argv[1:]), indent=2))
