#!/usr/bin/env python3
"""Seal the new verification-only package; preserve all historical packages."""
from pathlib import Path
import argparse
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile

ARTIFACT = 'rk325-native-qualification-kit-v5p1'
MAX_PAYLOAD = 16 * 1024**2
FROZEN_MANIFEST_SHA = 'ff430542d8ab81078d522a2468a1ce1baa22bac501f922e534c839c7c47ee1bb'


def require(ok, message):
    if not ok:
        raise ValueError('native.package325.' + message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def safe(name):
    path = Path(name)
    require(type(name) is str and not path.is_absolute() and str(path) == name and
            '..' not in path.parts, 'unsafe_member')
    return path


def regular_bytes(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'regular_source:' + str(path))
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    require((before.st_size, before.st_mtime_ns, before.st_ino) ==
            (after.st_size, after.st_mtime_ns, after.st_ino) and len(data) == before.st_size,
            'source_changed_during_seal:' + str(path))
    return data


def frozen_payload(shared, source=None):
    """Return exact historical payload; never rewrite its source/manifest."""
    if source is not None:
        root = Path(source).resolve(strict=True)
        manifest_data = regular_bytes(root / 'files-sha256.json')
        require(sha(manifest_data) == FROZEN_MANIFEST_SHA, 'historical_manifest_identity')
        manifest = json.loads(manifest_data)
        actual = {str(path.relative_to(root)) for path in root.rglob('*')
                  if (path.is_file() or path.is_symlink()) and '__pycache__' not in path.parts
                  and path.name != 'files-sha256.json'}
        require(actual == set(manifest), 'frozen_inventory')
        result = {name: regular_bytes(root / safe(name)) for name in manifest}
        require(all(sha(result[name]) == value for name, value in manifest.items()), 'frozen_bytes')
        result['files-sha256.json'] = manifest_data
        return result
    root = shared / 'recovery_v4p1'
    recipe = json.loads(regular_bytes(root / 'package_recipe325.json'))
    manifest_data = regular_bytes(root / 'files-sha256.json')
    manifest = json.loads(manifest_data)
    require(sha(manifest_data) == recipe['sealed_manifest_sha256'] == FROZEN_MANIFEST_SHA and
            set(manifest) == set(recipe['sources']), 'frozen_recipe_manifest')
    result = {}
    for name, row in recipe['sources'].items():
        kind = row['source']
        require(kind in ('recovery', 'shared', 'shared-normalize-final-newlines'), 'frozen_source_kind')
        data = regular_bytes((root if kind == 'recovery' else shared) / safe(name))
        if kind == 'shared-normalize-final-newlines':
            require(name == 'verify322.py', 'frozen_newline_exception')
            data = data.rstrip(b'\n') + b'\n\n'
        require(sha(data) == row['sha256'] == manifest[name], 'historical_source_differs:' + name)
        result[name] = data
    result['files-sha256.json'] = manifest_data
    return result


def build(shared, output, *, frozen_package=None, run_tests=True):
    shared, output = Path(shared).resolve(strict=True), Path(output).resolve()
    here = Path(__file__).resolve().parent
    require(not output.exists(), 'fresh_output_required')
    required = ('resume_state325.py', 'resume_workflow325.py', 'resume_native325.py',
                'stream_compare325.py', 'native_runner325.py', 'test_native_runner325.py')
    require(all((shared / name).is_file() for name in required), 'native_runner_and_tests_required')
    payload = {}
    # Keep actual checker closure and complete offline test/data dependencies.
    # Old launch/sbatch shell entry points are excluded; the root wrappers below
    # invoke only this new verification controller.
    for path in sorted(shared.iterdir()):
        if path.is_file() and path.suffix in ('.py', '.json', '.f90', '.c', '.in'):
            payload[path.name] = regular_bytes(path)
    for folder in ('fixture', 'grids', 'validation'):
        for path in sorted((shared / folder).rglob('*')):
            if path.is_file():
                require(not path.is_symlink() and '__pycache__' not in path.parts,
                        'regular_fixture')
                payload[str(path.relative_to(shared))] = regular_bytes(path)
    frozen = frozen_payload(shared, frozen_package)
    for name, data in frozen.items():
        payload['recovery_v4p1/' + name] = data
        # The current worker-binding test has a local-project fallback, which
        # does not exist on a clean extracted native package. Supply its exact
        # inspected source fixture from the byte-preserved historical payload.
        if name.startswith('fixture/source-binding/'):
            require(name not in payload or payload[name] == data, 'source_fixture_conflict')
            payload[name] = data
    for path in sorted(here.iterdir()):
        if path.is_file() and path.suffix in ('.py', '.sh', '.sbatch', '.md'):
            require(path.name not in payload, 'new_controller_name_collision:' + path.name)
            payload[path.name] = regular_bytes(path)
    evidence = {
        'schema': 1, 'artifact': ARTIFACT, 'source_job': 55569,
        'scope': 'verification of existing captures; zero builds/solver/performance runs',
        'independent_allocation_seconds': 1800, 'automatic_retry': False,
        'predecessor_job': 55575, 'predecessor_elapsed_seconds': 0,
        'predecessor_failure': 'native.runner325.hard_and_soft_address_limit before workflow',
        'reviewed_user_approval_required_before_launch': True,
        'historical_payload_manifest_sha256': sha(payload['recovery_v4p1/files-sha256.json']),
        'checker_source_sha256': {name: sha(data) for name, data in sorted(payload.items())
                                 if name.endswith('.py') and '/' not in name},
        'evidence_sha256': {name: sha(data) for name, data in sorted(payload.items())
                           if name.startswith('validation/') or name in (
                               'RESUME_VALIDATION.json', 'PROVIDER_PREFLIGHT.json')},
    }
    payload['native-package325.json'] = (json.dumps(evidence, sort_keys=True, indent=2) + '\n').encode()
    require(sum(map(len, payload.values())) <= MAX_PAYLOAD, 'payload_cap')
    manifest = {name: sha(data) for name, data in sorted(payload.items())}
    payload['native-files-sha256.json'] = (json.dumps(manifest, sort_keys=True, indent=2) + '\n').encode()
    output.mkdir()
    kit = output / ARTIFACT
    kit.mkdir(mode=0o700)
    for name, data in sorted(payload.items()):
        path = kit / safe(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(0o755 if name.endswith('.sh') else 0o644)
    check = [sys.executable, '-B', '-c',
             'from control_native325 import package_check; print(package_check("."))']
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    subprocess.run(check, cwd=kit, env=environment, check=True, capture_output=True, text=True)
    if run_tests:
        with (output / 'native-package-tests.log').open('wb') as log:
            subprocess.run([sys.executable, '-B', '-m', 'unittest', 'discover', '-v'],
                           cwd=kit, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
        # Tests must not mutate sealed inputs or leave unlisted regular files.
        subprocess.run(check, cwd=kit, env=environment, check=True, capture_output=True, text=True)
    archive = output / (ARTIFACT + '.tar.gz')
    with archive.open('wb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as tar:
            for name, data in sorted(payload.items()):
                info = tarfile.TarInfo(ARTIFACT + '/' + name)
                info.size, info.mode, info.mtime = len(data), 0o755 if name.endswith('.sh') else 0o644, 0
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                tar.addfile(info, io.BytesIO(data))
    archive_sha = sha(archive.read_bytes())
    Path(str(archive) + '.sha256').write_text(archive_sha + '  ' + archive.name + '\n')
    result = dict(artifact=ARTIFACT, archive=str(archive), archive_sha256=archive_sha,
                  manifest_sha256=sha(payload['native-files-sha256.json']), files=len(payload),
                  payload_bytes=sum(map(len, payload.values())), tests_requested=run_tests,
                  solver_runs=0, compiler_runs=0, performance_runs=0)
    (output / 'native-package-build.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shared', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--frozen-package', type=Path)
    parser.add_argument('--skip-tests', action='store_true')
    args = parser.parse_args()
    print(json.dumps(build(args.shared, args.output, frozen_package=args.frozen_package,
                           run_tests=not args.skip_tests), indent=2))
