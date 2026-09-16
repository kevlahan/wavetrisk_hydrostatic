#!/usr/bin/env python3
"""Verify recorded Stage 205 build identities before entering the cluster campaign."""
import json
from pathlib import Path
import subprocess
import sys
from build_restriction_205 import digest, verify, verify_helpers


def submit(kit, tree, fixture, out):
    kit, tree = kit.resolve(strict=True), tree.resolve(strict=True)
    verify_helpers(kit)
    manifest = json.loads((kit / 'stage205-sources.json').read_text())
    identity_path = tree / 'restriction-205-builds.json'
    identity = json.loads(identity_path.read_text())
    if identity.get('status') != 'passed' or identity.get('fp_contract') != 'off':
        raise ValueError('Six successful contraction-off builds are required')
    if (identity.get('baseline_sources') != manifest['baseline_sources'] or
            identity.get('reference_sources') != manifest['reference_sources'] or
            identity.get('candidate_sources') != manifest['candidate_sources'] or
            identity.get('patch_sha256') != manifest['candidate_patch_sha256']):
        raise ValueError('Build record does not match this kit')
    verify(tree, manifest['candidate_sources'])
    # The final block build is the checked, uninstrumented Stage203 control.
    if digest(tree / 'Makefile') != manifest['candidate_sources']['Makefile']:
        raise ValueError('Makefile changed since the builds')
    command = [sys.executable, str(kit / 'cluster_205.py')]
    for kind in ('baseline','baseline_profile','candidate','checked','legacy','legacy_profile'):
        record = identity['builds'][kind]
        binary = tree / record['binary']
        if not kind.startswith('legacy'):
            source_key={'baseline':'reference_sources','baseline_profile':'baseline_profile_sources',
                        'candidate':'candidate_sources','checked':'candidate_sources'}[kind]
            if record.get('sources')!=manifest[source_key]:raise ValueError('Build source identity differs: '+kind)

        if kind == 'legacy_profile':
            if record.get('sources')!=manifest['legacy_profile_sources']:raise ValueError('Legacy observer identity differs')
            verify(Path(record['source_root']),manifest['legacy_profile_sources'])
        if kind == 'legacy':
            if record.get('sources') != manifest['legacy_io_sources']:
                raise ValueError('Legacy source identity differs from kit')
            if record.get('io_patch_sha256') != manifest['legacy_io_patch_sha256']:
                raise ValueError('Legacy I/O patch identity differs')
            verify(Path(record['source_root']), manifest['legacy_io_sources'])
        if digest(binary) != record['sha256']:
            raise ValueError('Executable changed since build: ' + kind)
        command.extend(['--' + kind, str(binary)])
    command.extend(['--fixture', str(fixture.resolve(strict=True)), '--out', str(out.resolve()),
                    '--build-notes-file', str(identity_path)])
    subprocess.run(command, check=True)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit('Usage: submit_restriction_205.py KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR')
    submit(*(Path(p) for p in sys.argv[1:]))
