"""Bind new candidate replays and exact reused tiled/main references to sealed sources."""
import json
import os
from pathlib import Path
import sys
from common import save, sha
from build import objects
from evidence import publication_evidence


def validate(kit, work, baseline, main_anchor):
    kit, work, baseline, main_anchor = [Path(p).resolve(strict=True) for p in (kit, work, baseline, main_anchor)]
    source = json.loads((kit / 'source-files.json').read_text())
    old = json.loads((baseline / 'local-validation.json').read_text())
    old_builds = json.loads((baseline / 'builds.json').read_text())
    main_validation = json.loads((main_anchor / 'local-validation.json').read_text())
    main_builds = json.loads((main_anchor / 'builds.json').read_text())
    if any(record['status'] != 'passed' for record in (old, old_builds, main_validation, main_builds)):
        raise ValueError('Require completed source-bound Stage265t and main J4 references')
    for family, tree in (('candidate', work/'source'), ('stage265t', baseline/'source'), ('main', main_anchor/'sources/main')):
        for name, expected in source[family].items():
            if sha(tree / name) != expected:
                raise ValueError('Validated/packaged source differs: ' + family + '/' + name)
    physics = objects(work / 'source')
    if physics != old_builds['physics_objects'] or physics != main_builds['physics_objects']:
        raise ValueError('Candidate physics differs from reference builds')
    notes = dict(status='passed', host=old_builds['host'], level='J4', builds={}, physics_objects=physics,
                 source_manifest_sha256=sha(kit / 'source-files.json'), identity=json.loads((kit / 'IDENTITY.json').read_text()),
                 baseline_build_manifest_sha256=sha(baseline / 'builds.json'),
                 baseline_validation_sha256=sha(baseline / 'local-validation.json'),
                 main_build_manifest_sha256=sha(main_anchor / 'builds.json'),
                 main_validation_sha256=sha(main_anchor / 'local-validation.json'),
                 limitation='Exact source-bound Stage265t/main reference builds reused; both candidate builds and three replays are new.')
    cases = {}
    for family, directory, build_notes, previous, key in (
            ('stage265t', baseline, old_builds, old, 'candidate'),
            ('main', main_anchor, main_builds, main_validation, 'main')):
        row = build_notes['builds'][key].copy()
        binary = directory / row['binary']
        if sha(binary) != row['binary_sha256']:
            raise ValueError('Reference binary changed')
        row['binary'] = os.path.relpath(binary, work)
        row['reused_reference'] = True
        notes['builds'][family] = row
        cases[family] = previous['cases'][key].copy()
        if family == 'stage265t':
            log = baseline / 'opt/run.log'
            if sha(log) != cases[family]['log_sha256']:
                raise ValueError('Baseline replay log changed')
            cases[family]['publication'] = publication_evidence(log.read_text(), family, False)
    for kind, suffix in (('candidate', 'opt'), ('checked', 'checked')):
        binary = work / 'source' / ('bin-rk266-' + suffix) / 'climate'
        command = ['make', '-f', 'Makefile', '-f', 'stage264-flags.mk', '-j1', 'TEST_CASE=climate',
                   'PARAM=param_J4', 'DEBUG=' + ('check' if kind == 'checked' else 'false'),
                   'BUILD_DIR=build-rk266-' + suffix, 'BIN_DIR=bin-rk266-' + suffix,
                   'FP_CONTRACT=off', 'STAGE202=on', 'STAGE203=on']
        notes['builds'][kind] = dict(command=command, binary=str(binary.relative_to(work)), binary_sha256=sha(binary),
                                   physics_objects=physics, log_sha256=sha(work / ('build-' + suffix + '.log')))
    for label, kind, oracle in (('opt', 'candidate', False), ('checked', 'checked', False), ('oracle', 'candidate', True)):
        row = json.loads((work / label / 'results.json').read_text())
        if row['status'] != 'passed' or row['ranks'] != 4 or row['oracles'] != oracle:
            raise ValueError('Incomplete real replay: ' + label)
        if row['binary_sha256'] != notes['builds'][kind]['binary_sha256']:
            raise ValueError('Replayed binary differs')
        log = work / label / 'run.log'
        if sha(log) != row['log_sha256']:
            raise ValueError('Replay log changed')
        row['publication'] = publication_evidence(log.read_text(), kind, oracle)
        cases['candidate-oracle' if oracle else kind] = row
    save(work / 'builds.json', notes)
    result = dict(status='passed', ranks=4, cases=cases, source_manifest_sha256=sha(kit / 'source-files.json'),
                  build_manifest_sha256=sha(work / 'builds.json'),
                  limitation='Local four-rank correctness only; cluster performance not yet measured.')
    save(work / 'local-validation.json', result)
    print('PASS: exact optimized/checked/oracle candidate replay, source identity and shared physics objects')
    return result


if __name__ == '__main__':
    validate(*sys.argv[1:])
