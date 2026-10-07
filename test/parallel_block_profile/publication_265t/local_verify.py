"""Bind real candidate replays and sealed reused Stage264 references to kit sources."""
import json
import os
from pathlib import Path
import sys
from common import save, sha
from build import objects
from evidence import publication_evidence


def validate(kit, work, anchor):
    kit, work, anchor = [Path(p).resolve(strict=True) for p in (kit, work, anchor)]
    source = json.loads((kit / 'source-files.json').read_text())
    previous = json.loads((anchor / 'local-validation.json').read_text())
    prior_builds = json.loads((anchor / 'builds.json').read_text())
    if previous['status'] != 'passed' or prior_builds['status'] != 'passed' or prior_builds['level'] != 'J4':
        raise ValueError('Require completed source-bound Stage264 J4 references')
    for family in ('stage230', 'main', 'candidate'):
        tree = work / 'source' if family == 'candidate' else anchor / 'sources' / family
        for name, expected in source[family].items():
            if sha(tree / name) != expected:
                raise ValueError('Validated/packaged source differs: ' + family + '/' + name)
    physics = objects(work / 'source')
    if physics != prior_builds['physics_objects']:
        raise ValueError('Candidate physics differs from reference builds')
    notes = dict(status='passed', host=prior_builds['host'], level='J4', builds={}, physics_objects=physics,
                 source_manifest_sha256=sha(kit / 'source-files.json'), identity=json.loads((kit / 'IDENTITY.json').read_text()),
                 reference_build_manifest_sha256=sha(anchor / 'builds.json'),
                 reference_validation_sha256=sha(anchor / 'local-validation.json'),
                 limitation='Sealed Stage264 reference builds reused; both candidate builds and three replays are new.')
    cases = {}
    for kind in ('stage230', 'main'):
        row = prior_builds['builds'][kind].copy()
        binary = anchor / row['binary']
        if sha(binary) != row['binary_sha256']:
            raise ValueError('Reference binary changed')
        row['binary'] = os.path.relpath(binary, work)
        row['reused_reference'] = True
        notes['builds'][kind] = row
        cases[kind] = previous['cases'][kind]
    for kind, suffix in (('candidate', 'opt'), ('checked', 'checked')):
        binary = work / 'source' / ('bin-rk265-' + suffix) / 'climate'
        command = ['make', '-f', 'Makefile', '-f', 'stage264-flags.mk', '-j1', 'TEST_CASE=climate',
                   'PARAM=param_J4', 'DEBUG=' + ('check' if kind == 'checked' else 'false'),
                   'BUILD_DIR=build-rk265-' + suffix, 'BIN_DIR=bin-rk265-' + suffix,
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
