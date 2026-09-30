"""Collect source-bound Stage231 local correctness evidence; no timing claim."""
import argparse
import hashlib
import json
from pathlib import Path
from capability_223 import digest, save

BASE = Path(__file__).resolve().parents[2]


def review(root):
    root = root.resolve()
    build = json.loads((root/'final-build-integrity.json').read_text())
    tree = root/'final-candidate'
    for name, expected in build['candidate_sources'].items():
        assert hashlib.sha256((tree/name).read_bytes()).hexdigest() == expected, name
    for name, expected in build['changed_sources'].items():
        assert digest(BASE/name) == expected, name
    for name, expected in build['binaries'].items():
        assert digest(tree/name) == expected, name
    for name, expected in build['build_logs'].items():
        assert digest(root/name) == expected, name
    cases = {}
    for name in ('production-4', 'oracles-4', 'production-41', 'checked-4',
                 'negative-cache', 'negative-domain', 'mixed-validation'):
        directory = root/('final-'+name)
        result = json.loads((directory/'results.json').read_text())
        assert result['status'] == 'passed', name
        for key, file in (('binary_sha256', 'climate'), ('log_sha256', 'run.log'),
                          ('input_sha256', 'simple.in'), ('injected_source_sha256', 'main.f90')):
            if key in result:
                assert digest(directory/file) == result[key], (name, file)
        if name.startswith('negative-'):
            assert result['expected_failure'] and result['returncode'] != 0
            assert result['diagnostic'] in (directory/'run.log').read_text()
        else:
            executable = 'bin-230-checked/climate' if name == 'checked-4' else 'bin-230/climate'
            assert result['binary_sha256'] == build['binaries'][executable]
            if name != 'mixed-validation':
                assert result['steps'] == 12 and result['coverage']['remaps'] == 2
                assert result['numerical']['exact_observed'] and result['step_states_exact']
                assert result['vtk_exact_members'] == 62
        cases[name] = result
    result = dict(stage=231, status='local-correctness-passed', build=build, cases=cases,
                  raw_root=str(root),
                  scope='Optimized production at 4 and 41 ranks (40 Domains), optimized full oracles and checked production at 4 ranks. Two explicit negative guards and one-step mixed-validation control.',
                  limitations='Local oversubscribed tests are not performance/scaling evidence or cluster qualification. Live Domain primitive, RK, wavelet and external-consumer dependencies remain.',
                  performance_claim=None)
    save(BASE/'docs/parallel-blocks-stage231-local-results.json', result)
    print('PASS: seven source-bound local cases; no performance claim')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=BASE/'.stage230-validation/stage231')
    review(parser.parse_args().root)
