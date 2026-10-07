"""Four-core source-bound correctness replay of the packaged Stage230/main build."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))
from common import runtime, save, sha
from capability_223 import coverage, numerical
from run_local_230 import run as block_run, vtk_members


def validate(kit, out, reference):
    kit, out, reference = [Path(p).resolve(strict=True) for p in (kit, out, reference)]
    notes = json.loads((out / 'builds.json').read_text())
    if notes['status'] != 'passed' or notes['level'] != 'J4':
        raise ValueError('Require the three real packaged J4 builds')
    sources = json.loads((kit / 'source-files.json').read_text())
    for family in ('stage230', 'main'):
        for name, expected in sources[family].items():
            if sha(out / 'sources' / family / name) != expected:
                raise ValueError('Local built source differs: ' + family + '/' + name)
    fixture = kit / 'fixture-J4'
    result = dict(status='running', ranks=4, build_manifest_sha256=sha(out / 'builds.json'),
                  source_manifest_sha256=sha(kit / 'source-files.json'), cases={})
    for kind, oracle in (('stage230', False), ('checked', False), ('stage230', True)):
        label = kind + ('-oracle' if oracle else '')
        binary = out / notes['builds'][kind]['binary']
        print('Local four-rank replay: ' + label, flush=True)
        case = block_run(binary, fixture, out / label, reference, ranks=4, oracles=oracle)
        case['binary_sha256'] = sha(binary)
        result['cases'][label] = case
        save(out / 'local-validation.json', result)
    directory = out / 'main'
    directory.mkdir()
    for name in ('simple.in', 'localJ4J6_checkpoint_0002.bin.zst'):
        shutil.copy2(fixture / name, directory / name)
    shutil.copytree(fixture / 'grids', directory / 'grids')
    binary = out / notes['builds']['main']['binary']
    shutil.copy2(binary, directory / 'climate')
    print('Local four-rank replay: main', flush=True)
    command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', '4', './climate', 'simple.in']
    with (directory / 'run.log').open('x') as log:
        code = subprocess.run(command, cwd=directory, env=runtime(), stdout=log,
                              stderr=subprocess.STDOUT, timeout=1800).returncode
    analysis, cov = coverage(directory, 'J4', 'localJ4J6')
    gate = numerical(directory, out / 'stage230', cov['new_checkpoint'], 40, False)
    states = [r['state'] for r in analysis['step_records']]
    from analyze import analyze
    expected = [r['state'] for r in analyze(out / 'stage230/run.log')['step_records']]
    if code or len(states) != 12 or states != expected or cov['remaps'] != 2 or not gate['passed']:
        raise ValueError('Local main-relative science gate failed')
    vtk = vtk_members(directory)
    if len(vtk) != 62:
        raise ValueError('Incomplete main VTK outputs')
    result['cases']['main'] = dict(returncode=code, states=states, remaps=2,
         vtk_members=62, numerical=gate, binary_sha256=sha(binary), checkpoint_sha256=sha(directory / cov['new_checkpoint']))
    result.update(status='passed', limitation='Four-rank local correctness only; no performance/scaling claim.')
    save(out / 'local-validation.json', result)
    print('PASS: optimized/checked/oracle Stage230 exact replay and unchanged main numerical gate', flush=True)
    return result


if __name__ == '__main__':
    validate(*sys.argv[1:])
