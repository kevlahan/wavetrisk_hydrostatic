"""Real MPI negative tests: reject unused cache creation and missing reference state.

Compile only an injected main module against the already-built candidate. The
production source and executable are never edited. This executable is test-only.
"""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess

from capability_223 import digest, save
from profile_protocol import checkpoint_name, environment
from run_local_220 import SELECTORS


def build(tree, build_log, out, mode="cache"):
    out.mkdir(exist_ok=False)
    (out/'mod').mkdir()
    source = (tree/'src/main.f90').read_text()
    anchor = '         error stop "production block timestep retained Domain tendency storage"\n'
    if source.count(anchor) != 1:
        raise ValueError('Expected one storage guard')
    if mode == 'cache':
        declaration='  use parallel_block_mod, only : ensure_local_block_hydrostatic_state'
        call='    call ensure_local_block_hydrostatic_state ! Test-only forbidden cache recreation\n'
    elif mode == 'domain':
        declaration='  use parallel_block_mpi_mod, only : capture_block_domain_multistage_candidate_tendency'
        call='    call capture_block_domain_multistage_candidate_tendency(1,4,compatibility_remainder=.false.)\n'
    else:
        raise ValueError('Unknown rejection mode')
    source=source.replace('module main_mod','module main_mod\n'+declaration,1)
    source=source.replace(anchor,anchor+call)
    (out/'main.f90').write_text(source)
    lines = build_log.read_text().splitlines()
    compile = shlex.split(next(l for l in lines if ' -c src/main.f90 -o ' in l))
    object_name = compile[compile.index('-o')+1]
    compile = [str(out/'main.f90') if x == 'src/main.f90' else
               str(out/'main.o') if x == object_name else
               '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link = shlex.split(next(l for l in reversed(lines)
                           if object_name in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1] = str(out/'climate')
    link = [str(out/'main.o') if x == object_name else x for x in link]
    with (out/'build.log').open('x') as log:
        for command in (compile, link):
            subprocess.run(command, cwd=tree, stdout=log, stderr=subprocess.STDOUT, check=True)
    save(out/'commands.json', [compile, link])


def run(tree, build_log, fixture, out, mode="cache"):
    tree, build_log, fixture, out = (p.resolve() for p in (tree, build_log, fixture, out))
    original = digest(tree/'src/main.f90')
    build(tree, build_log, out, mode)
    text = (fixture/'simple.in').read_text()
    for name in ('simple.in', checkpoint_name(text)):
        shutil.copy2(fixture/name, out/name)
    (out/'grids').symlink_to((fixture/'grids').resolve(), target_is_directory=True)
    env = environment(oracles=False, detail=False)
    env.update(SELECTORS)
    env['OMPI_MCA_mpi_yield_when_idle'] = '1'
    command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', '4', './climate', 'simple.in']
    with (out/'run.log').open('x') as log:
        code = subprocess.run(command, cwd=out, env=env, stdout=log,
                              stderr=subprocess.STDOUT, timeout=180).returncode
    log = (out/'run.log').read_text()
    diagnostic = ('block hydrostatic cache requested without a numerical consumer' if mode=='cache' else
                  'non-native tendency assembly requires Domain state')
    assert code != 0 and diagnostic in log, (code, log[-1500:])
    assert 'Domain tendency storage: absent (resident block dynamics)' in log
    assert digest(tree/'src/main.f90') == original
    result = dict(status='passed', expected_failure=True, returncode=code,
                  diagnostic=diagnostic, production_source_unchanged=True,
                  injected_source_sha256=digest(out/'main.f90'),
                  binary_sha256=digest(out/'climate'), log_sha256=digest(out/'run.log'),
                  command=command)
    save(out/'results.json', result)
    print('PASS: '+diagnostic, flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('tree', 'build-log', 'fixture', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--mode', choices=('cache','domain'), default='cache')
    run(**vars(parser.parse_args()))
