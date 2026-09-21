"""Validate cleanup trajectories, discarded-state poison, traces and ownership probe."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct
import tempfile

import numpy as np
from checkpoint_gate import accept, POLICY
from cluster_199 import workload_gate
from compare_restart_checkpoints import compare
from local_j4 import compare_windows
from run_local_220 import SELECTORS
from audit_cleanup_222 import instrument as instrument_cleanup
from probe_domainless_222 import instrument as instrument_domainless


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main_before_empty_cleanup(workspace):
    text = (workspace/'src/main.f90').read_text()
    anchor = ('    ! Empty-Domain ranks still own shared arrays and zero-sized field stores.\n'
              '    ! Their Domain loops are empty, but the common restart cleanup must run.\n')
    if text.count(anchor) != 1:
        raise ValueError('Unexpected empty-rank cleanup change')
    return text.replace(anchor, '    if (size(grid) <= 0) return\n')


def trace_inventory(root):
    rows = {}
    for p in sorted(root.glob('cleanup222-*.bin')):
        records = values = 0
        with p.open('rb') as stream:
            while prefix := stream.read(4):
                if len(prefix) != 4:
                    raise ValueError('Truncated trace count')
                n, = struct.unpack('=i', prefix)
                if n < 0:
                    raise ValueError('Invalid trace extent')
                raw = stream.read(8*n)
                if len(raw) != 8*n or not np.isfinite(np.frombuffer(raw, dtype='=f8')).all():
                    raise ValueError('Truncated/nonfinite trace: '+str(p))
                records += 1
                values += n
        rows[p.name] = dict(bytes=p.stat().st_size, records=records, values=values, sha256=digest(p))
    expected = {f'cleanup222-{label}-{rank}.bin' for rank in range(4) for label in
                ('compact_scalar', 'compact_vector', 'inverse', 'rk', 'forward',
                 'mass_tendency', 'velocity_tendency')}
    if rows.keys() != expected or any(not v['values'] for v in rows.values()):
        raise ValueError('Incomplete trace inventory')
    return rows


def probe_source_identity(workspace, build, domainless):
    """Reproduce the probe edits instead of trusting arbitrary instrumented sources."""
    names = ('parallel_block.f90', 'parallel_block_build.f90', 'parallel_block_mpi.f90',
             'parallel_block_inverse.f90', 'time_integr.f90', 'arch.f90', 'checkpoint.f90', 'main.f90')
    with tempfile.TemporaryDirectory(prefix='stage222-source-check-') as tmp:
        tree = Path(tmp)
        (tree/'src').mkdir()
        for name in names:
            shutil.copy2(workspace/'src'/name, tree/'src'/name)
        if domainless:
            instrument_domainless(tree)
        else:
            instrument_cleanup(tree, poison=True)
            (tree/'src/main.f90').write_text(main_before_empty_cleanup(workspace))
        for name in names:
            if name == 'checkpoint.f90' and not domainless:
                continue  # The bounded cleanup trace predates this separate fix.
            if digest(tree/'src'/name) != digest(build/'src'/name):
                raise ValueError('Instrumented source differs from reproducible probe: '+name)


def validate(root, workspace, output):
    baseline = Path('/private/tmp/stage215-local/optimized-final')
    ref = json.loads((baseline/'results.json').read_text())
    result = dict(stage=222, status='passed', performance_claim=None,
                  scope='Four-task J4/J6 Z30 adaptive/restart correctness; not an 83-task timing or >Domain-count release.',
                  baseline=str(baseline), runs={})
    result['baseline_record_sha256'] = digest(baseline/'results.json')
    result['validation_helpers'] = {n: digest(Path(__file__).parent/n) for n in
                                    ('audit_cleanup_222.py', 'probe_domainless_222.py',
                                     'run_local_222.py', 'validate_cleanup_222.py',
                                     'compare_restart_checkpoints.py', 'test_restart_checkpoints.py')}
    checked_logs = ('build-final-checked.log', 'build-domainless-release3-checked.log')
    for name in checked_logs:
        log = (root/name).read_text()
        if any(flag not in log for flag in
               ('-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffp-contract=off')):
            raise ValueError('Missing checked build flags: '+name)
    result['checked_build_logs'] = {n: digest(root/n) for n in checked_logs}
    source_names = ('src/parallel_block.f90', 'src/parallel_block_mpi.f90',
                    'src/parallel_block_build.f90', 'src/checkpoint.f90', 'src/main.f90')
    result['workspace_sources'] = {n: digest(workspace/n) for n in source_names}
    for name in ('release-optimized', 'final-checked', 'final-oracles', 'domainless-restart-final'):
        p = root/name
        run = json.loads((p/'results.json').read_text())
        log = (p/'run.log').read_text()
        if (run['status'] != 'passed' or run['ranks'] != 4 or
                run['reloads'] != [2, 3] or run['remaps'] != 2):
            raise ValueError('Incomplete run: '+name)
        if len(run['analysis']['step_records']) != 12:
            raise ValueError('Wrong step count: '+name)
        for key in ('parameters', 'input_sha256', 'seed_sha256', 'grids', 'ranks'):
            if run[key] != ref[key]:
                raise ValueError('Workload identity differs: '+name+' '+key)
        if any(run['environment'].get(k) != v for k, v in SELECTORS.items()):
            raise ValueError('Wrong accepted runtime selectors')
        if (run['environment']['WAVETRISK_VALIDATE_BLOCK_DYNAMICS'] == '1') != (name == 'final-oracles'):
            raise ValueError('Wrong oracle mode')
        if digest(p/'climateJ4') != run['binary_sha256'] or digest(Path(run['binary'])) != run['binary_sha256']:
            raise ValueError('Binary identity changed')
        for n, h in run['source_sha256'].items():
            if digest(Path(run['build'])/n) != h:
                raise ValueError('Run source changed: '+name+' '+n)
        if name in ('release-optimized', 'final-oracles'):
            if any(run['source_sha256'][n] != h for n, h in result['workspace_sources'].items()
                   if n != 'src/main.f90'):
                raise ValueError('Final solver source differs from tested source')
            if (Path(run['build'])/'src/main.f90').read_text() != main_before_empty_cleanup(workspace):
                raise ValueError('Unexpected normal-layout main source difference')
        eq = compare_windows(baseline, p, 'localJ4J6_checkpoint_0002.bin.zst', expected_new=[3])
        work = workload_gate(ref['analysis'], run['analysis'])
        if not work['passed']:
            raise ValueError('Adaptive trajectory differs: '+name)
        cp = compare(baseline/'localJ4J6_checkpoint_0003.bin.zst',
                     p/'localJ4J6_checkpoint_0003.bin.zst', expected_domains=40)
        failures = accept(cp, exact=name != 'domainless-restart-final')
        if failures or (name != 'domainless-restart-final' and not eq['passed']):
            raise ValueError('Checkpoint mismatch: '+name+' '+repr(failures))
        row = dict(path=str(p), binary_sha256=run['binary_sha256'],
                   record_sha256=digest(p/'results.json'), log_sha256=digest(p/'run.log'),
                   exact_checkpoint=not accept(cp, exact=True), checkpoint=cp,
                   workload=work, reloads=run['reloads'], remaps=run['remaps'],
                   sources={n: run['source_sha256'][n] for n in source_names})
        if 'postprocess_recovery' in run:
            row['postprocess_recovery'] = run['postprocess_recovery']
        if name == 'final-checked':
            probe_source_identity(workspace, Path(run['build']), domainless=False)
            row['scope'] = ('Checked cleanup and intermediate equivalence; predates the '
                            'empty-rank restart fixes, checked by domainless-restart-final.')
            inverse = re.findall(r'Stage222 poisoned inverse rank = (\d+) slots = (\d+)', log)
            staging = re.findall(r'Stage222 poisoned staging rank = (\d+)', log)
            grid = re.findall(r'Stage222 skipped grid publication verified rank = (\d+)', log)
            for rank in range(4):
                if sum(int(r) == rank and int(n) > 0 for r, n in inverse) != 36:
                    raise ValueError('Missing inverse poison coverage')
                if staging.count(str(rank)) != 48 or grid.count(str(rank)) != 12:
                    raise ValueError('Missing staging/grid proof coverage')
            row['deletion_checks'] = dict(inverse_publications_poisoned=len(inverse),
                                         temperature_transactions_poisoned=len(staging),
                                         omitted_grid_publications_verified=len(grid))
        if name == 'domainless-restart-final':
            probe_source_identity(workspace, Path(run['build']), domainless=True)
            if run['source_sha256']['src/checkpoint.f90'] != result['workspace_sources']['src/checkpoint.f90']:
                raise ValueError('Domainless run lacks current checkpoint fix')
            owned = re.findall(r'Stage222 domainless owner rank = (\d+) blocks = (\d+)', log)
            if len(owned) < 2 or any(int(r) != 3 or int(n) <= 0 for r, n in owned):
                raise ValueError('Missing domainless ownership/rebuild evidence')
            row['rank3_owned_block_counts'] = list(map(lambda x: int(x[1]), owned))
            grid = re.findall(r'Stage222 discarded grid buffers rank = (\d+)', log)
            if any(grid.count(str(rank)) != 12 for rank in range(4)):
                raise ValueError('Missing discarded grid-buffer poison coverage')
            row['discarded_grid_buffers_poisoned'] = len(grid)
            row['comparison_policy'] = POLICY
        result['runs'][name] = row
    left = trace_inventory(root/'reference-trace')
    right = trace_inventory(root/'final-checked')
    if left != right:
        raise ValueError('Intermediate streams differ from frozen Stage215')
    result['exact_intermediate_streams'] = right
    result['intermediate_trace_bytes'] = sum(v['bytes'] for v in right.values())
    result['notes'] = ['Startup ntasks<=N_GLO_DOMAIN guard remains unchanged.',
                       'Domainless probe changes geometry ownership only in a disposable test tree.',
                       'Checked cleanup traces precede the restart fix; the domainless checked run tests that fix.',
                       'Normal-layout runs precede the empty-rank early-return removal; this exact sole main.f90 difference is verified.',
                       'No cluster performance gain or >160-task capability is claimed.']
    output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: cleanup correctness, exact intermediates, discarded-state poison, and domainless owned-block trajectory.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--workspace', type=Path, default=Path.cwd())
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    validate(args.root, args.workspace, args.output)
