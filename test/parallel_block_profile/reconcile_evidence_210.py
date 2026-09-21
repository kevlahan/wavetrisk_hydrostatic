#!/usr/bin/env python3
"""Reconcile existing Stage206–210 evidence; never run or modify the solver.

Inputs are extracted review archives. Verify each manifest before consuming raw
logs. Keep historical checkpoint-only and current startup/checkpoint selections
separate. Never pool timesteps as independent repetitions or sum nested costs.
"""
import argparse
import csv
import hashlib
import importlib
import json
import math
import statistics as st
import subprocess
from pathlib import Path

from analyze import analyze
from analyze_operations_199 import common, CONTEXTS
from protocol_210 import ordinary

JOBS = {206: 54702, 207: 54704, 208: 54706, 209: 54708, 210: 54710}
GROUPS = {'wavelets': [1], 'RK tendencies': [2], 'adaptation/reconstruction': [3, 4],
          'other': [0, 5, 6, 7, 8]}
SELECTIONS = ('historical', 'screen_protocol', 'common')


def read(path):
    return json.loads(path.read_text())


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def close(a, b, tol=1e-6):
    if abs(a-b) > tol:
        raise ValueError(f'Accounting mismatch: {a} != {b}')


def verified(root, stage):
    manifest = read(root / 'collection.json')
    assert int(manifest['stage']) == stage and str(manifest['job']) == str(JOBS[stage])
    for name, expected in manifest['files'].items():
        rel = Path(name)
        assert not rel.is_absolute() and '..' not in rel.parts
        path = root / rel
        assert path.is_file() and not path.is_symlink()
        assert path.stat().st_size == expected['bytes'] and digest(path) == expected['sha256'], name
    return {'manifest_sha256': digest(root/'collection.json'),
            'verified_files': len(manifest['files']), 'root': str(root)}


def strict_selection(records):
    # In this archived workload write_checkpoint calls restart before returning
    # to print_log. The checkpoint row is not the next full integration step.
    excluded = {1}
    for n, row in enumerate(records, 1):
        if row['checkpoint']:
            excluded.update((n, n+1))
    return [n for n in range(1, len(records)+1) if n not in excluded]


def stats(values):
    return {'mean': st.mean(values), 'median': st.median(values), 'values': values}


def rounding_interval(values, statistic):
    # stdout uses ES8.2: three significant digits, not two decimal places.
    radius = [0.5 * 10 ** (math.floor(math.log10(x))-2) for x in values]
    return [statistic([x-d for x, d in zip(values, radius)]),
            statistic([x+d for x, d in zip(values, radius)])]


def workflow(directory, ranks, selections, analysis):
    data = {rank: common(directory/f'workflow195-rank-{rank}.txt') for rank in range(ranks)}
    assert all(list(rows) == list(range(1, 13)) for rows in data.values())
    for rows in data.values():
        for n, row in rows.items():
            assert bool(row['restart']) == analysis['step_records'][n-1]['checkpoint']
            assert bool(row['remap']) == analysis['step_records'][n-1]['remap']
    out = {}
    for key, indices in selections.items():
        steps = []
        for n in indices:
            wall = [data[r][n]['wall'] for r in range(ranks)]
            slow = max(range(ranks), key=lambda r: (wall[r], -r))
            contexts = {str(c): st.mean(data[r][n]['context'][c] for r in range(ranks)) for c in range(9)}
            close(sum(contexts.values()), st.mean(wall))
            steps.append({'sequence': n, 'rank_mean': st.mean(wall), 'rank0': wall[0],
                          'rank_max': wall[slow], 'slow_rank': slow, 'contexts': contexts})
        contexts = {str(c): st.mean(q['contexts'][str(c)] for q in steps) for c in range(9)}
        out[key] = {'steps': steps, 'contexts': contexts,
                    'rank_mean': st.mean(q['rank_mean'] for q in steps),
                    'rank0_mean': st.mean(q['rank0'] for q in steps),
                    'rank0_median': st.median(q['rank0'] for q in steps),
                    'mean_rank_max': st.mean(q['rank_max'] for q in steps)}
    return out, data


def operations(directory, stage, ranks, selections, wf):
    regions = importlib.import_module(f'instrument_mechanism_{stage}').REGIONS
    nops = len(regions)
    sums = {k: {} for k in selections}
    max_error = 0.
    for rank in range(ranks):
        steps = {}
        for line in (directory/f'operations{stage}-rank-{rank}.txt').read_text().splitlines():
            f = line.split()
            if f[0] == 'step':
                n = int(f[1]); steps[n] = {'self': {}, 'mpi': {}, 'inclusive': {}}
                close(float(f[4]), wf[rank][n]['wall'])
            elif f[0] == 'self':
                n, c, p = map(int, f[1:4]); key = (c, p)
                assert key not in steps[n]['self']
                steps[n]['self'][key] = float(f[5]); steps[n]['inclusive'][key] = float(f[6])
            elif f[0] == 'mpi':
                n, c, p, owner = map(int, f[1:5]); key = (c, p, owner)
                assert key not in steps[n]['mpi']
                steps[n]['mpi'][key] = float(f[5])
        assert list(steps) == list(range(1, 13))
        for n, row in steps.items():
            assert len(row['self']) == 9*(nops+1)
            for c in range(9):
                error = abs(sum(row['self'][c, p] for p in range(nops+1))-wf[rank][n]['context'][c])
                max_error = max(max_error, error); assert error <= 1e-6
                for p in (26, 27, 28, 29):
                    close(sum(v for (cc, pp, owner), v in row['mpi'].items() if (cc, pp) == (c, p)), row['self'][c, p])
            for selection, indices in selections.items():
                if n not in indices:
                    continue
                target = sums[selection]
                for (c, p), value in row['self'].items():
                    # MPI attributed once to its immediate owner below.
                    if p in (26, 27, 28, 29):
                        continue
                    q = target.setdefault((c, p), {'local': 0., 'mpi': 0., 'inclusive': 0.})
                    q['local'] += value; q['inclusive'] += row['inclusive'][c, p]
                for (c, p, owner), value in row['mpi'].items():
                    q = target.setdefault((c, owner), {'local': 0., 'mpi': 0., 'inclusive': 0.})
                    q['mpi'] += value
    out = {}
    for selection, data in sums.items():
        divisor = ranks*len(selections[selection])
        out[selection] = [{'context': c, 'operation': p, 'name': regions.get(p, 'unselected'),
                           **{k: v/divisor for k, v in row.items()},
                           'exclusive_with_owned_mpi': (row['local']+row['mpi'])/divisor}
                          for (c, p), row in sorted(data.items())]
        close(sum(q['exclusive_with_owned_mpi'] for q in out[selection]),
              st.mean(wf[r][n]['wall'] for r in range(ranks) for n in selections[selection]))
    return out, max_error


def phase_ledger(rows, phase, selection):
    result = {}
    for group, contexts in GROUPS.items():
        sides = {side: st.mean(sum(q['workflow'][selection]['contexts'][str(c)] for c in contexts)
                              for q in rows if q['side'] == side and q['phase'] == phase)
                 for side in ('baseline', 'legacy')}
        result[group] = {**sides, 'difference': sides['baseline']-sides['legacy']}
    result['total'] = {k: sum(q[k] for q in result.values()) for k in ('baseline', 'legacy', 'difference')}
    return result


def campaign(root, stage):
    integrity = verified(root, stage)
    r = read(root/'results/results.json')
    assert r['status'] == 'passed'
    records, csvrows = [], []
    for label, record in r['runs'].items():
        directory = root/'results'/label
        text = (directory/'run.log').read_text()
        a = analyze(directory/'run.log')
        stored = read(directory/'analysis.json')
        assert {k: v for k, v in a.items() if k != 'log'} == {k: v for k, v in stored.items() if k != 'log'}
        assert a['completed'] and len(a['step_records']) == 12
        current = ordinary(text, a)
        selections = {'historical': [i+1 for i, x in enumerate(a['step_records']) if not x['checkpoint']],
                      'screen_protocol': [x['index']+1 for x in current['included']],
                      'common': strict_selection(a['step_records'])}
        assert selections['historical'] == [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12]
        assert selections['screen_protocol'] == [2, 3, 4, 5, 6, 7, 9, 10, 11, 12]
        assert selections['common'] == [2, 3, 4, 5, 6, 7, 10, 11, 12]
        timing = {key: stats([a['step_records'][n-1]['seconds'] for n in indices]) for key, indices in selections.items()}
        close(timing['historical' if stage <= 208 else 'screen_protocol']['median'], record['ordinary_step_median'], 1e-12)
        side = record.get('side', record['kind'])
        q = {'label': label, 'phase': record['phase'], 'side': side, 'round': record.get('round'),
             'binary_sha256': record['binary_sha256'], 'timing': timing, 'selections': selections,
             'states': [x['state'] for x in a['step_records']], 'rounding_bounds': {},
             'run_log_sha256': digest(directory/'run.log')}
        for key in selections:
            q['rounding_bounds'][key] = {stat: rounding_interval(timing[key]['values'], getattr(st, stat)) for stat in ('median', 'mean')}
        if record['phase'] in ('coarse', 'operations'):
            q['workflow'], rawwf = workflow(directory, record['ranks'], selections, a)
            if record['phase'] == 'operations':
                q['operations'], q['maximum_conservation_error'] = operations(directory, stage, record['ranks'], selections, rawwf)
        for i, step in enumerate(a['step_records'], 1):
            csvrows.append({'stage': stage, 'run': label, 'phase': record['phase'], 'sequence': i,
                            'state': step['state'], 'printed_seconds': step['seconds'],
                            'historical': i in selections['historical'], 'screen_protocol': i in selections['screen_protocol'], 'common': i in selections['common'],
                            'checkpoint': step['checkpoint'], 'remap': step['remap']})
        records.append(q)
    assert len({tuple(q['states']) for q in records}) == 1
    out = {'stage': stage, 'integrity': integrity, 'baseline_sources': r['build_notes']['builds']['baseline']['sources'],
           'inputs': {k: r[k] for k in ('input_sha256', 'seed_sha256', 'grids')},
           'cases': r['cases'], 'runs': records, 'measurement_status': r['measurement_status']}
    if stage <= 208:
        out['ledgers'] = {key: {phase: phase_ledger(records, phase, key) for phase in ('coarse', 'operations')} for key in SELECTIONS}
        chain = []
        def link(name, func):
            sides = {side: func([q for q in records if q['side'] == side]) for side in ('baseline', 'legacy')}
            gap = sides['baseline']-sides['legacy']
            chain.append({'name': name, **sides, 'gap': gap, 'change_from_previous': gap-chain[-1]['gap'] if chain else 0.})
        link('historical bare median', lambda rows: st.median(q['timing']['historical']['median'] for q in rows if q['phase'] == 'bare'))
        link('previous screen-protocol bare median', lambda rows: st.median(q['timing']['screen_protocol']['median'] for q in rows if q['phase'] == 'bare'))
        link('common bare median', lambda rows: st.median(q['timing']['common']['median'] for q in rows if q['phase'] == 'bare'))
        link('common bare mean', lambda rows: st.mean(q['timing']['common']['mean'] for q in rows if q['phase'] == 'bare'))
        link('common detailed printed mean', lambda rows: st.mean(q['timing']['common']['mean'] for q in rows if q['phase'] == 'operations'))
        link('common detailed rank0 workflow mean', lambda rows: st.mean(q['workflow']['common']['rank0_mean'] for q in rows if q['phase'] == 'operations'))
        link('common detailed all-rank mean', lambda rows: st.mean(q['workflow']['common']['rank_mean'] for q in rows if q['phase'] == 'operations'))
        close(chain[-1]['gap'], out['ledgers']['common']['operations']['total']['difference'])
        out['statistic_bridge'] = chain
        out['operation_differences'] = {}
        for key in SELECTIONS:
            lookup = {}
            for side in ('baseline', 'legacy'):
                for q in records:
                    if q['phase'] != 'operations' or q['side'] != side:
                        continue
                    for row in q['operations'][key]:
                        cell = lookup.setdefault((row['context'], row['operation']), {'name': row['name'], 'baseline': [], 'legacy': []})
                        cell[side].append(row)
            differences = []
            for (c, p), row in lookup.items():
                sides = {side: {metric: st.mean(x[metric] for x in row[side]) if row[side] else 0.
                                for metric in ('local', 'mpi', 'exclusive_with_owned_mpi', 'inclusive')}
                         for side in ('baseline', 'legacy')}
                differences.append({'context': c, 'operation': p, 'name': row['name'], **sides,
                                    'difference': sides['baseline']['exclusive_with_owned_mpi']-sides['legacy']['exclusive_with_owned_mpi']})
            close(sum(x['difference'] for x in differences), out['ledgers'][key]['operations']['total']['difference'])
            out['operation_differences'][key] = sorted(differences, key=lambda x: -x['difference'])
        out['calibration_common'] = []
        for side in ('baseline', 'legacy'):
            for rep in (1, 2):
                control = next(q for q in records if (q['side'], q['phase'], q['round']) == (side, 'coarse', rep))
                detail = next(q for q in records if (q['side'], q['phase'], q['round']) == (side, 'operations', rep))
                for group, contexts in {'whole': list(range(9)), 'dynamics outside WT': [0, 2], **{k: v for k, v in GROUPS.items() if k != 'other'}, 'setup': [5], 'physics': [6], 'remap': [7], 'output': [8]}.items():
                    values = [sum(q['workflow']['common']['contexts'][str(c)] for c in contexts) for q in (control, detail)]
                    change = values[1]-values[0]
                    out['calibration_common'].append({'side': side, 'round': rep, 'group': group,
                                                     'control': values[0], 'detail': values[1], 'change': change,
                                                     'flag': abs(change) > max(.001, .02*values[0])})
        # Preserve the historical decision; changed selection cannot waive flags.
        out['historical_calibration'] = r['summary']['calibration']
        close(out['ledgers']['historical']['operations']['total']['difference'], r['summary']['instrumented_rank_mean_excess'])
        if stage == 208:
            out['migration'] = migration(root/'results', records)
            groups = {'rebuild infrastructure': {2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,24,30,31,32,33},
                      'dynamic representation and ownership': {1,13,21,22,23,25}}
            out['mechanism_groups_common'] = {name: {side: sum(row[side]['exclusive_with_owned_mpi'] for row in out['operation_differences']['common'] if row['operation'] in ops) for side in ('baseline','legacy')} for name, ops in groups.items()}
    else:
        a = next(q for q in records if q['label'].endswith('screen-A'))
        b = next(q for q in records if q['label'].endswith('screen-B'))
        out['original_decision'] = r['decisions']
        out['screen'] = {'baseline': a['timing']['common']['median'], 'candidate': b['timing']['common']['median'],
                         'gain': a['timing']['common']['median']-b['timing']['common']['median'],
                         'gain_rounding_interval': [a['rounding_bounds']['common']['median'][0]-b['rounding_bounds']['common']['median'][1],
                                                    a['rounding_bounds']['common']['median'][1]-b['rounding_bounds']['common']['median'][0]]}
    return out, csvrows


def migration(results, records):
    """Catalog signatures are replicated; sent payloads are rank-local sums."""
    repeats = []
    for record in records:
        if record['phase'] != 'operations' or record['side'] != 'baseline':
            continue
        steps = {}
        for rank in range(83):
            seen = set()
            for line in (results/record['label']/f'lifecycle208-rank-{rank}.txt').read_text().splitlines()[1:]:
                f = line.split(); n, kind = int(f[2]), int(f[6])
                if kind not in (4, 5, 50, 51) or n not in record['selections']['historical']:
                    continue
                assert (n, kind) not in seen; seen.add((n, kind))
                abc = list(map(int, f[9:12])); step = steps.setdefault(n, {})
                if kind in (50, 51):
                    if kind in step:
                        assert step[kind] == abc
                    else:
                        step[kind] = abc
                else:
                    target = step.setdefault(kind, [0, 0, 0])
                    for i, v in enumerate(abc):
                        target[i] += v
            assert len(seen) == len(record['selections']['historical'])*4
        rows = []
        for n, step in sorted(steps.items()):
            assert step[4][0] == step[4][1] and step[5][0] == step[5][1]
            rows.append({'sequence': n, 'common': n in record['selections']['common'],
                         'catalog_before_after_changed_keys': step[50],
                         'matching_keys_changed_owners_changed_weights': step[51],
                         'sent_blocks': step[4][0], 'sent_bytes': step[5][0]})
        repeats.append(rows)
    assert repeats[0] == repeats[1]
    selected = [x for x in repeats[0] if x['common']]
    return {'steps': repeats[0], 'repeats_identical': True,
            'common_mean_sent_blocks': st.mean(x['sent_blocks'] for x in selected),
            'common_mean_sent_bytes': st.mean(x['sent_bytes'] for x in selected),
            'scope': 'Aggregate send totals, not off-node traffic or savings; catalog equality does not establish full topology/value equality.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-root', type=Path, default=Path('/private/tmp'))
    parser.add_argument('--output', type=Path, default=Path('docs/parallel-blocks-evidence-reconciliation-210.json'))
    args = parser.parse_args()
    result, steps = {'campaigns': {}}, []
    for stage, job in JOBS.items():
        root = args.archive_root/f'stage{stage}-review-{job}{"-complete" if stage == 209 else ""}'
        q, rows = campaign(root, stage); result['campaigns'][str(stage)] = q; steps.extend(rows)
        print(f'PASS Stage{stage}: {q["integrity"]["verified_files"]} hashes, {len(q["runs"])} raw analyses', flush=True)
    campaigns = list(result['campaigns'].values())
    assert all(q['baseline_sources'] == campaigns[0]['baseline_sources'] for q in campaigns)
    assert all(q['inputs'] == campaigns[0]['inputs'] for q in campaigns)
    assert len({tuple(q['runs'][0]['states']) for q in campaigns}) == 1
    result['historical_report_only'] = {}
    for stage, job in [(200,54689),(201,54691),(202,54693),(203,54695),(204,54697),(205,54699)]:
        path = Path(f'docs/parallel-blocks-cluster-{stage}-{job}.json')
        q = read(path)
        result['historical_report_only'][str(stage)] = {'path': str(path), 'sha256': digest(path),
            'bare': q.get('bare', q.get('summary', {}).get('bare')),
            'scope': 'Prior reviewed record; underlying raw archive not available in the current extracted inputs. No common-step reanalysis claimed.'}
        if stage == 200:
            result['historical_report_only'][str(stage)]['cases'] = {name: {'bare': value['bare'], 'gap_seconds': value['gap_seconds']} for name, value in q['cases'].items()}
    result['same_baseline_source_inventory'] = True
    result['same_inputs_and_printed_trajectory'] = True
    result['main_commit'] = subprocess.check_output(['git', 'rev-parse', 'main'], text=True).strip()
    assert result['main_commit'] == 'b82647b36eb825aed0c0aae4e92adcf4ab46d608'
    result['scope'] = 'Raw Stage206–210 reconciliation; historical reports remain frozen. Rounded stdout is rank0 wall time, not CPU time, rank maximum or synchronized makespan. MPI-owned rows are additive; inclusive fields are not. Statistical bridges are descriptive identities, not causal overhead estimates.'
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    with args.output.with_suffix('.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(steps[0])); writer.writeheader(); writer.writerows(steps)
    print(args.output)


if __name__ == '__main__':
    main()
