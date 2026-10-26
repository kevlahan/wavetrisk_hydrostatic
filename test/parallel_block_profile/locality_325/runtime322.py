"""Strict actual-solver placement and shared Linux monotonic window checks."""
from pathlib import Path
import json
import re


def integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('runtime322 integer: ' + label)
    return value


def records(directory, n):
    directory = Path(directory)
    names = {f'runtime322-rank-{r}.json' for r in range(n)}
    if {p.name for p in directory.glob('runtime322-rank-*.json')} != names:
        raise ValueError('runtime322 rank coverage')
    rows = []
    expected_keys = {'rank','pid','host','clock_namespace','time_namespace_for_children',
                     'boot_id','begin_ns','end_ns','observations'}
    for r in range(n):
        value = json.loads((directory / f'runtime322-rank-{r}.json').read_text())
        if set(value) != expected_keys or integer(value['rank'], 'rank') != r:
            raise ValueError('runtime322 record schema/rank')
        integer(value['pid'], 'pid', 1)
        integer(value['begin_ns'], 'begin_ns')
        integer(value['end_ns'], 'end_ns')
        if not isinstance(value['host'], str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', value['host']):
            raise ValueError('runtime322 host')
        for key in ('clock_namespace','time_namespace_for_children'):
            if not isinstance(value[key], str) or not re.fullmatch(r'time:\[\d+\]', value[key]):
                raise ValueError('runtime322 clock namespace format')
        if not isinstance(value['boot_id'], str) or not re.fullmatch(r'[0-9a-f-]{36}', value['boot_id']):
            raise ValueError('runtime322 boot identity')
        rows.append(value)
    for key in ('host','clock_namespace','time_namespace_for_children','boot_id'):
        if len({v[key] for v in rows}) != 1:
            raise ValueError('runtime322 common clock/host: ' + key)
    if len({v['pid'] for v in rows}) != n:
        raise ValueError('runtime322 duplicate solver pid')
    return rows


def verify_solver_placement(directory, n, placement, qualification, steps):
    if steps != 4 or n not in (160, 256):
        raise ValueError('runtime322 configured extent')
    if isinstance(placement, (str, Path)):
        placement = json.loads(Path(placement).read_text())
    mapping = placement['mapping']
    if placement['ranks'] != n or len(mapping) != n:
        raise ValueError('runtime322 planned extent')
    if len({m['cpu'] for m in mapping}) != n or len({(m['socket'],m['core']) for m in mapping}) != n:
        raise ValueError('runtime322 oversubscription')
    rows = records(directory, n)
    expected_stages = list(range(7)) if qualification else [0,1,6]
    helpers = []
    for r, row in enumerate(rows):
        planned = mapping[r]
        if planned['rank'] != r or row['host'] != placement['host']:
            raise ValueError('runtime322 planned host/rank')
        observations = row['observations']
        if not isinstance(observations, list) or [o.get('stage') for o in observations] != expected_stages:
            raise ValueError('runtime322 boundary coverage')
        for observation in observations:
            if set(observation) != {'stage','main_tid','cpu','memory_policy','threads'}:
                raise ValueError('runtime322 boundary schema')
            if integer(observation['main_tid'], 'main_tid', 1) != row['pid']:
                raise ValueError('runtime322 wrong main thread')
            if integer(observation['cpu'], 'cpu') != planned['cpu']:
                raise ValueError('runtime322 current CPU')
            if integer(observation['memory_policy'], 'memory policy') != 4:
                raise ValueError('runtime322 require MPOL_LOCAL')
            threads = observation['threads']
            if not isinstance(threads, list) or not 1 <= len(threads) <= 64:
                raise ValueError('runtime322 thread coverage')
            tids = []
            main = []
            for thread in threads:
                if set(thread) != {'tid','mask'}:
                    raise ValueError('runtime322 thread schema')
                tids.append(integer(thread['tid'], 'tid', 1))
                mask = thread['mask']
                if not isinstance(mask, list) or not mask or mask != sorted(set(mask)):
                    raise ValueError('runtime322 CPU mask')
                if any(integer(c, 'thread CPU') >= 1024 for c in mask):
                    raise ValueError('runtime322 CPU bound')
                if thread['tid'] == row['pid']:
                    main.append(mask)
            if len(tids) != len(set(tids)) or main != [[planned['cpu']]]:
                raise ValueError('runtime322 main affinity')
            helpers.append(dict(rank=r,stage=observation['stage'],helper_threads=len(threads)-1))
        if qualification and (row['begin_ns'] != 0 or row['end_ns'] != 0):
            raise ValueError('runtime322 qualification timing enabled')
    return dict(passed=True,observed_ranks=n,stages=expected_stages,memory_policy='MPOL_LOCAL',
                host=rows[0]['host'],clock_namespace=rows[0]['clock_namespace'],
                helper_threads=helpers,scope='Boundary observations; transient affinity changes are not monitored.')


def read_window(directory, n):
    rows = records(directory, n)
    starts = [v['begin_ns'] for v in rows]
    ends = [v['end_ns'] for v in rows]
    if any(a <= 0 or b <= a for a,b in zip(starts,ends)):
        raise ValueError('runtime322 missing/reversed clock window')
    if min(ends) <= max(starts):
        raise ValueError('runtime322 windows have no common overlap')
    return dict(common_window_seconds=(max(ends)-min(starts))/(3*1e9),
                begin_spread_seconds=(max(starts)-min(starts))/1e9,
                end_spread_seconds=(max(ends)-min(ends))/1e9,
                total_window_seconds=(max(ends)-min(starts))/1e9,
                first_begin_rank=starts.index(min(starts)),last_end_rank=ends.index(max(ends)),
                host=rows[0]['host'],clock_namespace=rows[0]['clock_namespace'],
                boot_id=rows[0]['boot_id'],steps=[2,3,4],
                scope='Common-node elapsed window including existing interstep logging; not whole-job elapsed.')
