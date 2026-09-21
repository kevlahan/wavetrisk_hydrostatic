"""Correctness gates for the bounded Stage223 task-count capability campaign."""
import hashlib
import json
from pathlib import Path
import re
from analyze import analyze, STEP
from checkpoint_gate import accept, POLICY
from compare_restart_checkpoints import compare


def digest(path):
    if path.is_symlink():
        return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def verify_sources(root, sources):
    for name, expected in sources.items():
        if digest(root/name) != expected:
            raise ValueError('Source differs: '+name)


def verify_kit(kit):
    verify_sources(kit, json.loads((kit/'helpers-sha256.json').read_text()))


def schedule():
    for level, base, extra in (('J4',40,41), ('J5',160,167)):
        yield (level, base, 'baseline', False, 'native')
        yield (level, base, 'candidate', False, 'native')
        yield (level, extra, 'checked', False, 'native')
        yield (level, extra, 'candidate', False, 'native')
        yield (level, extra, 'candidate', True, 'native')
        if level == 'J4':
            yield (level, extra, 'checked', False, 'empty-root')


def coverage(root, level, run_id):
    domains, minimum, seed, final = (40,4,2,3) if level == 'J4' else (160,5,3,4)
    log = (root/'run.log').read_text()
    a = analyze(root/'run.log')
    if not a['completed'] or not a['step_records']:
        raise ValueError('Incomplete solver run')
    if not re.search(r'min_level\s*=\s*'+str(minimum)+r'\b', log):
        raise ValueError('Wrong compiled grid level')
    if not re.search(r'DOMAIN_LEVEL\s*=\s*'+str(minimum-3)+r'\b', log):
        raise ValueError('Wrong compiled Domain count')
    wanted = {f'{run_id}_checkpoint_{i:04d}.bin.zst' for i in (seed, final)}
    if {p.name for p in root.glob('*_checkpoint_*.bin.zst')} != wanted:
        raise ValueError('Unexpected checkpoint set')
    reloads = list(map(int, re.findall(r'Restarting from checkpoint\s+(\d+)', log)))
    if reloads != [seed, final]:
        raise ValueError('Missing checkpoint reload')
    after = log[re.search(r'Restarting from checkpoint\s+'+str(final)+r'\b', log).end():]
    if 'Remapping vertical coordinates' not in after or not any(STEP.match(s) for s in after.splitlines()):
        raise ValueError('Require evolution and vertical remap after reload')
    return a, dict(domains=domains, reloads=reloads, remaps=log.count('Remapping vertical coordinates'),
                   new_checkpoint=f'{run_id}_checkpoint_{final:04d}.bin.zst')


def events(root, ranks, domains, steps, layout):
    paths = sorted(root.glob('capability223-rank-*.txt'))
    if {p.name for p in paths} != {f'capability223-rank-{r}.txt' for r in range(ranks)}:
        raise ValueError('Incomplete rank capability records')
    rows = {}
    for rank in range(ranks):
        parsed = []
        for line in (root/f'capability223-rank-{rank}.txt').read_text().splitlines():
            f = line.split()
            if len(f) != 9 or f[0] not in ('restart','install','dynamics','step'):
                raise ValueError('Malformed capability event')
            seq, r, nd, nb, owner, iw, it = map(int, f[1:8])
            if seq != len(parsed)+1 or r != rank or nd < 0 or not 0 <= owner < ranks:
                raise ValueError('Invalid rank event sequence/ownership')
            if (nb != -1 if f[0] == 'restart' else nb < 1):
                raise ValueError('Rank did not own useful blocks')
            parsed.append(dict(phase=f[0], domains=nd, blocks=nb, source=owner,
                               output=iw, clock=it, time=f[8]))
        if [r['phase'] for r in parsed].count('step') != steps or [r['phase'] for r in parsed].count('dynamics') != steps:
            raise ValueError('Missing completed block dynamics/step evidence')
        if [r['phase'] for r in parsed].count('restart') != 2:
            raise ValueError('Missing rank restart records')
        rows[rank] = parsed
    first = rows[0]
    if any(len(v) != len(first) for v in rows.values()):
        raise ValueError('Rank lifecycle lengths differ')
    zero_sets = []
    for i, reference in enumerate(first):
        group = [rows[r][i] for r in range(ranks)]
        for row in group:
            if any(row[k] != reference[k] for k in ('phase','source','output','clock','time')):
                raise ValueError('Rank global metadata/lifecycle diverged')
        if sum(row['domains'] for row in group) != domains:
            raise ValueError('Domain ownership does not cover the mesh')
        zero = [r for r, row in enumerate(group) if row['domains'] == 0]
        if ranks > domains and len(zero) < ranks-domains:
            raise ValueError('Missing expected empty-Domain ranks')
        if layout == 'empty-root' and (0 not in zero or reference['source'] == 0):
            raise ValueError('Empty-root/nonzero restart source not exercised')
        zero_sets.append(zero)
    block_counts = [row['blocks'] for group in rows.values() for row in group if row['phase'] != 'restart']
    return dict(ranks=ranks, steps_per_rank=steps, dynamics_per_rank=steps, reloads_per_rank=2,
                events_per_rank=len(first), zero_domain_rank_sets=zero_sets,
                owned_blocks_min=min(block_counts), owned_blocks_max=max(block_counts),
                restart_source_ranks=sorted({r['source'] for r in first if r['phase']=='restart'}),
                files={p.name:digest(p) for p in paths})


def numerical(left, right, filename, domains, exact):
    report = compare(left/filename, right/filename, expected_domains=domains)
    failures = accept(report, exact=exact)
    return dict(passed=not failures, exact_required=exact, exact_observed=not accept(report,exact=True),
                failures=failures, policy=POLICY, report=report)
