"""Stage273: balance contiguous block ownership; uninstrumented scaling against Stage266/main."""
import hashlib
import json
import os
from pathlib import Path
import statistics

BASELINE = 'ff740744608731a05d53771c32a3055309f3a586'
MAIN = '6a245db7f1632b0205c2bebe95e5f0cd535c2582'
J5_HASHES = {
    'simple.in': 'ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
    'test_checkpoint_0003.bin.zst': '8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d',
    'grids/HR95JT_WT_004': 'd4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d',
}
CHECKPOINTS = {
    'main': '5a60b279f97aa48098126201e4e5a549f759a04c8a9d0ef18e94648ccba49765',
    'stage266': '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53',
}
TIMING_RANKS = (83, 160, 256)
ORDER = tuple((n, kind) for n in TIMING_RANKS for kind in
              (('main',) if n <= 160 else ()) + ('stage266', 'candidate', 'candidate', 'stage266'))
VALIDATIONS = ((16,'main',False),(16,'stage266',False),(16,'candidate',False),
               (16,'checked',False),(16,'candidate',True),(256,'candidate',True))
PLACEMENT_RANKS = (16,) + TIMING_RANKS
STEP_RESTRICTIONS = ('SLURM_CPUS_PER_TASK', 'SRUN_CPUS_PER_TASK', 'SLURM_EXACT',
                     'SLURM_EXCLUSIVE', 'SLURM_TRES_PER_TASK', 'SLURM_HINT',
                     'SLURM_THREADS_PER_CORE')


def sha(path):
    path = Path(path)
    if path.is_symlink():
        return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def verify(root, filename='files-sha256.json'):
    root = Path(root)
    manifest = json.loads((root / filename).read_text())
    for name, expected in manifest.items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or sha(root / name) != expected:
            raise ValueError('Sealed file differs: ' + name)
    return manifest


def preflight(kit, fixture):
    kit, fixture = Path(kit).resolve(strict=True), Path(fixture).resolve(strict=True)
    verify(kit)
    identity = json.loads((kit / 'IDENTITY.json').read_text())
    if identity['stage'] != 273 or identity.get('revision') != 'balanced-partition' or identity['stage266_commit'] != BASELINE or identity['main_commit'] != MAIN:
        raise ValueError('Wrong reset baseline source identities')
    for name, expected in J5_HASHES.items():
        if sha(fixture / name) != expected:
            raise ValueError('Established J5 fixture differs: ' + str(fixture / name))
    # Exercise the complete analysis import graph before submitting a job.
    from capability_223 import coverage, numerical
    from protocol_229 import ordinary
    from run_local_230 import vtk_members
    from placement_273 import plan
    source_scope(kit)
    print('PASS: sealed Stage266/main/candidate sources, bounded source delta and established J5 fixture', flush=True)
    return identity


def runtime(oracles=False, catalogue=False):
    from profile_protocol import environment
    from run_local_220 import SELECTORS
    env = environment(oracles=oracles, detail=False)
    env.update(SELECTORS)
    env.update(MV2_ENABLE_AFFINITY='0', OMPI_MCA_hwloc_base_binding_policy='none',
               WAVETRISK_SAMPLE_RANKS201='[]', WAVETRISK_COST_AUDIT='0',
               WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS='1' if catalogue else '0')
    return env


def summarize(rows):
    if [(r['ranks'],r['kind']) for r in rows] != list(ORDER):
        raise ValueError('Require complete fixed 83/160/256 scaling schedule')
    counts={}
    for n in TIMING_RANKS:
        batch=[r for r in rows if r['ranks']==n]
        groups={k:[r['ordinary']['median'] for r in batch if r['kind']==k] for k in ('stage266','candidate')}
        means={k:statistics.mean(v) for k,v in groups.items()}
        pairs=[groups['candidate'][i]-groups['stage266'][i] for i in range(2)]
        means['main']=next((r['ordinary']['median'] for r in batch if r['kind']=='main'),None)
        gap=means['candidate']-means['stage266']
        counts[str(n)]=dict(means_seconds=means,candidate_minus_stage266_seconds=gap,
            paired_candidate_gaps_seconds=pairs,
            candidate_minus_main_seconds=None if means['main'] is None else means['candidate']-means['main'],
            stage266_minus_main_seconds=None if means['main'] is None else means['stage266']-means['main'],
            independent_runs=dict(stage266=2,candidate=2,main=1 if means['main'] is not None else 0),
            spans_seconds={k:max(v)-min(v) for k,v in groups.items()},
            assessment='mixed_pair_signs' if min(pairs)<0<max(pairs) else ('both_pairs_faster' if max(pairs)<0 else 'no_consistent_gain'))
    scaling={}
    for kind in ('stage266','candidate'):
        scaling[kind]={}
        for anchor in (83,160):
            base=counts[str(anchor)]['means_seconds'][kind]
            scaling[kind][str(anchor)]={str(n):dict(speedup=base/counts[str(n)]['means_seconds'][kind],
                relative_efficiency=(base/counts[str(n)]['means_seconds'][kind])*anchor/n)
                for n in TIMING_RANKS if n>=anchor}
    return dict(by_ranks=counts,scaling=scaling,geometric_domains=160,strategic_target_ranks_minimum=640,
                automatic_promotion=False,scope='Two independent uninstrumented candidate/baseline runs per count; one main anchor at 83/160. Ordinary steps 2-7 and 10-12. No main beyond its supported Domain limit. No 640-rank claim.')


def source_scope(kit):
    manifest = json.loads((Path(kit) / 'source-files.json').read_text())
    if set(manifest) != {'stage266','main','candidate'}:
        raise ValueError('Require exactly three source families')
    a, b = manifest['stage266'], manifest['candidate']
    changed = sorted(name for name in set(a) | set(b) if a.get(name) != b.get(name))
    expected = ['src/parallel_block_mpi.f90']
    identity = json.loads((Path(kit) / 'IDENTITY.json').read_text())
    if changed != expected or identity['changed_source'] != expected:
        raise ValueError('Candidate source exceeds the one-file block partition change')
    if any(b[name] != identity['candidate_source_sha256'][name] for name in expected):
        raise ValueError('Candidate source identity differs')
    if any(a[name] != identity['baseline_source_sha256'][name] for name in expected):
        raise ValueError('Accepted baseline source identity differs')
    return changed
