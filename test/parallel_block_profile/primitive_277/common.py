"""Stage277: cached fused production kernel and uninstrumented legacy/scaling controls."""
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
ORDER = ((83,'main'),(83,'stage266'),(83,'candidate'),(160,'main'),(160,'stage266'),(160,'candidate'),(256,'stage266'),(256,'candidate'))
VALIDATIONS = ((16,'main',False),(16,'stage266',False),(16,'candidate',False),
               (16,'checked',False),(16,'candidate',True),
               (83,'candidate',False),(160,'candidate',False),(256,'candidate',True))
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


def verify_packaged_grid(kit):
    """Bind the fixture checksum to actual sealed grid bytes before delivery."""
    name='grids/HR95JT_WT_004'
    actual=sha(Path(kit)/name)
    if actual!=J5_HASHES[name]:
        raise ValueError(f'Packaged J5 grid checksum differs: expected {J5_HASHES[name]}, observed {actual}')
    return actual


def preflight(kit, fixture):
    kit, fixture = Path(kit).resolve(strict=True), Path(fixture).resolve(strict=True)
    verify(kit)
    verify_packaged_grid(kit)
    identity = json.loads((kit / 'IDENTITY.json').read_text())
    if identity['stage'] != 277 or identity.get('revision') != 'primitive-production' or identity['stage266_commit'] != BASELINE or identity['main_commit'] != MAIN:
        raise ValueError('Wrong reset baseline source identities')
    for name, expected in J5_HASHES.items():
        actual=sha(fixture/name)
        if actual != expected:
            raise ValueError(f'Established J5 fixture differs: {fixture/name}; expected {expected}, observed {actual}')
    # Exercise the complete analysis import graph before submitting a job.
    from capability_223 import coverage, numerical
    from protocol_229 import ordinary
    from run_local_230 import vtk_members
    from placement_277 import plan
    from contract import evidence
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
               WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS='0', WAVETRISK_PRIMITIVE277_CHECK='1' if catalogue else '0')
    return env


def summarize(rows):
    if [(r['ranks'],r['kind']) for r in rows]!=list(ORDER):raise ValueError('Timing schedule differs')
    medians={f"{r['ranks']}-{r['kind']}":r['ordinary']['median'] for r in rows}
    by_rank={}
    for n in TIMING_RANKS:
        baseline=medians[f'{n}-stage266']; candidate=medians[f'{n}-candidate']
        main=medians.get(f'{n}-main')
        by_rank[str(n)]=dict(stage266_seconds=baseline,candidate_seconds=candidate,main_seconds=main,
            candidate_minus_stage266_percent=100*(candidate/baseline-1),
            candidate_minus_main_percent=None if main is None else 100*(candidate/main-1))
    return dict(automatic_promotion=False,control_step_seconds=medians,by_ranks=by_rank,
                speedup_160_to_256={k:medians[f'160-{k}']/medians[f'256-{k}'] for k in ('stage266','candidate')},
                scope='Single uninstrumented run per configuration. No automatic promotion; resolve a repeatable regression before stacking another change. Production primitive kernel integration remains on Domain owners, so this is not a beyond-Domain ownership migration.')


def source_scope(kit):
    manifest = json.loads((Path(kit) / 'source-files.json').read_text())
    if set(manifest) != {'stage266','main','candidate'}:
        raise ValueError('Require exactly three source families')
    a, b = manifest['stage266'], manifest['candidate']
    changed = sorted(name for name in set(a) | set(b) if a.get(name) != b.get(name))
    expected = ['src/multi_level.f90','src/ops.f90']
    identity = json.loads((Path(kit) / 'IDENTITY.json').read_text())
    if changed != expected or identity['changed_source'] != expected:
        raise ValueError('Candidate source exceeds the two-file production overlay')
    if any(b[name] != identity['candidate_source_sha256'][name] for name in expected):
        raise ValueError('Candidate source identity differs')
    if any(a[name] != identity['baseline_source_sha256'][name] for name in expected):
        raise ValueError('Accepted baseline source identity differs')
    return changed
