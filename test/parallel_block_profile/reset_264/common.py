"""Stage264: sealed Stage230/main baseline with no prior campaign dependencies."""
import hashlib
import json
import os
from pathlib import Path
import statistics

STAGE230 = '4ea64012'  # Expanded to the full commit in IDENTITY.json when sealed.
MAIN = '6a245db7f1632b0205c2bebe95e5f0cd535c2582'
J5_HASHES = {
    'simple.in': 'ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
    'test_checkpoint_0003.bin.zst': '8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d',
    'grids/HR95JT_WT_004': 'd4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d',
}
CHECKPOINTS = {
    'main': '5a60b279f97aa48098126201e4e5a549f759a04c8a9d0ef18e94648ccba49765',
    'stage230': '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53',
}
ORDER = ('main', 'stage230', 'stage230', 'main') * 2
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
    if not identity['stage230_commit'].startswith(STAGE230) or identity['main_commit'] != MAIN:
        raise ValueError('Wrong reset baseline source identities')
    for name, expected in J5_HASHES.items():
        if sha(fixture / name) != expected:
            raise ValueError('Established J5 fixture differs: ' + str(fixture / name))
    # Exercise the complete analysis import graph before submitting a job.
    from capability_223 import coverage, numerical
    from protocol_229 import ordinary
    from run_local_230 import vtk_members
    from placement_230 import plan
    print('PASS: self-contained Stage230/main sources, scripts and established J5 fixture', flush=True)
    return identity


def runtime(oracles=False):
    from profile_protocol import environment
    from run_local_220 import SELECTORS
    env = environment(oracles=oracles, detail=False)
    env.update(SELECTORS)
    env.update(MV2_ENABLE_AFFINITY='0', OMPI_MCA_hwloc_base_binding_policy='none',
               WAVETRISK_SAMPLE_RANKS201='[]')
    return env


def summarize(rows):
    if [r['kind'] for r in rows] != list(ORDER):
        raise ValueError('Require two independent main/Stage230 ABBA batches')
    means = {kind: statistics.mean(r['ordinary']['median'] for r in rows if r['kind'] == kind)
             for kind in ('main', 'stage230')}
    pairs = []
    for start in (0, 4):
        batch = rows[start:start + 4]
        pairs += [batch[1]['ordinary']['median'] - batch[0]['ordinary']['median'],
                  batch[2]['ordinary']['median'] - batch[3]['ordinary']['median']]
    return dict(means_seconds=means, stage230_minus_main_seconds=means['stage230'] - means['main'],
                stage230_to_main_ratio=means['stage230'] / means['main'],
                paired_gaps_seconds=pairs, independent_runs_per_version=4,
                decision='baseline_measured_no_promotion',
                scope='Mean of independent ordinary-step run medians; validation runs excluded.')
