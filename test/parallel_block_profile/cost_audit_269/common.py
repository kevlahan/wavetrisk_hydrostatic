"""Stage269 immutable controls, bounded instrumentation and explicit calibration."""
import hashlib,json,os,statistics
from pathlib import Path
BASELINE='ff740744608731a05d53771c32a3055309f3a586'
MAIN='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
J5_HASHES={'simple.in':'ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d','test_checkpoint_0003.bin.zst':'8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d','grids/HR95JT_WT_004':'d4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d'}
CHECKPOINTS={'main':'5a60b279f97aa48098126201e4e5a549f759a04c8a9d0ef18e94648ccba49765','stage266':'1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'}
ORDER=('main','main_observer','main_observer','main','stage266','observer','observer','stage266')
VALIDATIONS=((16,'main',False),(16,'stage266',False),(16,'main_observer',False),(16,'observer',False),(16,'checked',False),(16,'observer',True),(167,'observer',True))
STEP_RESTRICTIONS=('SLURM_CPUS_PER_TASK','SRUN_CPUS_PER_TASK','SLURM_EXACT','SLURM_EXCLUSIVE','SLURM_TRES_PER_TASK','SLURM_HINT','SLURM_THREADS_PER_CORE')

def sha(path):
    path=Path(path)
    if path.is_symlink():return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def save(path,value):
    path=Path(path);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def verify(root,filename='files-sha256.json'):
    manifest=json.loads((Path(root)/filename).read_text())
    for name,expected in manifest.items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha(Path(root)/name)!=expected:raise ValueError('Sealed file differs: '+name)
    return manifest

def source_scope(kit):
    kit=Path(kit);m=json.loads((kit/'source-files.json').read_text());identity=json.loads((kit/'IDENTITY.json').read_text())
    if set(m)!={'stage266','main','observer','main_observer'}:raise ValueError('Require four source families')
    for observer,control in (('observer','stage266'),('main_observer','main')):
        changed=sorted(name for name in set(m[observer])|set(m[control]) if m[observer].get(name)!=m[control].get(name))
        if changed!=identity['changed_source'][observer]:raise ValueError('Unsealed diagnostic overlay')
        if any(not(name=='Makefile' or name.startswith('src/') or name=='test/climate/climate.f90') for name in changed):raise ValueError('Unexpected observer scope')
    return identity['changed_source']

def preflight(kit,fixture):
    kit,fixture=map(lambda p:Path(p).resolve(strict=True),(kit,fixture));verify(kit)
    identity=json.loads((kit/'IDENTITY.json').read_text())
    if (identity['stage'],identity['revision'],identity['stage266_commit'],identity['main_commit'])!=(269,'bounded-cost-audit',BASELINE,MAIN):raise ValueError('Wrong campaign identity')
    for name,expected in J5_HASHES.items():
        if sha(fixture/name)!=expected:raise ValueError('Established fixture differs: '+name)
    from capability_223 import coverage,numerical
    from protocol_229 import ordinary
    from run_local_230 import vtk_members
    from placement_230 import plan
    from cost_analysis import read_costs,calibrate
    source_scope(kit)
    print('PASS: sealed untouched Stage266/main controls, separate observer overlays and established J5 fixture',flush=True)
    return identity

def runtime(oracles=False,observer=False):
    from profile_protocol import environment
    from run_local_220 import SELECTORS
    env=environment(oracles=oracles,detail=False);env.update(SELECTORS)
    env.update(MV2_ENABLE_AFFINITY='0',OMPI_MCA_hwloc_base_binding_policy='none',WAVETRISK_SAMPLE_RANKS201='[]',WAVETRISK_COST_AUDIT='1' if observer else '0')
    return env

def summarize(rows):
    if [r['kind'] for r in rows]!=list(ORDER):raise ValueError('Require two bare-observer-observer-bare calibration batches')
    from cost_analysis import calibrate,compare_costs
    calibration={family:calibrate(rows[start:start+4]) for family,start in (('main',0),('stage266',4))}
    means={family:calibration[family]['control_seconds'] for family in calibration}
    return dict(means_seconds=means,candidate_minus_main_seconds=means['stage266']-means['main'],
                calibration=calibration,attribution=compare_costs(rows,calibration),automatic_promotion=False,
                limitation='Diagnostic attribution only; no intervention tested and no recoverable speedup established.')
