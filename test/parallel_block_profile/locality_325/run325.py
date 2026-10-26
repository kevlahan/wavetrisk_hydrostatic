"""Bounded init/four-step correctness and changed-work gate. No timing path."""
from pathlib import Path
from datetime import datetime,timezone
import json,os,shutil,socket,subprocess,sys,time,signal
from common import preflight,save,sha
from config325 import SCHEDULE,environment,STEP_RESTRICTIONS
from continuation_binding import verify_binding,recheck_artifacts
from build import build
from placement_278 import topology,plan,verify
from validation322 import (SEED_NAME,private_copy,checkpoint_shape,comparison,snapshots,census,
                           trajectory,log_evidence,saved_environment,activation)
from runtime322 import verify_solver_placement

def utc():return datetime.now(timezone.utc).isoformat()
def require_allocation():
 if not os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOB_PARTITION')!='nk' or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or socket.gethostname().split('.')[0]!='bb13' or os.environ.get('SLURM_NTASKS')!='256':raise ValueError('Require exclusive256-task one-node bb13/nk allocation')
 for key in ('LD_PRELOAD','LD_AUDIT'):
  if os.environ.get(key):raise ValueError('Unsupported loader override: '+key)
 for cmd in ('srun','numactl','lscpu','zstd','mpif90','ldd','readelf'):
  if not shutil.which(cmd):raise ValueError('Missing '+cmd)

def launch(kit,out,*args):
 return ['srun','--nodes=1','--ntasks=256','--cpu-bind=none','--kill-on-bad-exit=1',sys.executable,str(kit/'rank_launch_223.py'),str(out/'placement-256.json'),*map(str,args)]

CAPTURE_RANK_BYTES=64*1024**2
CAPTURE_RUN_BYTES=2*1024**3
CAPTURE_CAMPAIGN_BYTES=3*1024**3

def capture_usage(out,active):
 files=list(Path(out).glob('*/execution316-rank-*.bin.gz'))
 sizes=[]
 for p in files:
  if p.is_symlink():raise ValueError('Capture became a symlink')
  sizes.append((p,p.stat().st_size))
 total=sum(n for _,n in sizes);current=sum(n for p,n in sizes if p.parent==Path(active))
 if any(n>CAPTURE_RANK_BYTES for _,n in sizes) or current>CAPTURE_RUN_BYTES or total>CAPTURE_CAMPAIGN_BYTES:raise ValueError('Live compact capture budget exceeded; stopping solver and retaining partial evidence')
 return dict(campaign_stored_bytes=total,run_stored_bytes=current,files=len(files))

def run_bounded(command,directory,env,log,out,timeout=600):
 """Monitor file sizes without reading streams; SIGTERM propagates through srun.

 Aggregate caps are monitored stop thresholds, not atomic filesystem quotas.
 The native per-rank sink cap is the independent hard bound; files can exceed
 an aggregate threshold between half-second polls and during termination.
 """
 capture_usage(out,directory)
 process=subprocess.Popen(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 started=time.monotonic()
 try:
  while process.poll() is None:
   capture_usage(out,directory)
   if time.monotonic()-started>timeout:raise TimeoutError('Bounded qualification solver exceeded600seconds')
   time.sleep(0.5)
  usage=capture_usage(out,directory)
  return subprocess.CompletedProcess(command,process.returncode),usage
 except BaseException:
  if process.poll() is None:
   try:os.killpg(process.pid,signal.SIGTERM)
   except ProcessLookupError:pass
   try:process.wait(timeout=10)
   except subprocess.TimeoutExpired:
    try:os.killpg(process.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    process.wait(timeout=10)
  raise

def bind_prior_toolchain(current,binding):
 old=binding['original_toolchain']
 for key in ('wrapper_resolved','wrapper_sha256','compiler_resolved','compiler_sha256','compiler_version','c_compiler_resolved','c_compiler_sha256','c_compiler_version','target','wrapper_command','expected','library_version'):
  if current[key]!=old[key]:raise ValueError('Current toolchain differs from bound baseline: '+key)
 return dict(passed=True,fields_checked=12,old_job=55561)

def run(kit,out):
 kit=Path(kit).resolve(strict=True);out=Path(out).resolve()
 identity=preflight(kit);require_allocation()
 if out.exists() or out.is_relative_to(kit) or kit.is_relative_to(out):raise ValueError('Fresh separate results directory required')
 out.mkdir();start=time.monotonic()
 result=dict(stage=325,status='running',identity=identity,job_id=os.environ['SLURM_JOB_ID'],started_utc=utc(),runs={},schedule=[],performance_measurement=False,performance_qualified=False,qualification=None)
 def persist():result['elapsed_seconds']=time.monotonic()-start;save(out/'results.json',result)
 persist()
 try:
  print('Binding prior source, native toolchain, private fresh-restart seed and legacy oracle; no raw capture reads',flush=True)
  binding=verify_binding(kit);save(out/'prior-binding.json',binding);result['prior_binding_sha256']=sha(out/'prior-binding.json');persist()
  prior=Path(binding['original_results']);seed=Path(binding['private_seed'])
  for protected in (prior,Path(binding['original_kit'])):
   if out.is_relative_to(protected) or protected.is_relative_to(out):raise ValueError('Output overlaps preserved baseline')
  result['removed_step_restrictions']={k:os.environ.pop(k) for k in STEP_RESTRICTIONS if k in os.environ}
  result['builds']=build(kit,out);result['toolchain_binding']=bind_prior_toolchain(result['builds']['toolchain'],binding);persist()
  top=topology();placements=plan(top);save(out/'topology.json',top);save(out/'placement-256.json',placements[256])
  prior_place=json.loads((prior/'placement-256.json').read_text())
  if placements[256]!=prior_place:raise ValueError('Native256 placement differs from audited map placement')
  p=subprocess.run(launch(kit,out,'--probe'),text=True,capture_output=True,timeout=90)
  (out/'placement-256-probe.log').write_text(p.stdout+p.stderr)
  if p.returncode:raise RuntimeError('Placement probe failed')
  rows=[json.loads(s) for s in p.stdout.splitlines() if s.startswith('{')];verify(rows,placements[256]);save(out/'placement-256-verified.json',rows)
  old_result=json.loads((prior/'results.json').read_text());legacy=old_result['runs']['legacy-four']
  for ordinal,spec in enumerate(SCHEDULE,1):
   label=spec['label'];kind=spec['kind'];d=out/label;d.mkdir();row=dict(spec,ordinal=ordinal,status='preparing');result['runs'][label]=row;result['schedule'].append(label);result['active']=label;persist()
   shutil.copy2(kit/'simple.in',d/'simple.in');shutil.copytree(kit/'grids',d/'grids')
   if sha(seed)!=binding['private_seed_sha256']:raise ValueError('Bound private seed changed')
   row['seed_sha256']=private_copy(seed,d/SEED_NAME)
   record=result['builds']['builds'][kind];binary=out/record['binary']
   if sha(binary)!=record['binary_sha256']:raise ValueError('Compiled executable changed')
   (d/'climate').symlink_to(binary)
   env=environment(kit,kind,spec['owner']);save(d/'environment.json',saved_environment(env))
   cmd=launch(kit,out,d,'./climate','simple.in');row.update(command=cmd,binary_sha256=sha(binary),input_sha256=sha(d/'simple.in'),environment_sha256=sha(d/'environment.json'),status='running',started_utc=utc(),start_monotonic_ns=time.monotonic_ns());persist()
   print(f'Starting {ordinal}/3 {label}: level5 tol0,256 ranks,4 steps; correctness/work only',flush=True)
   with (d/'run.log').open('x') as log:
    try:p,row['capture_budget']=run_bounded(cmd,d,env,log,out)
    finally:
     row.update(end_monotonic_ns=time.monotonic_ns(),ended_utc=utc());row['solver_launch_elapsed_seconds']=(row['end_monotonic_ns']-row['start_monotonic_ns'])/1e9;persist()
   row['returncode']=p.returncode;persist()
   if p.returncode:raise RuntimeError(label+': solver failed; full log retained')
   if sha(d/SEED_NAME)!=row['seed_sha256'] or sha(seed)!=row['seed_sha256']:raise ValueError('Read-only seed was modified')
   row['activation']=activation(d);row['snapshots']=snapshots(d,256,4);row['census']=raw=census(d,256,4,True)
   resources=[json.loads((d/f'rank-{r}.json').read_text()) for r in range(256)];verify(resources,placements[256])
   if any(r.get('returncode')!=0 for r in resources):raise ValueError('Rank child failure')
   row['solver_placement']=verify_solver_placement(d,256,placements[256],env.get('WAVETRISK_EXECUTION_AUDIT316')=='1',4)
   row['events']=events=log_evidence((d/'run.log').read_text(),4,True,256)
   if trajectory(raw)!=trajectory(legacy['census']) or events!=legacy['events']:raise ValueError('Prior independent legacy trajectory/events differ')
   refs=binding['reference_checkpoints']['legacy-four']
   row['science']={
    'initial-legacy-exact':comparison(Path(refs['initial']),d/'test_checkpoint_9900.bin.zst',True),
    'final-legacy-tolerance':comparison(Path(refs['final']),d/'test_checkpoint_9901.bin.zst',False)}
   if label=='qual-C0':
    ref=binding['reference_checkpoints']['qual-C0'];row['science']['control-prior-exact']=comparison(Path(ref['final']),d/'test_checkpoint_9901.bin.zst',True)
   if label=='qual-L256':
    for index in ('9900','9901'):row['science']['paired-'+index+'-exact']=comparison(out/'qual-C0'/f'test_checkpoint_{index}.bin.zst',d/f'test_checkpoint_{index}.bin.zst',True)
   if not all(x['passed'] for x in row['science'].values()):raise ValueError('Science rejection: '+label)
   save(d/'science325.json',row['science'])
   row['checkpoints']={str(i):sha(d/f'test_checkpoint_{i}.bin.zst') for i in (9900,9901)}
   from account325 import validate_run
   row['work_ledger']=validate_run(d,nrank=256,arm='C0' if label=='qual-C0' else 'L256')
   from review_work325 import project
   row['hardware_work_review']=project(d,placements[256])
   row['capture_budget']=capture_usage(out,d)
   if row['work_ledger'].get('passed') is not True:raise ValueError('Per-run work qualification rejected')
   row.update(status='passed',checks_ended_utc=utc());save(d/'run-summary325.json',row);persist()
   print('PASS '+label+'; science/placement/work, semantic verification pending',flush=True)
  from verify325 import qualify
  from account325 import qualify as qualify_work
  semantic=qualify(out,kit,placements);save(out/'semantic-qualification325.json',semantic)
  work=qualify_work(out,kit);save(out/'work-qualification325.json',work)
  if semantic.get('passed') is not True or work.get('passed') is not True:raise ValueError('Complete semantic/work gate rejected')
  result['prior_preservation']=recheck_artifacts(binding)
  if result['schedule']!=[s['label'] for s in SCHEDULE] or any(r['status']!='passed' for r in result['runs'].values()):raise ValueError('Incomplete sealed schedule')
  result.update(status='qualification_complete_requires_independent_review',decision='HOLD_PERFORMANCE_PENDING_INDEPENDENT_QUALIFICATION_REVIEW',qualification=dict(passed=True,semantic_sha256=sha(out/'semantic-qualification325.json'),work_sha256=sha(out/'work-qualification325.json')),ended_utc=utc());result.pop('active',None);persist()
  print('COMPLETE Stage325 qualification; no performance runs scheduled or authorized by this result',flush=True)
 except BaseException as exc:
  if result.get('active') in result['runs'] and result['runs'][result['active']]['status']!='passed':result['runs'][result['active']]['status']='failed'
  result.update(status='failed',decision='STOP_RETAIN_PARTIAL_EVIDENCE_NO_PERFORMANCE_INFERENCE',error=repr(exc),ended_utc=utc());persist();raise
if __name__=='__main__':
 if len(sys.argv)!=3:raise SystemExit('Usage: run325.py KIT NEW_RESULTS')
 run(*sys.argv[1:])
