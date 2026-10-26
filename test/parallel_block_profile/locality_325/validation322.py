"""Raw-record and checkpoint gates for the fresh single-level campaign."""
from pathlib import Path
import hashlib,json,math,re,shutil
import numpy as np
from analyze322 import census as raw_census
from checkpoint_gate import accept
from compare_restart_checkpoints import compare,read

SEED_NAME='test_checkpoint_0003.bin.zst'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical_sha(value):return hashlib.sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def private_copy(source,destination):
 from common import private_copy_checkpoint
 result=private_copy_checkpoint(source,destination)
 Path(destination).chmod(0o444)
 return result['sha256']

def checkpoint_shape(path):
 records=read(path,-10,30);headers=[];nonzero=0
 if len(records)!=160:raise ValueError('Require160 Domains')
 for header,threshold,poles,coarse,waves,required in records:
  headers.append(header)
  if any(not np.isfinite(x).all() for x in (threshold,poles,coarse,waves)):raise ValueError('Nonfinite checkpoint')
  # Source test_case_module sets the soil sentinel separately from atmosphere.
  levels=threshold.reshape(41,3)
  if np.any(levels[:11]!=1e16) or np.any(levels[11:]!=0):raise ValueError('Soil sentinel or tol0 atmospheric thresholds differ')
  if required.shape!=(4,4) or np.any(required!=0):raise ValueError('Require exactly4 level5 root patches per Domain and no children')
  nonzero+=int(np.count_nonzero(waves))
 if any(h!=headers[0] for h in headers):raise ValueError('Domain checkpoint time metadata differs')
 h=headers[0]
 if not math.isfinite(h[1]) or h[1]<0:raise ValueError('Invalid checkpoint physical time')
 return dict(domains=160,patches_by_level={'5':640},maximum_level=5,header=list(h),nonzero_wavelet_values=nonzero)

def comparison(left,right,exact):
 # Decode both even for identical compressed bytes: equality does not prove finite valid data.
 a=checkpoint_shape(left);b=checkpoint_shape(right)
 same=sha(left)==sha(right)
 report=None if same else compare(left,right)
 failures=[] if same else accept(report,exact=exact)
 return dict(passed=not failures,sha256_exact=same,exact_fields_required=exact,failures=failures,
             left_shape=a,right_shape=b,report=report)

def snapshots(directory,n,steps):
 root=Path(directory);indices=[9900]+([9901] if steps else [])
 if {p.name for p in root.glob('snapshot312-rank-*.txt')}!={f'snapshot312-rank-{r}.txt' for r in range(n)}:raise ValueError('Snapshot preservation rank extent')
 totals={i:[0,0] for i in indices}
 for r in range(n):
  rows=[list(map(int,line.split())) for line in (root/f'snapshot312-rank-{r}.txt').read_text().splitlines()]
  if len(rows)!=len(indices) or any(len(v)!=3 for v in rows) or [v[0] for v in rows]!=indices or any(x<0 for v in rows for x in v):raise ValueError('Snapshot preservation extent/order')
  for index,a,b in rows:totals[index][0]+=a;totals[index][1]+=b
 return {str(i):dict(restored_deletion_flags=x[0],restored_level_entries_or_lengths=x[1]) for i,x in totals.items()}

def stage_work(directory,n,steps):
 result=[]
 for r in range(n):
  rows=[];previous=None;old=[0]*5
  for line in (Path(directory)/f'campaign312-rank-{r}.txt').read_text().splitlines():
   f=line.split()
   if not f or f[0]!='R':continue
   if len(f)!=13:raise ValueError('Malformed RHS record')
   step,gen,mask,nd,groups,prim,scalar,vel,rv,rm,nnode,call=map(int,f[1:]);cur=[prim,scalar,vel,rv,rm];key=(gen,mask)
   delta=[a-b for a,b in zip(cur,old)] if key==previous else cur
   if min(delta)<0:raise ValueError('Negative work delta')
   rows.append([step,nd,groups,nnode,*delta]);previous=key;old=cur
  if [x[0] for x in rows]!=[s for s in range(1,steps+1) for _ in range(4)]:raise ValueError('Incomplete per-RK work')
  result.append(rows)
 return result

def physical_maps(directory):
 return {phase:[list(map(int,f[1:])) for line in (Path(directory)/f'groups312-{phase}.txt').read_text().splitlines() if (f:=line.split()) and f[0]=='E'] for phase in ('initial','final')}

def census(directory,n,steps,block):
 for rank in range(n):
  first=(Path(directory)/f'campaign312-rank-{rank}.txt').read_text().splitlines()[0].split()
  if len(first)!=5 or not math.isfinite(float(first[4])) or float(first[4])<=0:raise ValueError('Invalid per-rank initialization duration')
 value=raw_census(directory,n,steps,block)
 if not math.isfinite(value['initialization_max_seconds']) or value['initialization_max_seconds']<=0:raise ValueError('Invalid initialization duration')
 if block:
  for phase,g in value['groups'].items():
   if g['original_roots']!=640 or g['execution_groups']!=640 or g['fine_units']!=640:raise ValueError('Single-level640-group extent differs')
 return value

def trajectory(raw):return [{k:s[k] for k in ('physical_time','dt','active_wavelets')} for s in raw['steps']]

def require_bare(directory,env,log):
 for key,value in env.items():
  if key.startswith('WAVETRISK_') and any(x in key for x in ('PROFILE','VALIDATE','AUDIT','_CHECK','DIAGNOSTICS')) and value!='0':raise ValueError('Timing observer enabled: '+key)
 if env.get('WAVETRISK_SAMPLE_RANKS201')!='[]':raise ValueError('Sampler enabled')
 for pattern in ('execution316-*','semantic322-*','boundary310-rank-*','rhs278-*','cost-rank-*','trace314-*','resident315-*','balance-rank-*'):
  if any(Path(directory).glob(pattern)):raise ValueError('Unexpected observer output: '+pattern)
 if 'validation oracle enabled' in log.lower():raise ValueError('Oracle enabled in bare run')

def log_evidence(log,steps,block,ranks,bootstrap=False):
 if len(re.findall(r'^\d{8}\s+.* cpu =',log,re.M))!=steps:raise ValueError('Wrong completed step count')
 if block!=('Stage278 connected resident RHS = on' in log):raise ValueError('Connected/legacy path differs')
 if ranks>160 and f'MPI tasks = {ranks}, geometric Domains = 160' not in log:raise ValueError('Above-Domain path absent')
 if not bootstrap and not re.search(r'Restarting from checkpoint\s+3\b',log):raise ValueError('Private checkpoint3 restart absent')
 parts=log.split('Start simulation run')
 if len(parts)!=2:raise ValueError('Missing or repeated simulation boundary')
 step=1;remaps=[];outputs=[];triggers=[]
 for line in parts[1].splitlines():
  if 'Remapping vertical coordinates' in line:remaps.append(step)
  if 'Saving checkpoint ' in line:outputs.append(step)
  if 'Vertical-remap trigger =' in line:triggers.append(dict(step=step,reason=line.split('=',1)[1].strip()))
  if re.match(r'^\d{8}\s+.* cpu =',line):step+=1
 if any(x>steps for x in remaps+outputs):raise ValueError('Event outside timed trajectory')
 return dict(remap_steps=remaps,ordinary_checkpoint_steps=outputs,initial_active_wavelets=re.findall(r'time \[d\].*?dof =\s+(\d+)',log))

def work_anchor(directory,n,steps,raw,events):
 return dict(stage_work_sha256=canonical_sha(stage_work(directory,n,steps)),
             physical_maps_sha256={k:canonical_sha(v) for k,v in physical_maps(directory).items()},
             groups=raw['groups'],trajectory=trajectory(raw),events=events)

def saved_environment(env):
 return {k:v for k,v in env.items() if k.startswith(('WAVETRISK_','MV2_','OMPI_','GFORTRAN_')) or k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'}

def activation(directory):
 rows=[line.split() for line in (Path(directory)/'tolerance312.txt').read_text().splitlines()]
 if len(rows)!=2 or any(len(x)!=6 for x in rows):raise ValueError('Malformed tolerance activation record')
 tol=float(rows[0][0]);iteration,missing,low,high,maximum=map(int,rows[0][1:]);bounds=list(map(float,rows[1]))
 if tol!=0 or iteration!=2 or missing!=0 or (low,high,maximum)!=(5,5,5) or any(v!=0 for v in bounds):raise ValueError('Runtime level5/tol0 activation differs')
 return dict(tol=tol,iterations=iteration,missing=missing,level_start=low,level_end=high,max_level=maximum,atmospheric_threshold_bounds=bounds)
