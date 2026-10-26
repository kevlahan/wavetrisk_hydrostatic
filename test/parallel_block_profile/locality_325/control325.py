"""One user-launched Stage325 allocation; failures remain collectable, no resubmit."""
from pathlib import Path
import getpass,hashlib,json,os,shutil,subprocess,sys,tarfile,tempfile
from common import verify,sha,save,preflight

def queued(job):
 if isinstance(job,bool) or not str(job).isdigit() or int(job)<=0:raise ValueError('Invalid scheduler job identity')
 q=subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%A'],capture_output=True,text=True)
 if q.returncode==0 and str(job) in q.stdout.split():return True
 a=subprocess.run(['sacct','-n','-P','-X','-j',str(job),'-o','JobIDRaw,State'],capture_output=True,text=True)
 states=[]
 for line in a.stdout.splitlines():
  f=line.strip().split('|')
  if len(f)>1 and f[0]==str(job) and f[1].split():states.append(f[1].split()[0].rstrip('+'))
 terminal={'BOOT_FAIL','CANCELLED','COMPLETED','DEADLINE','FAILED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','REVOKED','SPECIAL_EXIT','TIMEOUT'}
 active={'PENDING','RUNNING','CONFIGURING','COMPLETING','SUSPENDED','RESIZING','REQUEUED','REQUEUE_FED','REQUEUE_HOLD','RESV_DEL_HOLD','SIGNALING','STAGE_OUT','STOPPED'}
 if a.returncode==0 and states and all(s in terminal for s in states):return False
 if a.returncode==0 and any(s in active for s in states):return True
 raise RuntimeError('Job state unknown; collection stopped. squeue: '+q.stderr+' sacct: '+a.stderr+' '+a.stdout)

def verify_archive(path):
 with tarfile.open(path) as t:
  members=t.getmembers();names=[m.name for m in members]
  if len(names)!=len(set(names)):raise ValueError('Duplicate review archive member')
  if any(m.issym() or m.islnk() or m.name.startswith('/') or '..' in Path(m.name).parts or not(m.isfile() or m.isdir()) for m in members):raise ValueError('Unsafe review archive member')
  m=json.load(t.extractfile('review/files-sha256.json'))
  if not isinstance(m,dict):raise ValueError('Invalid review manifest')
  files={v.name for v in members if v.isfile()}
  if files!={'review/files-sha256.json',*('review/'+n for n in m)}:raise ValueError('Review archive inventory differs')
  for n,h in m.items():
   if Path(n).is_absolute() or '..' in Path(n).parts:raise ValueError('Unsafe review inventory path')
   if hashlib.file_digest(t.extractfile('review/'+n),'sha256').hexdigest()!=h:raise ValueError('Review archive mismatch '+n)
 return len(m)

def state_path():return Path.home()/'rk325-v1-latest.json'
def _state(state):
 save(state_path(),state)
 if state.get('work'):save(Path(state['work'])/'submission.json',state)
 if state.get('job'):save(Path.home()/f'rk325-v1-job-{state["job"]}.json',state)

def submit(kit):
 kit=Path(kit).resolve();lock=Path.home()/'.rk325-v1-launch.lock';lock.mkdir()
 try:
  if state_path().exists():raise ValueError('Stage325 already has saved state; collect/review it before another launch')
  intent=Path.home()/'rk325-v1-launch-attempt.json'
  if intent.exists():raise ValueError('Prior submission attempt exists; inspect scheduler and saved intent before retrying')
  work=Path(tempfile.mkdtemp(prefix='rk325-v1-owner-intervention-',dir=Path.home()));before=work/'preflight';before.mkdir()
  state=dict(stage=325,job=None,status='preflight',run_id=work.name,kit=str(kit),work=str(work),out=str(work/'results'),log=None,preflight=str(before),intent=str(intent),external_fixture=None)
  _state(state)
  try:
   from build_environment313 import controlled_environment
   from toolchain313 import inspect
   preflight(kit)
   with controlled_environment(kit,before):inspect(kit,before)
   from continuation_binding import verify_binding
   binding=verify_binding(kit);save(before/'prior-binding.json',binding)
  except BaseException as exc:
   state.update(status='preflight_failed',error=repr(exc));_state(state)
   raise RuntimeError('No job submitted; preflight failed. Run collect.sh to recover evidence: '+str(before)) from exc
  cmd=['sbatch','--parsable','--chdir='+str(Path.home()),'--output='+str(Path.home()/'rk325-v1-%j.log'),str(kit/'qualification.sbatch'),str(kit),str(work/'results')]
  attempt=dict(command=cmd,status='submitting',out=state['out'],run_id=state['run_id']);save(intent,attempt)
  state['status']='submitting';_state(state)
  try:p=subprocess.run(cmd,text=True,capture_output=True)
  except BaseException as exc:
   attempt.update(status='outcome_unknown',error=repr(exc));save(intent,attempt);state.update(status='submission_unknown',error=repr(exc));_state(state);raise
  attempt.update(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr);save(intent,attempt)
  answer=p.stdout.strip().split(';')[0]
  if p.returncode or not answer.isdigit() or int(answer)<=0:
   state.update(status='submission_unknown',error='sbatch result did not establish a successful job identity');_state(state)
   raise RuntimeError('Submission outcome unconfirmed; inspect '+str(intent)+' and scheduler before retrying. '+p.stderr)
  state.update(status='submitted',job=int(answer),log=str(Path.home()/f'rk325-v1-{answer}.log'));_state(state)
  attempt.update(status='submitted',job=state['job']);save(intent,attempt)
  print('Submitted '+answer+': exact native gates, bound private seed reuse, bounded Stage325 qualification; collect even if the job fails',flush=True)
 finally:lock.rmdir()

def status(state):
 print('controller:',state.get('status','submitted'))
 if state.get('job') is None:
  print('No confirmed job ID; evidence:',state.get('work'));print(state.get('error',''));return
 print('active:',queued(state['job']))
 p=subprocess.run(['sacct','-j',str(state['job']),'-X','-o','JobID,State,ExitCode,Start,End,Elapsed,NodeList'],capture_output=True,text=True);print(p.stdout)
 if state.get('log'):
  log=Path(state['log'])
  if log.exists():print('\n'.join(log.read_text(errors='replace').splitlines()[-30:]))

MAX_MEMBER=16*1024**2
MAX_TOTAL=64*1024**2

def export_files(state):
 """Inventory selected small metadata only; never read capture/checkpoint bytes."""
 entries=[];excluded=[];total=0
 def add(path,name):
  nonlocal total
  path=Path(path)
  if path.is_symlink():excluded.append(dict(path=str(path),reason='symlink_not_followed'));return
  st=path.stat()
  if st.st_size>MAX_MEMBER:raise ValueError('Metadata member exceeds16MiB; original retained: '+str(path))
  total+=st.st_size
  if total>MAX_TOTAL:raise ValueError('Metadata collection exceeds64MiB; originals retained; do not full-copy results')
  entries.append((path,name,st.st_size,st.st_mtime_ns))
 roots=[(Path(state['kit']),'kit'),(Path(state['preflight']),'preflight'),(Path(state['out']),'results')]
 for root,label in roots:
  if not root.exists():continue
  for path in sorted(root.rglob('*')):
   rel=path.relative_to(root)
   if '__pycache__' in rel.parts:continue
   if label=='results' and rel.parts[0]=='sources':continue
   if path.is_symlink():excluded.append(dict(path=str(path),reason='symlink_not_followed'));continue
   if not path.is_file():continue
   if label=='results' and (path.name.startswith(('execution316-','test_checkpoint_')) or (path.name.startswith('work325-rank-') and path.suffix=='.tsv') or path.name=='work325-rank-vectors.json' or 'grids' in rel.parts or path.suffix not in ('.json','.txt','.log','.tsv')):
    st=path.stat();excluded.append(dict(path=str(path),reason='native_artifact_retained_no_read',stored_bytes=st.st_size,mtime_ns=st.st_mtime_ns));continue
   add(path,label+'/'+str(rel))
 for key,name in [('log','slurm.log'),('intent','launch-attempt.json')]:
  if state.get(key) and Path(state[key]).is_file():add(Path(state[key]),name)
 return entries,excluded,total

def collect(state):
 import io
 job=state.get('job')
 if job is not None:
  if queued(job):raise ValueError('Job still active')
 elif state.get('status')!='preflight_failed':raise RuntimeError('No confirmed terminal job; inspect saved intent/scheduler')
 token=str(job) if job is not None else 'preflight-'+state['run_id']
 if '/' in token or '..' in token:raise ValueError('Unsafe archive identity')
 archive=Path.home()/f'rk325-review-{token}.tar.gz'
 if not archive.exists():
  entries,excluded,total=export_files(state)
  if job is not None:
   p=subprocess.run(['sacct','-j',str(job),'-X','--parsable2','-o','JobID,JobName,State,ExitCode,Start,End,Elapsed,ElapsedRaw,NodeList'],capture_output=True,text=True)
   accounting=dict(command=p.args,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
  else:accounting=dict(status='not_submitted',reason='preflight failed')
  try:verify(Path(state['kit']));seal=dict(passed=True)
  except Exception as exc:seal=dict(passed=False,error=repr(exc))
  generated={'submission.json':state,'accounting.json':accounting,'kit-seal-check.json':seal,
   'collection-notes.json':dict(policy='Compact metadata only; checkpoint/capture/binary/raw work TSV/full peer-vector bytes remain native. Bound summaries and separate compact review counters are exported. No staging tree.',selected_bytes=total,limit_bytes=MAX_TOTAL,excluded=excluded)}
  temp=archive.with_suffix('.tar.tmp')
  if temp.exists():raise ValueError('Prior collection temporary exists; inspect before retry: '+str(temp))
  manifest={}
  def encoded(obj):return (json.dumps(obj,indent=2,allow_nan=False)+'\n').encode()
  try:
   with tarfile.open(temp,'x:gz') as t:
    def add_bytes(name,data):
     info=tarfile.TarInfo('review/'+name);info.size=len(data);info.mode=0o644;t.addfile(info,io.BytesIO(data));manifest[name]=hashlib.sha256(data).hexdigest()
    for name,obj in generated.items():add_bytes(name,encoded(obj))
    for i,(path,name,size,mtime) in enumerate(entries,1):
     if path.is_symlink():raise ValueError('Collection source became symlink')
     data=path.read_bytes();st=path.stat()
     if len(data)!=size or st.st_size!=size or st.st_mtime_ns!=mtime:raise ValueError('Metadata changed during collection: '+str(path))
     add_bytes(name,data)
     if i%500==0:print(f'Collected{i}/{len(entries)} compact files',flush=True)
    payload=encoded(manifest);info=tarfile.TarInfo('review/files-sha256.json');info.size=len(payload);t.addfile(info,io.BytesIO(payload))
   verify_archive(temp);temp.rename(archive)
  except BaseException:
   temp.unlink(missing_ok=True);raise
 print('Verified',verify_archive(archive),'files; compact archive',archive,'SHA256',sha(archive),flush=True)
 try:subprocess.run(['scp',str(archive),'kevlahan@murier.math.mcmaster.ca:~/Desktop/'],check=True)
 except (OSError,subprocess.CalledProcessError) as exc:raise RuntimeError('Archive retained. Rerun collect.sh to retry transfer; do not resubmit.') from exc
 return archive

if __name__=='__main__':
 if len(sys.argv)!=3 or sys.argv[1] not in ('submit','status','collect'):raise SystemExit('Usage: control325.py submit|status|collect KIT')
 action=sys.argv[1];kit=Path(sys.argv[2]).resolve()
 if action=='submit':submit(kit)
 else:(status if action=='status' else collect)(json.loads(state_path().read_text()))
