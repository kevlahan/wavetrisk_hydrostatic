"""Single verification-only recovery submission and compact local archive."""
from pathlib import Path
import getpass,hashlib,io,json,subprocess,sys,tarfile,tempfile
from recover325 import package_check,require,save,sha

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

def state_path():return Path.home()/'rk325-recovery-v4p1-latest.json'
def record(state):
 save(state_path(),state)
 if state.get('work'):save(Path(state['work'])/'submission.json',state)
 if state.get('job'):save(Path.home()/f'rk325-recovery-v4p1-job-{state["job"]}.json',state)

def submit(kit):
 kit=Path(kit).resolve(strict=True);package_check(kit)
 require(not queued(55569),'Original job55569 still active')
 require(not queued(55571),'Earlier recovery55571 still active')
 require(not queued(55572),'Earlier recovery55572 still active')
 require(not queued(55573),'Predecessor recovery55573 still active')
 lock=Path.home()/'.rk325-recovery-v4p1-launch.lock';lock.mkdir()
 try:
  require(not state_path().exists(),'Recovery already has state; inspect it, no automatic retry')
  intent=Path.home()/'rk325-recovery-v4p1-launch-attempt.json';require(not intent.exists(),'Prior recovery submission outcome requires inspection')
  work=Path(tempfile.mkdtemp(prefix='rk325-verification-recovery-v4p1-',dir=Path.home()))
  state=dict(stage=325,revision='verification-recovery-v4p1',source_job=55569,predecessor_job=55573,status='submitting',job=None,kit=str(kit),work=str(work),out=str(work/'results'),intent=str(intent),log=None)
  record(state)
  cmd=['sbatch','--parsable','--chdir='+str(Path.home()),'--output='+str(Path.home()/'rk325-recovery-v4p1-%j.log'),str(kit/'verification.sbatch'),str(kit),state['out']]
  save(intent,dict(status='submitting',command=cmd))
  try:p=subprocess.run(cmd,text=True,capture_output=True)
  except BaseException as exc:
   state.update(status='submission_unknown',error=repr(exc));record(state);raise
  answer=p.stdout.strip().split(';')[0];save(intent,dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
  if p.returncode or not answer.isdigit() or int(answer)<=0:
   state.update(status='submission_unknown',error='No confirmed job ID');record(state);raise RuntimeError('Submission unconfirmed; inspect intent and scheduler; no retry')
  state.update(status='submitted',job=int(answer),log=str(Path.home()/f'rk325-recovery-v4p1-{answer}.log'));record(state)
  print('Submitted recovery',answer,'one task16CPUs,5minute cap,zero builds/solver launches',flush=True)
 finally:lock.rmdir()

def status(state):
 print('controller:',state['status'],'job:',state.get('job'))
 if state.get('job'):
  print('active:',queued(state['job']))
  p=subprocess.run(['sacct','-j',str(state['job']),'-X','-o','JobID,State,ExitCode,Start,End,Elapsed,NodeList'],capture_output=True,text=True);print(p.stdout)
 for name in ('recovery-report.json','cross-arm.json','qual-C0-provider-progress.json','qual-L256-provider-progress.json'):
  path=Path(state['out'])/name
  if path.is_file():
   obj=json.loads(path.read_text());print(name,{k:obj[k] for k in ('status','phase','active_arm','active_rank','active_family','elapsed_seconds','passed','error') if k in obj})
 if state.get('log') and Path(state['log']).is_file():print('\n'.join(Path(state['log']).read_text(errors='replace').splitlines()[-25:]))

def selected_files(state):
 files=[];total=0
 for root,label in ((Path(state['kit']),'kit'),(Path(state['out']),'results')):
  if not root.exists():continue
  for p in sorted(root.rglob('*')):
   rel=p.relative_to(root)
   if '__pycache__' in rel.parts or 'scratch' in rel.parts or p.is_symlink() or not p.is_file():continue
   if p.name=='work325-rank-vectors.json' or p.suffix=='.tsv':continue
   size=p.stat().st_size;require(size<=16*1024**2,'Recovery return member exceeds16MiB')
   total+=size;require(total<=64*1024**2,'Recovery compact return exceeds64MiB; original evidence retained')
   files.append((p,label+'/'+str(rel),size,p.stat().st_mtime_ns))
 if state.get('log') and Path(state['log']).is_file():
  p=Path(state['log']);s=p.stat();require(s.st_size<=16*1024**2,'Recovery log exceeds16MiB');total+=s.st_size;files.append((p,'slurm.log',s.st_size,s.st_mtime_ns))
 require(total<=64*1024**2,'Recovery compact return exceeds64MiB')
 return files,total

def verify_archive(path):
 with tarfile.open(path) as t:
  members=t.getmembers();names=[m.name for m in members]
  require(len(names)==len(set(names)) and all(m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in members),'Unsafe review archive')
  m=json.load(t.extractfile('review/files-sha256.json'))
  require(set(names)=={'review/files-sha256.json',*('review/'+n for n in m)},'Review inventory differs')
  for n,h in m.items():require(hashlib.file_digest(t.extractfile('review/'+n),'sha256').hexdigest()==h,'Review archive changed: '+n)
 return len(m)

def collect(state):
 require(state.get('job') and not queued(state['job']),'Confirmed terminal recovery job required')
 archive=Path.home()/f'rk325-replay-v4p1-review-{state["job"]}.tar.gz'
 if not archive.exists():
  files,total=selected_files(state)
  p=subprocess.run(['sacct','-j',str(state['job']),'-X','--parsable2','-o','JobID,JobName,State,ExitCode,Start,End,Elapsed,ElapsedRaw,NodeList'],capture_output=True,text=True)
  generated={'submission.json':state,'accounting.json':dict(command=p.args,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr),'collection.json':dict(selected_bytes=total,maximum_bytes=64*1024**2,policy='Recovery reports/helpers only. Original captured bytes, checkpoints and full work vectors are omitted and remain native.')}
  temp=Path(str(archive)+'.tmp');require(not temp.exists(),'Prior recovery archive temporary exists; inspect it')
  manifest={}
  try:
   with tarfile.open(temp,'x:gz') as t:
    def add(name,data):
     info=tarfile.TarInfo('review/'+name);info.size=len(data);t.addfile(info,io.BytesIO(data));manifest[name]=hashlib.sha256(data).hexdigest()
    for n,o in generated.items():add(n,(json.dumps(o,indent=2)+'\n').encode())
    for path,name,size,mtime in files:
     data=path.read_bytes();st=path.stat();require(not path.is_symlink() and st.st_size==size and len(data)==size and st.st_mtime_ns==mtime,'Recovery report changed during collection');add(name,data)
    data=(json.dumps(manifest,indent=2)+'\n').encode();info=tarfile.TarInfo('review/files-sha256.json');info.size=len(data);t.addfile(info,io.BytesIO(data))
   verify_archive(temp);temp.rename(archive)
  except BaseException:temp.unlink(missing_ok=True);raise
 print('Verified',verify_archive(archive),'files;',archive,'SHA256',sha(archive),flush=True)
 return archive

if __name__=='__main__':
 require(len(sys.argv)==3 and sys.argv[1] in ('submit','status','collect'),'Usage: control_recovery325.py submit|status|collect KIT')
 if sys.argv[1]=='submit':submit(Path(sys.argv[2]))
 else:(status if sys.argv[1]=='status' else collect)(json.loads(state_path().read_text()))
