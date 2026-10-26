"""Read-only replay of job55569 captures; never build or execute a solver."""
from pathlib import Path
from contextlib import contextmanager
import gzip,hashlib,io,json,os,resource,signal,socket,sys,time,platform
LABELS=('checked-L256','qual-C0','qual-L256')
MAX_RETURN=64*1024**2

def require(ok,message):
 if not ok:raise ValueError(message)
def save(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('w') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 tmp.replace(path)
def sha(path,links=False):
 path=Path(path)
 if path.is_symlink():
  require(links,'Unexpected symbolic link: '+str(path));return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
 with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def safe(root,name):
 rel=Path(name);require(not rel.is_absolute() and '..' not in rel.parts and name==str(rel),'Unsafe relative path')
 root=Path(root).resolve(strict=True);p=root/rel
 require(p.parent.resolve(strict=True).is_relative_to(root),'Escaping parent path')
 return p

def package_check(kit):
 kit=Path(kit).resolve(strict=True);manifest=json.loads((kit/'files-sha256.json').read_text())
 actual={str(p.relative_to(kit)) for p in kit.rglob('*') if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts and p.name!='files-sha256.json'}
 require(actual==set(manifest),'Recovery package inventory differs')
 for name,h in manifest.items():require(sha(safe(kit,name))==h,'Recovery package changed: '+name)
 return sha(kit/'files-sha256.json')

def metadata_check(binding):
 root=Path(binding['original_results']).resolve(strict=True);kit=Path(binding['original_kit']).resolve(strict=True)
 require(sha(kit/'files-sha256.json')==binding['original_kit_manifest_sha256'],'Original kit manifest changed')
 km=json.loads((kit/'files-sha256.json').read_text())
 for name,h in km.items():require(sha(safe(kit,name),True)==h,'Original kit changed: '+name)
 require(sha(kit/'source-files.json')==binding['source_manifest_sha256'],'Original source manifest changed')
 count=0
 for family,files in json.loads((kit/'source-files.json').read_text()).items():
  for name,h in files.items():require(sha(safe(root/'sources'/family,name),True)==h,'Original source changed: '+name);count+=1
 for name,h in binding['metadata_hashes'].items():
  family,relative=name.split('/',1);base=root if family=='results' else root.parent/'preflight'
  require(sha(safe(base,relative))==h,'Original recorded metadata changed: '+name)
 ledger=json.loads((root/'results.json').read_text())
 require(ledger['status']=='failed' and ledger['error']==binding['expected_failure'] and ledger['qualification'] is None,'Original failure ledger changed')
 require(ledger['performance_measurement'] is False and ledger['performance_qualified'] is False and ledger['schedule']==list(LABELS),'Original qualification scope changed')
 builds=json.loads((root/'builds.json').read_text());require(builds['status']=='passed' and builds==ledger['builds'],'Original build gate changed')
 require(builds['toolchain']['status']=='passed' and builds['toolchain']['expected']['gcc_version']=='13.2.0' and builds['toolchain']['expected']['mpi_version']=='4.0','Recorded native toolchain gate did not pass')
 binaries=[];checkpoints=[]
 for kind,row in builds['builds'].items():
  path=safe(root,row['binary']);require(sha(path)==row['binary_sha256'],'Original executable changed: '+kind)
  binaries.append(dict(kind=kind,path=str(path),sha256=row['binary_sha256']))
 for label in LABELS:
  row=ledger['runs'][label]
  require(row['status']=='passed' and row['ranks']==256 and row['steps']==4,'Incomplete original solver run: '+label)
  require(row['work_ledger']['passed'] is True and row['solver_placement']['passed'] is True and all(x['passed'] is True for x in row['science'].values()),'Original science/placement/work gate failed: '+label)
  require(row['binary_sha256']==builds['builds'][row['kind']]['binary_sha256'],'Run executable not bound')
  require(row['owner']==('off' if label=='qual-C0' else 'L256'),'Run owner selector changed')
  for index,h in row['checkpoints'].items():
   path=root/label/f'test_checkpoint_{index}.bin.zst';require(sha(path)==h,'Original checkpoint changed')
   checkpoints.append(dict(path=str(path),sha256=h))
  seed=root/label/'test_checkpoint_0003.bin.zst';require(sha(seed)==row['seed_sha256'],'Original private seed changed')
  vectors=root/label/'work325-rank-vectors.json';require(sha(vectors)==row['work_ledger']['per_rank_artifact']['sha256'],'Original work vector certificate changed')
 return dict(passed=True,original_kit_files=len(km),source_files=count,metadata_files=len(binding['metadata_hashes']),binaries=binaries,checkpoints=checkpoints,science_policy='Reuses byte-bound original native decoded science/placement checks; no new solver or independent checkpoint decode.')

def stat_match(path,row):
 require(path.is_file() and not path.is_symlink(),'Expected regular native capture/work file')
 st=path.stat();require(st.st_size==row['stored_bytes'] and st.st_mtime_ns==row['mtime_ns'],'Native capture/work inventory changed: '+str(path))
 return st

def load_history(kit,binding):
 history=json.loads((kit/'recovery-history325.json').read_text())
 require(history['source_job']==55569 and history['predecessor_job']==55571 and history['predecessor_elapsed_seconds']==13,'Unexpected predecessor identity or cost')
 require(history['original_elapsed_seconds']==392 and history['recovery_walltime_cap_seconds']==2100 and history['cumulative_maximum_seconds']==2505,'Recovery cost bound changed')
 prior=history['predecessor_report']
 require(prior['status']=='failed' and prior['error']==history['expected_failure']=="AuditError('physical.transform_owner_range')" and prior['binding_before']['passed'] is True,'Predecessor failure or input binding differs')
 require(prior['verification_job']=='55571' and prior['builds']==prior['solver_runs']==prior['performance_runs']==0,'Predecessor scope differs')
 anchor_path=safe(kit,history['anchors_file']);require(sha(anchor_path)==history['anchors_sha256'],'Historical anchor certificate changed')
 anchors=json.loads(anchor_path.read_text());expected={}
 for family in ('capture_inventory','work_inventory'):
  for row in binding[family]:expected[row['relative']]={**row,'family':family}
 require(len(expected)==1536 and set(anchors)==set(expected),'Historical anchor inventory differs')
 for name,row in anchors.items():
  require(set(row)==set(expected[name])|{'sha256'} and all(row[k]==v for k,v in expected[name].items()),'Historical anchor metadata differs: '+name)
  h=row['sha256'];require(isinstance(h,str) and len(h)==64 and all(c in '0123456789abcdef' for c in h),'Invalid historical SHA256: '+name)
 return history,anchors

def hash_bound_input(path,row,historical):
 stat_match(path,row);h=sha(path);stat_match(path,row)
 name=row['relative'];require(name in historical and historical[name]['sha256']==h,'Native bytes differ from job55571 historical anchor: '+name)
 return h

def anchor_inputs(binding,out,historical):
 root=Path(binding['original_results']).resolve(strict=True);anchors={};started=time.monotonic()
 for inventory,stem in [('capture_inventory','execution316-rank-*.bin.gz'),('work_inventory','work325-rank-*.tsv')]:
  expected={row['relative'] for row in binding[inventory]}
  actual={str(p.relative_to(root)) for label in LABELS for p in (root/label).glob(stem)}
  require(expected==actual and len(expected)==768,'Native input inventory differs: '+inventory)
  for i,row in enumerate(binding[inventory],1):
   path=safe(root,row['relative']);h=hash_bound_input(path,row,historical)
   anchors[row['relative']]={**row,'sha256':h,'family':inventory}
   if i%64==0:print(f'Binding {inventory}: {i}/768; {time.monotonic()-started:.1f}s',flush=True)
 save(out/'native-input-anchors.json',anchors)
 return anchors

class HashingInput:
 def __init__(self,raw):self.raw=raw;self.digest=hashlib.sha256();self.count=0
 def read(self,n=-1):
  data=self.raw.read(n);self.digest.update(data);self.count+=len(data);return data
 def seekable(self):return False
 def tell(self):return self.count

@contextmanager
def verified_stream(path,root,anchors,observed):
 path=Path(path);name=str(path.resolve(strict=True).relative_to(Path(root).resolve(strict=True)))
 require(name in anchors and anchors[name]['family']=='capture_inventory','Unbound semantic input')
 row=anchors[name];stat_match(path,row)
 with path.open('rb') as raw:
  wrapped=HashingInput(raw)
  with gzip.GzipFile(fileobj=wrapped,mode='rb') as stream:
   yield stream
  require(wrapped.count==row['stored_bytes'] and wrapped.digest.hexdigest()==row['sha256'],'Capture changed between binding and verification: '+name)
 stat_match(path,row)
 observed[name]=dict(compressed_sha256=row['sha256'],stored_bytes=row['stored_bytes'],complete=True)

def shadow_work(binding,anchors,out):
 root=Path(binding['original_results']).resolve(strict=True);shadow=out/'work-replay';shadow.mkdir()
 for label in LABELS:(shadow/label).mkdir()
 for row in binding['work_inventory']:
  original=safe(root,row['relative']);stat_match(original,row)
  require(sha(original)==anchors[row['relative']]['sha256'],'Work bytes changed before replay')
  (shadow/row['relative']).symlink_to(original)
 return shadow

def require_allocation():
 require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_PARTITION')=='nk' and os.environ.get('SLURM_JOB_NUM_NODES')=='1' and os.environ.get('SLURM_NTASKS')=='1' and socket.gethostname().split('.')[0]=='bb13','Require one-task bb13/nk allocation')
 require(resource.getrlimit(resource.RLIMIT_AS)==(16*1024**3,16*1024**3),'Require 16-GiB hard and soft virtual-memory limits')

def run(kit,out):
 sys.dont_write_bytecode=True
 kit=Path(kit).resolve(strict=True);package_sha=package_check(kit);binding=json.loads((kit/'binding325.json').read_text());history,historical=load_history(kit,binding);require_allocation()
 out=Path(out).resolve();require(not out.exists(),'Fresh recovery output required')
 original=Path(binding['original_results']).resolve(strict=True);oldkit=Path(binding['original_kit']).resolve(strict=True)
 for protected in (original,oldkit,kit):require(not out.is_relative_to(protected) and not protected.is_relative_to(out),'Recovery output overlaps protected input')
 out.mkdir();scratch=out/'scratch';scratch.mkdir();started=time.monotonic();observed={}
 import numpy as np
 report=dict(stage=325,revision='verification-recovery-v2',source_job=55569,verification_job=os.environ['SLURM_JOB_ID'],status='running',builds=0,solver_runs=0,performance_runs=0,performance_qualified=False,package_manifest_sha256=package_sha,original_elapsed_seconds=392,predecessor_elapsed_seconds=13,recovery_walltime_cap_seconds=2100,cumulative_maximum_seconds=2505)
 report['predecessor_binding']=dict(job=55571,review_archive_sha256=history['predecessor_review_archive_sha256'],anchors_sha256=history['anchors_sha256'],anchor_policy=history['historical_anchor_policy'])
 report['verifier_environment']=dict(python=platform.python_version(),executable=str(Path(sys.executable).resolve()),numpy=np.__version__,hostname=socket.gethostname(),allowed_cpus=sorted(os.sched_getaffinity(0)),requested_tasks=1,scheduler_memory='--mem=0 (site MEMORY=1 placeholder)',process_virtual_memory_limit_bytes=list(resource.getrlimit(resource.RLIMIT_AS)))
 def persist():
  report['elapsed_seconds']=time.monotonic()-started;report['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;save(out/'recovery-report.json',report)
 def stop(signum,frame):raise TimeoutError('Bounded recovery allocation ending; preserve partial evidence')
 signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGALRM,stop);signal.alarm(34*60);persist()
 try:
  report['phase']='bind_original_metadata_source_binaries';persist();report['binding_before']=metadata_check(binding);persist()
  report['phase']='bind_native_capture_and_work_bytes';persist();anchors=anchor_inputs(binding,out,historical);report['historical_input_hashes_matched']=len(anchors);persist()
  report['phase']='corrected_semantic_replay';persist()
  import verify325
  verify325.v.open_stream=lambda path:verified_stream(path,original,anchors,observed)
  semantic=verify325.qualify(original,kit,None,report_path=out/'semantic-replay.json',scratch_root=scratch)
  require(semantic['passed'] is True and len(observed)==768,'Incomplete corrected semantic replay')
  save(out/'capture-read-checks.json',observed);report['semantic_replay']=dict(passed=True,streams=768,report_sha256=sha(out/'semantic-replay.json'),capture_reads_sha256=sha(out/'capture-read-checks.json'));persist()
  report['phase']='work_replay_new_output_only';persist();shadow=shadow_work(binding,anchors,out)
  from account325 import qualify
  work=qualify(shadow,kit);require(work.get('passed') is True,'Work replay failed');save(out/'work-replay.json',work)
  work_reads=0
  for label in LABELS:
   vectors=json.loads((shadow/label/'work325-rank-vectors.json').read_text())['per_rank']
   require([x['rank'] for x in vectors]==list(range(256)),'Incomplete replayed work rank census')
   for rank,row in enumerate(vectors):
    anchor=anchors[f'{label}/work325-rank-{rank}.tsv']
    require(row['sha256']==anchor['sha256'] and row['bytes']==anchor['stored_bytes'],'Work changed during actual replay')
    work_reads+=1
  require(work_reads==768,'Incomplete actual work read bindings')
  report['work_actual_read_hashes_matched']=work_reads
  report['work_replay']=dict(passed=True,report_sha256=sha(out/'work-replay.json'));persist()
  report['phase']='preservation_after';persist();report['binding_after']=metadata_check(binding)
  for row in binding['capture_inventory']+binding['work_inventory']:stat_match(safe(original,row['relative']),row)
  for row in binding['work_inventory']:require(sha(safe(original,row['relative']))==anchors[row['relative']]['sha256'],'Original work changed after replay')
  report.update(status='passed',decision='CORRECTED_QUALIFICATION_PASSED_REQUIRES_INDEPENDENT_REVIEW_NO_PERFORMANCE',phase='complete',original_ledger_preserved=True,capture_compressed_hashes_matched_during_reads=True)
  persist();print('PASS corrected verification-only recovery; original failed ledger preserved, zero solver/build work',flush=True)
 except BaseException as exc:
  save(out/'capture-read-checks-partial.json',observed);report.update(status='failed',error=repr(exc),decision='STOP_RETAIN_PARTIAL_EVIDENCE_NO_AUTOMATIC_RETRY');persist();raise
 finally:
  signal.alarm(0)
  # Only this run's fresh temporary SQLite tree is eligible for cleanup.
  import shutil
  shutil.rmtree(scratch,ignore_errors=True)
if __name__=='__main__':
 require(len(sys.argv)==3,'Usage: recover325.py KIT NEW_OUTPUT');run(*sys.argv[1:])
