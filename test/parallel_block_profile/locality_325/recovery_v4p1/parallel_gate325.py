"""Bounded process supervision and deterministic cross-arm reconciliation.

The parent enforces the aggregate budget; child groups can be terminated without
signalling original solver jobs. This module never launches a solver or compiler.
"""
from pathlib import Path
import hashlib,json,os,signal,subprocess,sys,time
CAP_SCRATCH=4*1024**3
CAP_MEMORY=16*1024**3

def require(ok,message):
 if not ok:raise ValueError(message)
def save(path,obj):
 path=Path(path);p=path.with_suffix(path.suffix+'.tmp');p.write_text(json.dumps(obj,indent=2)+'\n');p.replace(path)
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def process_memory(root_pid):
 """Conservative sum: shared pages/address space can be counted more than once."""
 rows={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   fields={line.split(':',1)[0]:line.split(':',1)[1].strip() for line in (p/'status').read_text().splitlines() if ':' in line}
   rows[int(p.name)]=(int(fields['PPid']),int(fields.get('VmRSS','0 kB').split()[0])*1024,int(fields.get('VmSize','0 kB').split()[0])*1024)
  except (OSError,ValueError,KeyError):continue
 owned={root_pid};prior=0
 while len(owned)!=prior:
  prior=len(owned);owned.update(pid for pid,v in rows.items() if v[0] in owned)
 return {'pids':len(owned),'rss_bytes':sum(rows.get(p,(0,0,0))[1] for p in owned),'virtual_bytes':sum(rows.get(p,(0,0,0))[2] for p in owned)}
def resource_sample(scratch,stats):
 size=sum(p.stat().st_size for p in Path(scratch).rglob('*') if p.is_file())
 mem=process_memory(os.getpid())
 stats['sampled_peak_scratch_bytes']=max(stats.get('sampled_peak_scratch_bytes',0),size)
 stats['sampled_peak_aggregate_rss_bytes']=max(stats.get('sampled_peak_aggregate_rss_bytes',0),mem['rss_bytes'])
 stats['sampled_peak_aggregate_virtual_bytes']=max(stats.get('sampled_peak_aggregate_virtual_bytes',0),mem['virtual_bytes'])
 require(size<=CAP_SCRATCH,'Aggregate scratch cap exceeded')
 require(mem['virtual_bytes']<=CAP_MEMORY,'Aggregate conservative virtual-memory cap exceeded')
 return mem

def supervise(commands,scratch,stats,notify=None):
 """Each command is (argv, logpath); no shell, no detached descendants."""
 processes=[];streams=[]
 try:
  for argv,logpath in commands:
   f=Path(logpath).open('wb');streams.append(f)
   processes.append(subprocess.Popen(argv,stdout=f,stderr=subprocess.STDOUT,start_new_session=True))
  previous=time.monotonic()
  while True:
   states=[p.poll() for p in processes]
   require(all(s is None or s==0 for s in states),'Verification child failed; inspect retained phase log')
   if all(s==0 for s in states):break
   resource_sample(scratch,stats)
   if notify and time.monotonic()-previous>=5:notify();previous=time.monotonic()
   time.sleep(.25)
  resource_sample(scratch,stats)
 finally:
  # Kill the exact groups we created, including workers of a failed coordinator.
  for p in processes:
   try:os.killpg(p.pid,signal.SIGTERM)
   except ProcessLookupError:pass
  deadline=time.monotonic()+2
  while any(p.poll() is None for p in processes) and time.monotonic()<deadline:time.sleep(.05)
  for p in processes:
   try:os.killpg(p.pid,signal.SIGKILL)
   except ProcessLookupError:pass
   p.wait()
  for f in streams:f.close()

def compare_payloads(left,right,report_path):
 """Retain a full mismatch census and both proof digests even on failure."""
 arms=[json.loads(Path(p).read_text()) for p in (left,right)]
 require([a['arm'] for a in arms]==['qual-C0','qual-L256'],'Unexpected optimized pair identities')
 report={'passed':False,'arms':{a['arm']:{k:a[k] for k in ('counts','feature_count','features_sha256','proofs_sha256','rank_count')} for a in arms}}
 fa,fb=[a['features'] for a in arms];pa,pb=[a['proofs'] for a in arms]
 bridge_errors=[];bridges=[]
 # Every native arm must supply the complete1024row bridge census.
 if True:
  import inverse_identity325 as inv
  ba,bb=[a['bridges'] for a in arms]
  expected={f'{rank}:{step}' for rank in range(256) for step in range(1,5)}
  require(ba.keys()==bb.keys()==expected,'Physical inverse bridge census incomplete')
  for key in sorted(expected):
   try:bridges.append(inv.compare_inverse_bridge325(ba[key],bb[key]))
   except Exception as exc:bridge_errors.append(dict(key=key,error=repr(exc)))
  report['inverse_bridge']=dict(passed=not bridge_errors,contexts=len(expected),proved=len(bridges),moved_physical_nodes=sum(x['moved_physical_nodes'] for x in bridges),errors=bridge_errors)

 missing=sorted(fa.keys()-fb.keys());extra=sorted(fb.keys()-fa.keys());bad=sorted(k for k in fa.keys()&fb.keys() if fa[k]!=fb[k])
 pmissing=sorted(pa.keys()-pb.keys());pextra=sorted(pb.keys()-pa.keys());pbad=sorted(k for k in pa.keys()&pb.keys() if pa[k]!=pb[k])
 report.update(missing_features=missing,extra_features=extra,unequal_features=bad,missing_providers=pmissing,extra_providers=pextra,unequal_providers=pbad)
 for a in arms:
  require(a['counts']['families']==254 and len(a['proofs'])==254,'Incomplete provider-family census')
  require(a['rank_count']==256 and a['feature_count']==len(a['features']),'Incomplete provider-arm payload')
  require(digest(a['features'])==a['features_sha256'] and digest(a['proofs'])==a['proofs_sha256'],'Arm payload digest mismatch')
 report['passed']=not any((missing,extra,bad,pmissing,pextra,pbad,bridge_errors));save(report_path,report)
 require(report['passed'],'Cross-arm physical feature/provider mismatch; complete census retained')
 return report

def provider_main(database,arm,output,progress):
 import verify325 as q
 import inverse_identity325 as inv
 kit=Path(__file__).resolve().parent;binding=json.loads((kit/'binding325.json').read_text())
 manifest=Path(binding['original_kit'])/'source-files.json'
 require(hashlib.sha256(manifest.read_bytes()).hexdigest()==binding['source_manifest_sha256'],'Original source manifest changed')
 inv.bind_inverse_bridge_sources325(json.loads(manifest.read_text()),source_root=Path(binding['original_results'])/'sources/candidate')
 c=q.sqlite3.connect(database)
 def emit(stage,arm,completed,total,detail=None,persist=False):
  if persist:save(progress,dict(stage=stage,arm=arm,completed=completed,total=total,detail=detail))
 try:
  ranks=[x[0] for x in c.execute('select rank from domain where arm=? order by rank',(arm,))]
  require(ranks==list(range(256)),'Incomplete provider rank census')
  proofs,counts=q.route_certificates(c,arm,emit)
  require(counts['families']==len(proofs)==254,'Incomplete provider family census')
  features=dict(c.execute('select key,digest from feature where arm=?',(arm,)))
  bridges={f'{rank}:{key}':json.loads(payload) for rank,key,payload in c.execute("select rank,key,payload from context where arm=? and kind='inverse-bridge'",(arm,))}
  require(len(bridges)==1024,'Incomplete inverse physical bridge contexts')
  result=dict(bridges=bridges,arm=arm,rank_count=len(ranks),features=features,feature_count=len(features),proofs=proofs,counts=counts,features_sha256=digest(features),proofs_sha256=digest(proofs))
  save(output,result)
 finally:c.close()
if __name__=='__main__':
 require(len(sys.argv)==6 and sys.argv[1]=='providers','Usage: parallel_gate325.py providers DB ARM OUTPUT PROGRESS')
 provider_main(*sys.argv[2:])
