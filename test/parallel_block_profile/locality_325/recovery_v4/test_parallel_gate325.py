from pathlib import Path
import copy,json,os,sys,tempfile,time,unittest
import parallel_gate325 as p

class PairTests(unittest.TestCase):
 def payload(self,arm):
  f={'physical-key':'value-digest'};r={'provider-family':'proof-digest',**{'family'+str(i):'proof-digest' for i in range(253)}}
  import inverse_identity325 as inv
  from test_inverse_bridge_adversary325 import schedule
  inv.bind_inverse_bridge_sources325({'candidate':inv.EXPECTED_SOURCE_HASHES325})
  scopes,plan=schedule();context=inv.inverse_bridge_context325(scopes,plan)
  # Aggregator tests use compact synthetic publication tokens. The source-bound
  # context builder itself is covered by the separate adversarial bridge tests.
  context['publication']={k:p.digest(v) for k,v in context['publication'].items()}
  bridges={}
  for rank in range(256):
   for step in range(1,5):
    c=copy.deepcopy(context);c['step']=c['epoch']=step;bridges[f'{rank}:{step}']=c
  return dict(bridges=bridges,arm=arm,rank_count=256,feature_count=1,features=f,proofs=r,counts=dict(families=254),features_sha256=p.digest(f),proofs_sha256=p.digest(r))
 def run_pair(self,a,b,expect=None):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);left=root/'a.json';right=root/'b.json';out=root/'result.json'
   p.save(left,a);p.save(right,b)
   if expect:
    with self.assertRaisesRegex(ValueError,expect):p.compare_payloads(left,right,out)
    return json.loads(out.read_text()) if out.exists() else None
   return p.compare_payloads(left,right,out)
 def test_equal(self):self.assertTrue(self.run_pair(self.payload('qual-C0'),self.payload('qual-L256'))['passed'])
 def test_value_difference_retained(self):
  a=self.payload('qual-C0');b=self.payload('qual-L256');b['features']['physical-key']='wrong';b['features_sha256']=p.digest(b['features'])
  r=self.run_pair(a,b,'physical feature/provider');self.assertEqual(r['unequal_features'],['physical-key'])
 def test_provider_difference_retained(self):
  a=self.payload('qual-C0');b=self.payload('qual-L256');b['proofs']['provider-family']='wrong';b['proofs_sha256']=p.digest(b['proofs'])
  r=self.run_pair(a,b,'physical feature/provider');self.assertEqual(r['unequal_providers'],['provider-family'])
 def test_missing_and_extra_coverage(self):
  a=self.payload('qual-C0');b=self.payload('qual-L256');b['features']={'other-key':'other'};b['features_sha256']=p.digest(b['features'])
  r=self.run_pair(a,b,'physical feature/provider');self.assertEqual(r['missing_features'],['physical-key']);self.assertEqual(r['extra_features'],['other-key'])
 def test_inverse_bridge_payload_is_required(self):
  a=self.payload('qual-C0');b=self.payload('qual-L256');del b['bridges']['0:1']
  self.run_pair(a,b,'bridge census incomplete')
 def test_missing_provider_family_rejected(self):
  a=self.payload('qual-C0');b=self.payload('qual-L256');del b['proofs']['family0'];b['proofs_sha256']=p.digest(b['proofs'])
  self.run_pair(a,b,'provider-family census')
 def test_tamper_and_rank_census(self):
  for key,value in [('rank_count',255),('feature_count',0),('proofs_sha256','tampered'),('features_sha256','tampered')]:
   with self.subTest(key=key):
    a=self.payload('qual-C0');b=self.payload('qual-L256');b[key]=value
    self.run_pair(a,b,'Incomplete|digest mismatch')

class SupervisionTests(unittest.TestCase):
 @unittest.skipUnless(Path('/proc/self/status').exists(),'Linux target resource accounting')
 def test_memory_accounts_root(self):
  r=p.process_memory(os.getpid());self.assertGreater(r['rss_bytes'],0);self.assertGreaterEqual(r['virtual_bytes'],r['rss_bytes']);self.assertGreaterEqual(r['pids'],1)
 def test_success_and_child_failure_cleanup(self):
  # Use exact PID evidence to establish descendant cleanup, not just parent exit.
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);childpid=root/'descendant.json';log=root/'long.log'
   descendant="import signal,sys,time;from pathlib import Path;marker=Path(sys.argv[1]);signal.signal(signal.SIGTERM,lambda a,b:(marker.with_suffix('.terminated').write_text('SIGTERM'),sys.exit(0)));marker.write_text('ready');time.sleep(60)"
   longcode="import subprocess,sys,json,time; from pathlib import Path; c=subprocess.Popen([sys.executable,'-c',sys.argv[2],sys.argv[1]+'.ready']);Path(sys.argv[1]).write_text(json.dumps(c.pid));deadline=time.monotonic()+5\nwhile not Path(sys.argv[1]+'.ready').exists() and time.monotonic()<deadline:time.sleep(.01)\ntime.sleep(60)"
   failcode="from pathlib import Path;import sys,time;p=Path(sys.argv[1]);deadline=time.monotonic()+5\nwhile not Path(str(p)+'.ready').exists() and time.monotonic()<deadline:time.sleep(.01)\nsys.exit(7)"
   # resource_sample is Linux only; the portable group-cleanup test replaces only
   # that OS accounting boundary, never the supervisor or spawned process calls.
   from unittest.mock import patch
   with patch.object(p,'resource_sample',return_value={}):
    with self.assertRaisesRegex(ValueError,'child failed'):
     p.supervise([([sys.executable,'-c',longcode,str(childpid),descendant],log),([sys.executable,'-c',failcode,str(childpid)],root/'fail.log')],root,{},None)
   self.assertTrue(childpid.exists())
   marker=Path(str(childpid)+'.terminated');deadline=time.monotonic()+2
   while not marker.exists() and time.monotonic()<deadline:time.sleep(.02)
   self.assertEqual(marker.read_text(),'SIGTERM','A worker descendant did not receive process-group termination')
 def test_scratch_limit(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);f=root/'large';f.touch()
   with f.open('r+b') as out:out.truncate(p.CAP_SCRATCH+1)
   from unittest.mock import patch
   with patch.object(p,'process_memory',return_value={'rss_bytes':1,'virtual_bytes':1}):
    with self.assertRaisesRegex(ValueError,'scratch cap'):p.resource_sample(root,{})
 def test_aggregate_memory_limit(self):
  with tempfile.TemporaryDirectory() as tmp:
   from unittest.mock import patch
   with patch.object(p,'process_memory',return_value={'rss_bytes':1,'virtual_bytes':p.CAP_MEMORY+1}):
    with self.assertRaisesRegex(ValueError,'virtual-memory cap'):p.resource_sample(Path(tmp),{})

if __name__=='__main__':unittest.main()
