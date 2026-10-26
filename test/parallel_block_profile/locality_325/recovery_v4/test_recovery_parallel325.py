from pathlib import Path
import copy,json,shutil,tempfile,unittest
from unittest.mock import patch
import recover_parallel325 as r

class HistoryTests(unittest.TestCase):
 def fixture(self,root):
  kit=Path(__file__).parent
  for name in ('parallel-history325.json','recovery-history325.json','prior-native-input-anchors.json','predecessor-v3'):
   src=kit/name;dst=root/name
   if src.is_dir():shutil.copytree(src,dst)
   else:shutil.copy2(src,dst)
  h=json.loads((root/'parallel-history325.json').read_text());h['remote_predecessor_files']={str((root/'predecessor-v3'/Path(name).name).resolve()):digest for name,digest in h['remote_predecessor_files'].items()}
  r.save(root/'parallel-history325.json',h)
  return json.loads((kit/'binding325.json').read_text())
 def test_complete_historical_binding(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);binding=self.fixture(root);h,anchors,s=r.history_check(root,binding)
   self.assertEqual(h['cumulative_maximum_seconds'],2692);self.assertEqual(len(anchors),1536);self.assertEqual(s['runs']['checked-L256']['status'],'complete')
 def test_bad_cumulative_budget_rejected(self):
  for key,value in [('walltime_cap_seconds',600),('prior_elapsed_seconds',881),('predecessor_job',55572),('predecessor_elapsed_seconds',1510)]:
   with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);binding=self.fixture(root);h=json.loads((root/'parallel-history325.json').read_text());h[key]=value;r.save(root/'parallel-history325.json',h)
    with self.assertRaises(ValueError):r.history_check(root,binding)
 def test_checked_reuse_requires_complete_actual_certificate(self):
  for mutation in ('streams','providers','reads','scope'):
   with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);binding=self.fixture(root)
    name='predecessor-v3/'+('semantic-replay.json' if mutation in ('streams','providers') else 'capture-read-checks-partial.json' if mutation=='reads' else 'recovery-report.json')
    p=root/name;j=json.loads(p.read_text())
    if mutation=='streams':j['runs']['checked-L256']['streams'].pop()
    elif mutation=='providers':j['runs']['checked-L256']['families']=253
    elif mutation=='reads':j.pop(next(iter(j)))
    else:j['solver_runs']=1
    r.save(p,j);h=json.loads((root/'parallel-history325.json').read_text());h['certificate_sha256'][name]=r.sha(p);h['remote_predecessor_files'][str(p.resolve())]=r.sha(p);r.save(root/'parallel-history325.json',h)
    with self.assertRaises(ValueError):r.history_check(root,binding)
 def test_changed_native_predecessor_certificate_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);binding=self.fixture(root);h=json.loads((root/'parallel-history325.json').read_text());key=next(iter(h['remote_predecessor_files']));h['remote_predecessor_files'][key]='0'*64;r.save(root/'parallel-history325.json',h)
   with self.assertRaisesRegex(ValueError,'Original55573certificate'):r.history_check(root,binding)
 def test_exact_allocation(self):
  env=dict(SLURM_JOB_ID='99999',SLURM_JOB_PARTITION='nk',SLURM_JOB_NUM_NODES='1',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='16')
  with patch.dict(r.os.environ,env,clear=True),patch.object(r.socket,'gethostname',return_value='bb13.rhpcs.internal'),patch.object(r.os,'sched_getaffinity',create=True,return_value=set(range(16))),patch.object(r.resource,'getrlimit',return_value=(16*1024**3,16*1024**3)):
   r.allocation()
   with patch.object(r.os,'sched_getaffinity',create=True,return_value=set(range(8))):
    with self.assertRaisesRegex(ValueError,'16allocated'):r.allocation()

if __name__=='__main__':unittest.main()
