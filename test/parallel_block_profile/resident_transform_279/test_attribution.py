import importlib.util,json,tempfile,unittest,contextlib,io
from pathlib import Path
T=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('attribution',T/'attribution.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class AttributionTests(unittest.TestCase):
 def fixture(self,p):
  v={'status':'passed','oracles':False,'ranks':2,'binary_sha256':'fixture','ordinary':{'included':[{'index':i} for i in [1,2,3,4,5,6,9,10,11]],'median':3}}
  (p/'validation.json').write_text(json.dumps(v));lines=['  detail dropped step records = 0']
  for step in range(1,13):
   for rank in range(2):
    # Slowest whole-step rank alternates; summing independent region maxima is wrong.
    wall=3 if rank==step%2 else 2;values=[2,1] if wall==3 else [0,2]
    lines += [f'  detail step {step} rank {rank} restarts 0 remaps 0 wall {wall} CPU 0',f'  detail step {step} rank {rank} self-wall = '+' '.join(map(str,values+[0]*65))]
  (p/'run.log').write_text('\n'.join(lines)+'\n')
 def run_case(self,mutate=None):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);self.fixture(p)
   if mutate:mutate(p)
   with contextlib.redirect_stdout(io.StringIO()):mod.analyze(p)
   return json.loads((p/'attribution.json').read_text())
 def test_same_rank_partition(self):
  r=self.run_case();self.assertEqual(r['mean_slowest_step_seconds'],3);self.assertEqual(r['regions']['1']['mean_on_slowest_step_rank_seconds'],2);self.assertEqual(r['regions']['2']['mean_on_slowest_step_rank_seconds'],1)
 def test_overflow(self):
  with self.assertRaises(ValueError):self.run_case(lambda p:(p/'run.log').write_text((p/'run.log').read_text().replace('records = 0','records = 1')))
 def test_missing_rank(self):
  with self.assertRaises(ValueError):self.run_case(lambda p:(p/'run.log').write_text('\n'.join(l for l in (p/'run.log').read_text().splitlines() if 'rank 1' not in l)))
 def test_non_partition(self):
  with self.assertRaises(ValueError):self.run_case(lambda p:(p/'run.log').write_text((p/'run.log').read_text().replace('self-wall = 2 1','self-wall = 9 1')))
 def test_nan(self):
  with self.assertRaises(ValueError):self.run_case(lambda p:(p/'run.log').write_text((p/'run.log').read_text().replace('self-wall = 2 1','self-wall = NaN 1')))
 def test_oracle_rejected(self):
  def change(p):
   v=json.loads((p/'validation.json').read_text());v['oracles']=True;(p/'validation.json').write_text(json.dumps(v))
  with self.assertRaises(ValueError):self.run_case(change)
if __name__=='__main__':unittest.main()
