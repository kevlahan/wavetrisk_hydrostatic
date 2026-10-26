import json,tempfile,unittest
from pathlib import Path
from review_work325 import project

class HardwareTests(unittest.TestCase):
 def fixture(self,p):
  rows=[]
  for rank,peer in [(0,1),(1,0)]:
   rows.append({'rank':rank,'sha256':'0'*64,'records':[{'tag':'route.exchange','phase':700,'calls':1,'sum':[3],'max':[3],'peers':{str(peer):[24,24,24,24]}}]})
  (p/'work325-rank-vectors.json').write_text(json.dumps({'per_rank':rows}))
  return {'mapping':[{'rank':i,'cpu':i,'core':i,'socket':i,'node':i,'llc':[i]} for i in range(2)]}
 def test_cross_socket_exact_and_compact(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);placement=self.fixture(p);result=project(p,placement)
   self.assertTrue(result['passed']);r=json.loads((p/'work325-hardware-review.json').read_text())
   self.assertEqual(r['endpoint_burdens'][0]['global_byte_totals'],[48,48,0,0,48,48,48,48])
   rank=json.loads((p/'work325-review-rank-0.json').read_text());self.assertNotIn('peers',rank)
 def test_wrong_rank_and_nonconservation_fail(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);placement=self.fixture(p);placement['mapping'][0]['rank']=1
   with self.assertRaises(ValueError):project(p,placement)
   placement=self.fixture(p);f=p/'work325-rank-vectors.json';o=json.loads(f.read_text());o['per_rank'][0]['records'][0]['peers']['1'][0]=25;f.write_text(json.dumps(o))
   with self.assertRaises(ValueError):project(p,placement)

if __name__=='__main__':unittest.main()
