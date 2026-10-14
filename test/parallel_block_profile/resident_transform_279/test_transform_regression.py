import tarfile,tempfile,unittest
from pathlib import Path
from transform_contract import evidence,check_global_census
class RealCampaignRegression(unittest.TestCase):
    def test_completed_16_and_83_rank_censuses_match(self):
        archive=Path(__file__).parent/'qualification/job55506-evidence.tar.gz'
        with tempfile.TemporaryDirectory() as tmp:
            with tarfile.open(archive) as t:
                members=[m for m in t.getmembers() if m.isfile() and '/transform279-' in m.name]
                t.extractall(tmp,members=members,filter='data')
            root=Path(tmp)/'results';reports=[]
            for n,index,oracle in [(16,3,False),(16,5,True),(83,6,False)]:
                reports.append(evidence(root/f'J5-n{n}-validation-{index}-candidate-off',n,True,oracle))
            check_global_census(reports)
            p=root/'J5-n83-validation-6-candidate-off/transform279-production-rank-42.txt'
            row=list(map(int,p.read_text().splitlines()[24].split()))
            self.assertEqual(row[5:8],[0,0,573]);self.assertGreater(row[4],0)
