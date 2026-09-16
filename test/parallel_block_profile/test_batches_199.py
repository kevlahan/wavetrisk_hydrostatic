import os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_batches_199 import collect,campaign_summary
from test_operations_199 import DRIVER

class Batches199(unittest.TestCase):
 def test_comparison_includes_steps_with_no_batch(self):
  import json
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);folder=root/'candidate-work-1';folder.mkdir()
   row=dict(step=1,context=4,operation=120,layer=0,level=0,position=0,cpu_mean=0.02,mpi_mean=0.,cpu_max_minus_mean=0.01)
   (folder/'batches-summary.json').write_text(json.dumps(dict(rows=[row],ordinary_steps=[1,2])))
   r=campaign_summary(root,{'candidate-work-1':dict(phase='diagnostic-work',side='block')})
   self.assertEqual(r['rows'][0]['cpu_excess'],0.01)

 @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
 def test_compiled_nested_batches_census_and_corruption(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);src=Path(__file__).parent
   body=DRIVER.replace(' call op_enter(26)', ''' call batch_begin199(101,3,6,work=10)
 call batch_begin199(102,3,6,work=2)
 if(count199) leaf199(36)=leaf199(36)+5_int64
 call batch_end199(102)
 call batch_begin199(105,3,6,work=10)
 call op_enter(39)
 call op_leave(39)
 call batch_end199(105)
 call batch_end199(101)
 call batch_begin199(109,3,6,work=7)
 call batch_end199(109)
 call op_enter(26)''')
   (root/'driver.f90').write_text(body)
   subprocess.run(['cc','-O2','-c',str(src/'workflow_clock_199.c'),'-o',str(root/'clock.o')],check=True)
   subprocess.run(['gfortran','-O0','-ffp-contract=off','-fcheck=all','-Wall','-Wextra','-Werror',str(src/'workflow_profile_199.f90'),str(root/'driver.f90'),str(root/'clock.o'),'-o',str(root/'probe')],cwd=root,check=True)
   for mode in range(4):
    out=root/str(mode);out.mkdir()
    env={**os.environ,'WAVETRISK_WORKFLOW195':'1','WAVETRISK_OPERATIONS199':str(mode),'WAVETRISK_PROFILE_BLOCK_DETAIL':'0','WAVETRISK_PROFILE_PARALLEL_BLOCKS':'0'}
    subprocess.run([str(root/'probe')],cwd=out,env=env,check=True)
    r=collect(out,1,mode)
    self.assertEqual(r['records'],6 if mode>=2 else 0)
    if mode<2:continue
    self.assertEqual(r['preceding_producers'][0]['producer_operation'],101)
    path=out/'batches199-rank-0.txt';original=path.read_text();rows=original.splitlines()
    f=rows[2].split();self.assertEqual(int(f[3]),1)
    self.assertEqual(int(rows[1].split()[10]),5)
    # Reject both a broken parent and an interval exceeding the full step.
    f[3]='999';rows[2]=' '.join(f);path.write_text('\n'.join(rows)+'\n')
    with self.assertRaises(ValueError):collect(out,1,mode)
    rows=original.splitlines();f=rows[1].split();f[-2]='100';rows[1]=' '.join(f);path.write_text('\n'.join(rows)+'\n')
    with self.assertRaises(ValueError):collect(out,1,mode)
    path.write_text(original)
    with self.assertRaises(ValueError):collect(out,2,mode)

if __name__=='__main__':unittest.main()
