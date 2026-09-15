import copy,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_workflow_195 import collect as collect_workflow
from analyze_operations_196 import collect,GROUPS
from cluster_196 import schedule,observer_mode,operation_mode
DRIVER='''program probe
use workflow_profile_195_mod
use iso_fortran_env, only: int64
implicit none
integer :: s,p,i
real :: x
call wf_init
do s=1,2
 call wf_begin
 do p=1,5
  call wf_phase(p)
  if(p==1) then
   call op_enter(2)
   call wf_enter(2)
   call op_enter(26)
   call op_enter(28)
   call op_enter(58)
   x=0
   do i=1,10000
    x=x+sin(real(i))
   end do
   if(x>huge(x)) stop 1
   call op_count(9,100_int64)
   call op_leave(58)
   call op_enter(38)
   call op_leave(38)
   call op_leave(28)
   call op_leave(26)
   call wf_leave(2)
   call wf_enter(1)
   call op_enter(65)
   call op_enter(70)
   call op_enter(71)
   call op_count(23,42_int64)
   call op_enter(72)
   call op_leave(72)
   call op_leave(71)
   call op_leave(70)
   call op_enter(73)
   call op_leave(73)
   call op_enter(74)
   call op_enter(75)
   call op_leave(75)
   call op_leave(74)
   call op_enter(76)
   call op_leave(76)
   call op_enter(77)
   call op_leave(77)
   call op_enter(16)
   call op_volume(1,2_int64,80_int64)
   call op_leave(16)
   call op_enter(66)
   call op_leave(66)
   call op_leave(65)
   call wf_leave(1)
   call op_leave(2)
  end if
  if(p==4) call wf_flag(2)
  if(p==5.and.s==2) call wf_flag(1)
 end do
 call wf_end
end do
call wf_finish(0)
end program
'''

class Operations196(unittest.TestCase):
 def test_schedule(self):
  jobs=schedule(6,True)
  self.assertEqual(len(jobs),58)
  self.assertTrue(all(not j[3] for j in jobs))
  for mode,m in [('bare',0),('minimal',0),('coarse',1),('detail',2)]:
   rows=[j for j in jobs if j[4]=='timing-'+mode];self.assertEqual(len(rows),12)
   self.assertTrue(all(observer_mode(j[0],j[4])==1 and operation_mode(j[0],j[4])==m for j in rows))
  self.assertEqual(sorted(i for v in GROUPS.values() for i in v),list(range(78)))

 @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
 def test_compiled_observer_and_rejected_corruption(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);driver=root/'probe.f90';driver.write_text(DRIVER)
   source=Path(__file__).with_name('workflow_profile_196.f90')
   for opt in ('-O0','-O2'):
    subprocess.run(['gfortran',opt,'-fcheck=all','-Wall','-Wextra','-Werror','-ffp-contract=off',str(source),str(driver),'-o',str(root/'probe')],cwd=root,check=True,capture_output=True)
    for mode in (0,1,2):
     d=root/(opt+str(mode));d.mkdir()
     env={**os.environ,'WAVETRISK_WORKFLOW195':'1','WAVETRISK_OPERATIONS196':str(mode),'WAVETRISK_PROFILE_BLOCK_DETAIL':'0','WAVETRISK_PROFILE_PARALLEL_BLOCKS':'0'}
     subprocess.run([str(root/'probe')],cwd=d,env=env,check=True)
     collect_workflow(d,1,1);r=collect(d,1,mode);self.assertEqual(len(r['steps']),2)
     self.assertEqual(r['steps'][1]['restart'],1)
     c=r['steps'][0]['contexts']['outside_WT']['counters_global'];self.assertEqual(c[8],100 if mode else 0)
     self.assertEqual(r['steps'][0]['contexts']['WT']['counters_global'][22],42 if mode else 0)
     f=d/'operations196-rank-0.txt';raw=f.read_text()
     with self.assertRaises(ValueError):collect(d,2,mode)
     with self.assertRaises(ValueError):collect(d,1,(mode+1)%3)
     f.write_text(raw+'schema 1 0 0\n')
     with self.assertRaises(ValueError):collect(d,1,mode)
     f.write_text('\n'.join(l for l in raw.splitlines() if not l.startswith('self 1 1 16 '))+'\n')
     with self.assertRaises(ValueError):collect(d,1,mode)
     lines=raw.splitlines()
     idx=next(i for i,l in enumerate(lines) if l.startswith('self 1 1 16 '));fields=lines[idx].split();fields[-1]='10.0';lines[idx]=' '.join(fields);f.write_text('\n'.join(lines)+'\n')
     with self.assertRaises(ValueError):collect(d,1,mode)
     f.write_text(raw)

if __name__=='__main__':unittest.main()
