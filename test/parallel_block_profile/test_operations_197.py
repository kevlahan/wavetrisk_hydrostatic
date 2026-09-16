import copy,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_workflow_195 import collect as collect_workflow
from analyze_operations_197 import collect,GROUPS
from cluster_197 import schedule,observer_mode,operation_mode
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
   call op_enter(78)
   call op_count(28,128_int64)
   call op_leave(78)
   call op_dependency(1,1,1,5)
   call op_post(2,100,100,197)
   call op_wait_begin
   call op_wait_end
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

class Operations197(unittest.TestCase):
 def test_schedule(self):
  jobs=schedule(2,True)
  self.assertEqual(len(jobs),14)
  self.assertEqual([j[0] for j in jobs if j[4]=='timing-bare'],['baseline-1','candidate-1','legacy-1','legacy-2','candidate-2','baseline-2'])
  self.assertEqual(len([j for j in jobs if j[4]=='focused-profile']),3)
  self.assertTrue(all(not j[3] for j in jobs))
  self.assertEqual(sorted(i for v in GROUPS.values() for i in v),list(range(83)))

 def test_reject_oracle_only_divergence_counts(self):
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   (root/'dependencies197-rank-0.txt').write_text('schema 1 0\n1 1 1 1 1 5 5 2 100 100 197 1 0.01 0.02\n')
   contexts={c:{'calls_global':[0]*83,'counters_global':[0]*35} for c in ('WT','outside_WT')}
   contexts['WT']['calls_global'][19]=1
   contexts['outside_WT']['calls_global'][81]=1
   work=contexts['outside_WT']['counters_global'];work[34]=10
   fake={'steps':[{'restart':0,'contexts':contexts}]}
   with patch('analyze_operations_197._base_collect',return_value=fake):
    with self.assertRaisesRegex(ValueError,'Production divergence sample counter missing'):collect(root,1,2,'block')
    work[32]=10;work[33]=11
    with self.assertRaisesRegex(ValueError,'evaluations exceed samples'):collect(root,1,2,'block')
    work[33]=8
    self.assertEqual(collect(root,1,2,'block')['ordinary_counters_per_rank_step']['outside_WT'][33],8)

 @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
 def test_compiled_observer_and_rejected_corruption(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);driver=root/'probe.f90';driver.write_text(DRIVER)
   source=Path(__file__).with_name('workflow_profile_197.f90')
   for opt in ('-O0','-O2'):
    subprocess.run(['gfortran',opt,'-fcheck=all','-Wall','-Wextra','-Werror','-ffp-contract=off',str(source),str(driver),'-o',str(root/'probe')],cwd=root,check=True,capture_output=True)
    for mode in (0,1,2):
     d=root/(opt+str(mode));d.mkdir()
     env={**os.environ,'WAVETRISK_WORKFLOW195':'1','WAVETRISK_OPERATIONS197':str(mode),'WAVETRISK_PROFILE_BLOCK_DETAIL':'0','WAVETRISK_PROFILE_PARALLEL_BLOCKS':'0'}
     subprocess.run([str(root/'probe')],cwd=d,env=env,check=True)
     collect_workflow(d,1,1);r=collect(d,1,mode);self.assertEqual(len(r['steps']),2)
     self.assertEqual(r['steps'][1]['restart'],1)
     c=r['steps'][0]['contexts']['outside_WT']['counters_global'];self.assertEqual(c[8],100 if mode else 0)
     self.assertEqual(r['steps'][0]['contexts']['WT']['counters_global'][22],42 if mode else 0)
     f=d/'operations197-rank-0.txt';raw=f.read_text()
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
     dep=d/'dependencies197-rank-0.txt';saved=dep.read_text()
     if mode:
      lines=saved.splitlines();fields=lines[1].split();fields[7]='-1';lines[1]=' '.join(fields)
      dep.write_text('\n'.join(lines)+'\n')
      with self.assertRaises(ValueError):collect(d,1,mode)
      dep.write_text(saved)
      self.assertEqual(r['dependencies']['records'],2)

if __name__=='__main__':unittest.main()
