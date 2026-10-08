"""Real Fortran clock-stack checks, serial compilation and no MPI."""
import json,os,subprocess,tempfile,unittest
from pathlib import Path
from cost_analysis import read_costs
from common import sha,save
HERE=Path(__file__).resolve().parent
PROGRAM='''program check_observer
 use cost_observer_mod
 implicit none
 integer :: s,k,j,mode
 real,volatile :: value
 character(20) :: arg
 call get_command_argument(1,arg)
 read(arg,*) mode
 call cost_enter(C_RK)
 do s=1,12
  call cost_begin_step
  call cost_enter(C_RK)
  if(mode==1) call cost_leave(C_WT)
  if(mode==2) then
   do j=1,40
    call cost_enter(C_WT)
   end do
  end if
  do k=1,4
   call cost_enter(C_TENDENCY)
   value=0.
   do j=1,10000
    value=value+sqrt(real(j))
   end do
   call cost_enter(C_WT)
   value=value+1.
   call cost_leave(C_WT)
   call cost_leave(C_TENDENCY)
  end do
  call cost_leave(C_RK)
  call cost_end_step
 end do
 if(mode==3) then
  do s=1,53
   call cost_begin_step
   call cost_end_step
  end do
 end if
 call cost_flush(0)
end program
'''

class Observer(unittest.TestCase):
 def test_clock_stack_and_disabled_and_invalid_states(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);(root/'test.f90').write_text(PROGRAM)
   for opt in ('-O0','-O2'):
    build=root/opt[1:];build.mkdir()
    subprocess.run(['gfortran',opt,'-fcheck=all',str(HERE/'cost_observer.f90'),str(root/'test.f90'),'-o','observer-test'],cwd=build,check=True,capture_output=True)
    for selector,mode,success in (('1',0,True),('0',0,True),('invalid',0,False),('1',1,False),('1',2,False),('1',3,False)):
     case=build/(selector+'-'+str(mode));case.mkdir()
     env=dict(os.environ,WAVETRISK_COST_AUDIT=selector)
     result=subprocess.run([str(build/'observer-test'),str(mode)],cwd=case,env=env,capture_output=True,text=True)
     self.assertEqual(result.returncode==0,success,result.stderr)
     if success and selector=='1':
      costs=read_costs(case,1,False)
      self.assertEqual(costs['phases']['tendency context']['calls'],4)
      self.assertGreater(costs['mean_rank_step_seconds'],0.)
     else:self.assertEqual(list(case.glob('cost-rank-*.txt')),[])

if __name__=='__main__':unittest.main()
