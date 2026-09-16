"""Exercise actual sparse canonical-cache helpers with overlapping alias slots."""
from pathlib import Path
import subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
def source():
 s=(ROOT/'src/parallel_block_inverse.f90').read_text()
 parts=[]
 for name in ('store203','remember_outer203','install_inverse203'):
  a=s.index('  subroutine '+name);b=s.index('  end subroutine',a)+len('  end subroutine')
  parts.append(s[a:b])
 return '''module testmod
 use kind_mod
 use iso_fortran_env,only:int64
 implicit none
 integer,parameter :: nscalar=2,EDGE=3,zlevels=2,BLOCK_PAYLOAD_SOL=1
 logical,parameter :: ledger203=.false.,resident203=.true.
 type Transfer_Plan
 integer,allocatable :: slot(:)
 logical :: edgewise=.false.
 end type
 type(Transfer_Plan) :: interior_plan,boundary_plan
 real(dp) :: value(5,2,4,2),restore203(5,2,1),saved203(5,2,1)
 real(dp) :: installed(5,2,2),boundary(5,2,2)
 integer :: restore203_slot(1)=[1],restore203_index(4)=[1,0,0,0],calls=0
 integer(int64) :: work203(6)=0
 contains
 subroutine die(message)
 character(*),intent(in)::message
 print *,message
 stop 99
 end subroutine
 subroutine transfer_nodes(plan,c,f,install)
 type(Transfer_Plan),intent(in)::plan
 integer,intent(in)::c,f
 logical,intent(in)::install
 integer::lo,hi
 if(.not.install) stop 98
 lo=merge(1,3,c==1);hi=merge(2,5,c==1)
 calls=calls+1
 if(calls<=2) then
 installed(lo:hi,:,:)=value(lo:hi,:,plan%slot,f)
 else
 boundary(lo:hi,:,:)=value(lo:hi,:,plan%slot,f)
 end if
 end subroutine
'''+ '\n'.join(parts)+'''
end module
program test
 use testmod
 implicit none
 integer(int64)::bits,zero
 interior_plan%slot=[1,2];boundary_plan%slot=[1,3]
 value=999.0_dp;value(:,:,2,1)=222.0_dp;restore203=111.0_dp
 call store203(1,1,1,-0.0_dp)
 value(:,:,1,1)=999.0_dp
 value(3:5,:,1,1)=333.0_dp
 call remember_outer203(interior_plan)
 call install_inverse203
 zero=transfer(-0.0_dp,zero);bits=transfer(installed(1,1,1),bits)
 if(bits/=zero) stop 1
 if(installed(2,1,1)/=111.0_dp) stop 2
 if(any(installed(3:5,:,1)/=333.0_dp)) stop 3
 if(any(installed(:,:,2)/=222.0_dp)) stop 4
 if(any(boundary(1:2,:,1)/=999.0_dp)) stop 5
 if(any(boundary(3:5,:,1)/=333.0_dp)) stop 6
 if(any(value(1:2,:,1,1)/=999.0_dp)) stop 7
 if(calls/=4) stop 8
end program
'''
class Workspace203(unittest.TestCase):
 def test_canonical_and_boundary_epochs_remain_distinct(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'test.f90').write_text(source())
   for opt in ('-O0','-O2'):
    subprocess.run(['gfortran',opt,'-fcheck=all','-J'+d,'-I'+d,str(ROOT/'src/kind.f90'),str(p/'test.f90'),'-o',str(p/'test')],check=True)
    subprocess.run([str(p/'test')],check=True)
if __name__=='__main__':unittest.main()
