"""Exercise actual Stage202 kernels against their retained pre-prototype paths."""
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from test_scalar_storage import DRIVER
ROOT=Path(__file__).resolve().parents[2]
RECORDS='''program test
 use kind_mod
 use parallel_block_scalar_storage_mod
 implicit none
 type(Scalar_Record_Storage) :: a
 logical :: rebuilt
 integer :: start,n,i,mode,extent
 real(dp),allocatable :: x(:)
 open(10,file='result.bin',access='stream',form='unformatted',status='replace')
 do mode=0,1
 call scalar_allocate(a,2*3*6*2,3,6,2,-2,3,mode==1,rebuilt)
 extent=scalar_extent(a)
 call scalar_fill(a,-0.0_dp)
 do start=1,extent,7
   n=min(106,extent-start+1)
   allocate(x(n))
   do i=1,n
     x(i)=real(10000*start+i,dp)
     if(mod(i,5)==0)x(i)=-0.0_dp
   end do
   call scalar_write_range(a,start,start+n-1,x)
   ! Snapshot after each operation, including shared geometry alias effects.
   write(10)scalar_read_range(a,1,extent)
   do i=1,extent
     write(10)scalar_read(a,i)
   end do
   call scalar_write_range(a,start,start+n-1,-0.0_dp)
   write(10)scalar_read_range(a,1,extent)
   call scalar_write_range(a,start,start-1,x(:0))
   write(10)scalar_read_range(a,start,start-1)
   deallocate(x)
 end do
 end do
 close(10)
end program
'''

def alias_driver():
    s=(ROOT/'src/parallel_block_inverse.f90').read_text()
    routines=s[s.index('  subroutine apply_aliases('):s.index('  subroutine compare_native_inverse(')]
    return '''module harness
 use kind_mod
 implicit none
 integer,parameter :: nscalar=2,EDGE=3,zlevels=4
 integer,parameter :: DP_INVERSE_PACK=1,DP_INVERSE_LOCAL=2,DP_INVERSE_INSTALL=3
 type Alias_Plan
 integer,allocatable :: source(:),dest(:),sc(:),sd(:),rc(:),rd(:),install_order(:)
 end type
 type local_plan
 integer,allocatable :: source(:),dest(:)
 end type
 type(local_plan) :: local_alias(2)
 real(dp) :: value(5,4,8,2),alias_send(128),alias_recv(128)
 logical :: timing=.false.
 contains
 subroutine detail_enter(n)
 integer :: n
 end subroutine
 subroutine detail_leave(n)
 integer :: n
 end subroutine
 subroutine record_traffic(phase,sc,n)
 integer :: phase,sc(:),n
 end subroutine
 subroutine begin_exchange_values(send,sc,sd,recv,rc,rd,n,tag,nrequest)
 real(dp) :: send(:),recv(:)
 integer :: sc(:),sd(:),rc(:),rd(:),n,tag,nrequest
 recv=send
 nrequest=0
 end subroutine
 subroutine finish_exchange_values(n)
 integer :: n
 end subroutine
'''+routines+'''
end module
program test
 use harness
 implicit none
 type(Alias_Plan) :: a
 integer :: component,family,i,j,k,trial
 open(10,file='result.bin',access='stream',form='unformatted',status='replace')
 a%sc=[3];a%sd=[0];a%rc=[3];a%rd=[0]
 a%install_order=[3,1,2]
 do component=1,2
 local_alias(component)%source=[1,2,3,3]
 if(component==1)then
 a%source=[1,2,3];a%dest=[3,2,4]
 local_alias(component)%dest=[2,3,3,4]
 else
 a%source=[1,5,12];a%dest=[-5,6,-6]
 local_alias(component)%dest=[-2,3,-3,4]
 endif
 do trial=1,3
 do family=1,2
 do i=1,8
 do j=1,4
 do k=1,5
 value(k,j,i,family)=real(1000*family+100*i+10*j+k,dp)
 if(mod(k+j+i+trial,4)==0)value(k,j,i,family)=-0.0_dp
 enddo
 enddo
 enddo
 enddo
 do family=1,2
 call apply_aliases(a,component,family)
 write(10)value,alias_send,alias_recv
 enddo
 enddo
 enddo
 close(10)
end program
'''.replace('call apply_aliases(a,component,family)','alias_send=0;alias_recv=0\n call apply_aliases(a,component,family)')

class Batching202(unittest.TestCase):
    def build_run(self,source,flags,modules=(),output=False):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'test.f90').write_text(source)
            cmd=['gfortran','-cpp','-O2','-g','-fcheck=all','-ffp-contract=off','-ffpe-trap=invalid,zero,overflow',
                 '-ffree-line-length-132',*flags,str(ROOT/'src/kind.f90'),*[str(ROOT/n) for n in modules],'test.f90','-o','test']
            p=subprocess.run(cmd,cwd=root,capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            p=subprocess.run(['./test'],cwd=root,capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            return hashlib.sha256((root/'result.bin').read_bytes()).hexdigest() if output else p.stdout
    @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
    def test_record_layout_and_capacity(self):
        self.assertIn('PASS',self.build_run(DRIVER,['-DSTAGE202_RECORDS'],['src/parallel_block_scalar_storage.f90']))
    @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
    def test_mixed_ranges_canonical_geometry_and_signed_zero(self):
        modules=['src/parallel_block_scalar_storage.f90']
        self.assertEqual(self.build_run(RECORDS,[],modules,True),self.build_run(RECORDS,['-DSTAGE202_RECORDS'],modules,True))
    @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
    def test_snapshot_order_alias_chains_signed_destinations(self):
        source=alias_driver()
        self.assertEqual(self.build_run(source,[],output=True),self.build_run(source,['-DSTAGE202_INVERSE'],output=True))
if __name__=='__main__':unittest.main()
