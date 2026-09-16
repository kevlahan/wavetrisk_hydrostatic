"""Actual-module field batches against the scalar reference accessors."""
import shutil,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
DRIVER='''program test
 use kind_mod
 use parallel_block_scalar_storage_mod
 implicit none
 type(Scalar_Record_Storage) :: a,b
 logical :: rebuilt
 integer :: mode,i,j,first,n,slots(4),gs(3),r
 real(dp) :: pack(64),geom(3,16),x
 integer(kind=8) :: left,right
 slots=[9,10,11,35];gs=[12,13,14]
 do mode=0,1
 call scalar_allocate(a,192,16,6,2,-2,3,mode==1,rebuilt)
 call scalar_allocate(b,192,16,6,2,-2,3,mode==1,rebuilt)
 call scalar_fill(a,-0.0_dp);call scalar_fill(b,777.0_dp)
 do i=1,scalar_extent(a)
   x=real(i,dp)
   if(mod(i,7)==0)x=-0.0_dp
   call scalar_write(a,i,x)
 end do
 do first=1,192,16
   call scalar_pack_fields203(a,first,16,slots,pack)
   do i=1,16
     do j=1,4
       x=scalar_read(a,(first+i-2)*50+slots(j))
       left=transfer(x,left);right=transfer(pack((i-1)*4+j),right)
       if(left/=right) stop 1
     end do
   end do
   call scalar_unpack_fields203(b,first,16,slots,pack)
   do i=1,16
     do j=1,4
       x=scalar_read(b,(first+i-2)*50+slots(j))
       left=transfer(x,left);right=transfer(pack((i-1)*4+j),right)
       if(left/=right) stop 2
     end do
     if(scalar_read(b,(first+i-2)*50+1)/=777.0_dp) stop 3
   end do
   call scalar_geometry_patch203(a,first,16,gs,geom)
   do i=1,16
     do j=1,3
       x=scalar_read(a,(first+i-2)*50+gs(j))
       left=transfer(x,left);right=transfer(geom(j,i),right)
       if(left/=right) stop 4
     end do
   end do
 end do
 call scalar_fill(b,888.0_dp)
 call scalar_copy_fields203(a,17,b,49,16,slots)
 do i=0,15
   do j=1,4
     x=scalar_read(a,(16+i)*50+slots(j));left=transfer(x,left)
     x=scalar_read(b,(48+i)*50+slots(j));right=transfer(x,right)
     if(left/=right) stop 5
   end do
 end do
 call scalar_pack_fields203(a,193,0,slots,pack(:0))
 call scalar_unpack_fields203(b,193,0,slots,pack(:0))
 call scalar_release(a);call scalar_release(b)
 end do
 print *, 'PASS field batches, full/compact, geometry classes, signed zero, empty and untouched fields'
end program
'''
class Fields203(unittest.TestCase):
 def test_actual_fields_against_independent_accessors(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'test.f90').write_text(DRIVER)
   for opt in ('-O0','-O2'):
    subprocess.run(['gfortran','-cpp',opt,'-fcheck=all','-ffp-contract=off','-DSTAGE202_RECORDS','-J'+d,'-I'+d,str(ROOT/'src/kind.f90'),str(ROOT/'src/parallel_block_scalar_storage.f90'),str(p/'test.f90'),'-o',str(p/'test')],check=True)
    subprocess.run([str(p/'test')],check=True,capture_output=True)
if __name__=='__main__':unittest.main()
