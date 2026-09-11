"""Check batched field writes against independent per-record oracle storage."""
from pathlib import Path
import shutil,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
DRIVER='''program test
use kind_mod, only: dp
use parallel_block_scalar_storage_mod
implicit none
type(Scalar_Record_Storage)::a,b,c
real(dp),allocatable::value(:,:)
integer::span,groups,nk,ns,g,f,q,s,j
logical::rebuilt
character(20)::mode
call get_command_argument(1,mode)
nk=6;ns=2
if(len_trim(mode)>0)then
 call scalar_allocate(a,16*nk*ns,16,nk,ns,-2,3,.true.,rebuilt)
 allocate(value(3,16));value=1.0_dp
 select case(trim(mode))
 case('extent')
  call scalar_write_field_records(a,16*nk*ns,9,value)
 case('boundary')
  call scalar_write_field_records(a,1,9,value)
 case('geometry')
  call scalar_write_field_records(a,0,7,value)
 case('mixed')
  call scalar_write_field_records(a,0,5,value)
 case('slots')
  call scalar_write_field_records(a,0,49,value)
 end select
 error stop 'invalid batch accepted'
end if
do groups=1,2
 do span=16,48,16
  call scalar_allocate(a,groups*span*nk*ns,span,nk,ns,-2,3,.true.,rebuilt)
  call scalar_allocate(b,groups*span*nk*ns,span,nk,ns,-2,3,.false.,rebuilt)
  call scalar_allocate(c,groups*span*nk*ns,span,nk,ns,-2,3,.false.,rebuilt)
  call scalar_fill(a,-991.0_dp);call scalar_fill(b,-991.0_dp);call scalar_fill(c,-991.0_dp)
  allocate(value(3,span))
  do g=0,groups-1
   do f=0,nk*ns-1
    s=(g*nk*ns+f)*span
    do q=1,span
     value(:,q)=[(real(10000*g+100*f+10*q+j,dp),j=1,3)]
    end do
    call scalar_write_field_records(a,s,9,value)
    call scalar_write_field_records(c,s,9,value)
    do q=1,span
     call scalar_write_range(b,50*(s+q-1)+9,50*(s+q-1)+11,value(:,q))
    end do
   end do
  end do
  do s=0,groups*span*nk*ns-1
   if(any(abs(scalar_read_range(a,50*s+1,50*s+50)-scalar_read_range(b,50*s+1,50*s+50))>0.0_dp)) &
        error stop 'compact batch changed reference or untouched fields'
   if(any(abs(scalar_read_range(c,50*s+1,50*s+50)-scalar_read_range(b,50*s+1,50*s+50))>0.0_dp)) &
        error stop 'full batch changed reference or untouched fields'
  end do
  deallocate(value)
 end do
end do
print *, 'PASS field batches'
end program
'''
class ScalarFieldBatch(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_record_equivalence_and_rejections(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'test.f90').write_text(DRIVER)
            for opt in ('-O0','-O2'):
                build=subprocess.run(['gfortran',opt,'-fcheck=all','-ffpe-trap=invalid,zero,overflow','-Werror=line-truncation',str(ROOT/'src/kind.f90'),str(ROOT/'src/parallel_block_scalar_storage.f90'),'test.f90','-o','test'],cwd=out,capture_output=True,text=True)
                self.assertEqual(build.returncode,0,build.stdout+build.stderr)
                run=subprocess.run(['./test'],cwd=out,capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                for mode,reason in [('extent','storage extent'),('boundary','crosses field boundary'),('geometry','requires contiguous field slots'),('mixed','requires contiguous field slots'),('slots','shape invalid')]:
                    run=subprocess.run(['./test',mode],cwd=out,capture_output=True,text=True)
                    self.assertNotEqual(run.returncode,0)
                    self.assertIn(reason,run.stderr)
if __name__=='__main__':unittest.main()
