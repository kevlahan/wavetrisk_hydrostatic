"""Compile the actual production storage implementation with bounds/FP checks."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
DRIVER='''program test
  use kind_mod, only : dp
  use parallel_block_scalar_storage_mod
  implicit none
  type(Scalar_Record_Storage) :: a,b
  integer :: span,groups,nf,nk,ns,s,f,node,zone,j,at,h
  integer,parameter :: shared(33)=[7,8,12,13,14,21,22,23,24,25,26,27,28,29,30,31,32,33, &
       36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]
  real(dp) :: record(50),geom(33,16)
  logical :: rebuilt
  nk=6;ns=2;nf=nk*ns
  do groups=1,2
    do span=16,48,16
      call scalar_allocate(a,groups*span*nf,span,nk,ns,-2,3,.true.,rebuilt)
      call scalar_allocate(b,groups*span*nf,span,nk,ns,-2,3,.false.,rebuilt)
      call scalar_fill(a,-999.0_dp);call scalar_fill(b,-999.0_dp)
      do s=0,groups*span*nf-1
        f=mod(s/span,nf);node=(s/(span*nf))*span+mod(s,span)+1
        zone=2
        if(mod(f,nk)-2>=1)zone=1
        do j=1,50
          record(j)=real(s*1000+j,dp)
        end do
        record(shared)=real(100*zone+10*node,dp)+real(shared,dp)
        call scalar_write_range(a,50*s+1,50*(s+1),record)
        call scalar_write_range(b,50*s+1,50*(s+1),record)
      end do
      do s=0,groups*span*nf-1
        if(any(scalar_read_range(a,50*s+1,50*(s+1))/=scalar_read_range(b,50*s+1,50*(s+1)))) &
             error stop 'compact record differs from independently expanded record'
      end do
      if(scalar_capacity(a)>=scalar_capacity(b)) error stop 'storage did not shrink'
      call scalar_allocate(a,groups*span*nf,span,nk,ns,-2,3,.true.,rebuilt)
      if(rebuilt) error stop 'unchanged allocation rebuilt'
    end do
  end do
  call scalar_allocate(a,2*16*nf,16,nk,ns,-2,3,.true.,rebuilt)
  geom=7.0_dp
  call scalar_seed_patch(a,1,geom)
  geom=8.0_dp
  call scalar_seed_patch(a,1+16*nf,geom)
  do s=0,2*16*nf-1
    f=mod(s/16,nf);h=s/(16*nf)
    record=scalar_read_range(a,50*s+1,50*(s+1))
    if(mod(f,nk)-2>=1)then
      if(any(record(shared)/=real(7+h,dp)))error stop 'physical geometry missing'
    else
      if(any(record(shared)/=0.0_dp))error stop 'inactive geometry activated'
    end if
    record(shared)=0.0_dp
    if(any(record/=0.0_dp))error stop 'field initialization differs'
  end do
  at=50*(3*16)+9
  call scalar_write(a,at,123.0_dp)
  if(scalar_read(a,at)/=123.0_dp)error stop 'field write missing'
  if(scalar_read(a,at+50)/=0.0_dp)error stop 'field write aliased'
  call scalar_write(a,[at,at+1,at+2,at+25],[1.0_dp,2.0_dp,3.0_dp,4.0_dp])
  if(any(scalar_read(a,[at,at+1,at+2,at+25])/=[1.0_dp,2.0_dp,3.0_dp,4.0_dp])) &
       error stop 'indexed field transfer differs'
  call scalar_share_inactive(a)
  if(any(scalar_read(a,shared)/=7.0_dp))error stop 'remote inactive geometry differs'
  geom=11.0_dp
  call scalar_install_geometry(a,17,geom)
  if(any(scalar_read(a,50*16*nf+shared)/=11.0_dp))error stop 'ghost geometry differs'
  call scalar_fill_fields(a,1,16,0.0_dp)
  if(any(scalar_read(a,shared)/=7.0_dp))error stop 'field fill erased geometry'
  call scalar_allocate(a,0,0,nk,ns,-2,3,.true.,rebuilt)
  if(scalar_extent(a)/=0.or.scalar_capacity(a)/=0)error stop 'empty storage invalid'
  call scalar_release(a)
  if(scalar_is_allocated(a))error stop 'release failed'
  print *, 'PASS scalar storage'
end program test
'''


class ScalarStorageTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
    def test_real_module_checked_layout_and_capacity(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'driver.f90').write_text(DRIVER)
            build=subprocess.run(['gfortran','-std=f2018','-Wall','-Wextra','-Werror=line-truncation',
                '-fcheck=all','-ffpe-trap=invalid,zero,overflow','-finit-real=snan',
                str(ROOT/'src/kind.f90'),str(ROOT/'src/parallel_block_scalar_storage.f90'),
                'driver.f90','-o','test'],cwd=out,capture_output=True,text=True)
            self.assertEqual(build.returncode,0,build.stdout+build.stderr)
            run=subprocess.run(['./test'],cwd=out,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertIn('PASS scalar storage',run.stdout)


if __name__=='__main__':unittest.main()
