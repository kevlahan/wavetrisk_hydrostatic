"""Exercise actual inner inverse kernel with geometry reuse on/off."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from test_scalar_inverse_193 import routine,ROOT


class InnerInverse194(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_kernel_and_geometry_lifecycle(self):
        text=(ROOT/'src/parallel_block_mpi.f90').read_text()
        kernel=routine(text,'reconstruct_block_inner_vector_level')
        setter=routine(text,'set_block_vector_sol_value')
        getter=re.search(r'^  real\(dp\) function block_vector_family_value\b.*?^  end function block_vector_family_value',text,re.M|re.S).group()
        source='''module actual
use iso_fortran_env, only: dp=>real64,int64
use ieee_arithmetic
implicit none
integer,parameter::PATCH_SIZE=4,N_CHDRN=4,zlevels=4,EDGE=3,RT=0,DG=1,UP=2,LORT=0,UPLT=1
integer,parameter::STORE_PATCH=1,STORE_BDRY=2,STORE_GHOST=3,BLOCK_PAYLOAD_SOL=1,BLOCK_PAYLOAD_WAV_COEFF=2
integer,parameter::BLOCK_SCALAR_EDGE_LENGTH_START=46,BLOCK_SCALAR_TRIANGLE_AREA_START=49
integer,parameter::chd_offs(2,4)=reshape([0,0,2,0,0,2,2,2],[2,4])
logical::stage194_inner_geometry=.false.,block_profile=.true.
integer(int64)::inner_geometry_work(4)=0,reads=0
real(dp)::geometry_epoch=0
character(16)::mode
type Patch
 integer::level=5,children(4)=0,elts_start
end type
type Block_Data
 type(Patch)::patch(5)
 integer::n_field_level=8,field_level=-2,vector_mult=3
 integer::node(80),bdry_node(64),ghost_node(48)
 real(dp)::vector(1920),wavelet_vector(1920),bdry_vector(1536),bdry_wavelet_vector(1536)
 real(dp)::ghost_vector(1152),ghost_wavelet_vector(1152)
end type
type Ready
 logical::ready=.true.
 integer::catalog_index=1
end type
type(Ready)::block_scalar_tendency(1)
type Block_Vector_Inverse_Context
 integer::target_level=4
 integer(int64)::cell_count=0,inner_edge_count=0
end type
contains
integer function catalog_local_block(i)
integer,intent(in)::i
catalog_local_block=i
end function
subroutine fail(message)
character(*),intent(in)::message
error stop message
end subroutine
subroutine locate_block_scalar_record(block,local_index,p,i,j,storage,record,node)
type(Block_Data),intent(in)::block
integer,intent(in)::local_index,p,i,j
integer,intent(out)::storage,record,node
record=p-1
if(i>=0.and.i<4.and.j>=0.and.j<4) then
 storage=STORE_PATCH;node=j*4+i
else if(i>=4.or.i<0) then
 storage=STORE_BDRY;node=modulo((j+2)*8+i+2,64)
else
 storage=STORE_GHOST;node=modulo((j+2)*8+i+2,48)
end if
end subroutine
real(dp) function block_scalar_record_value(block,local_index,p,s,l,i,j,slot) result(v)
type(Block_Data),intent(in)::block
integer,intent(in)::local_index,p,s,l,i,j,slot
reads=reads+1
if(block%field_level+l-1<1.or.block%field_level+l-1>zlevels) error stop 'read soil geometry'
v=real(100+p*3+i+j*2+slot,dp)/10003.0_dp+geometry_epoch
if(mode=='bad-length'.and.slot==46) v=0
if(mode=='bad-area'.and.slot==49) v=0
end function
'''+kernel+'\n'+getter+'\n'+setter+'''
end module
program check
use actual
implicit none
type(Block_Data)::a,b,seed
type(Block_Vector_Inverse_Context)::ca,cb
integer::n,p,trial,epoch,layers,levels(4)=[-2,1,3,5]
integer(int64)::old_reads,new_reads
call get_command_argument(1,mode)
do p=1,5
 seed%patch(p)%elts_start=(p-1)*16
end do
seed%patch(1)%level=4;seed%patch(1)%children=[1,2,3,4]
do n=1,size(seed%vector)
 seed%vector(n)=sin(real(n,dp));seed%wavelet_vector(n)=cos(real(n,dp))/7
end do
do n=1,size(seed%bdry_vector)
 seed%bdry_vector(n)=sin(real(n,dp)/3);seed%bdry_wavelet_vector(n)=cos(real(n,dp))/11
end do
do n=1,size(seed%ghost_vector)
 seed%ghost_vector(n)=sin(real(n,dp)/5);seed%ghost_wavelet_vector(n)=cos(real(n,dp))/13
end do
do epoch=1,2
 geometry_epoch=real(epoch-1,dp)/100
 do trial=1,4
  seed%field_level=levels(trial)
  layers=max(0,min(seed%n_field_level,zlevels-seed%field_level+1)-max(1,2-seed%field_level)+1)
  a=seed;b=seed;ca=Block_Vector_Inverse_Context();cb=ca
  stage194_inner_geometry=.false.;reads=0;inner_geometry_work=0
  if(mode=='') call reconstruct_block_inner_vector_level(1,a,ca)
  old_reads=reads
  stage194_inner_geometry=.true.;reads=0;inner_geometry_work=0
  call reconstruct_block_inner_vector_level(1,b,cb)
  new_reads=reads
  if(any(a%vector/=b%vector)) error stop 'cached vector differs'
  if(any(a%bdry_vector/=b%bdry_vector).or.any(a%ghost_vector/=b%ghost_vector)) error stop 'alias overwritten'
  if(any(b%wavelet_vector/=seed%wavelet_vector)) error stop 'wavelet overwritten'
  if(ca%cell_count/=cb%cell_count.or.ca%inner_edge_count/=cb%inner_edge_count) error stop 'coverage differs'
  if(layers>0) then
   if(new_reads*layers/=old_reads.or.old_reads/=32*ca%cell_count) error stop 'geometry work differs'
   if(inner_geometry_work(3)/=new_reads.or.inner_geometry_work(4)/=old_reads) error stop 'counters differ'
  else
   if(new_reads/=0) error stop 'empty column read geometry'
  end if
 end do
end do
print *, 'PASS exact values, two geometry epochs, partial/empty layers and geometry work'
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'kernel.f90').write_text(source)
            for opt in ('-O0','-O2'):
                result=subprocess.run(['gfortran',opt,'-fcheck=all','-ffp-contract=off','-ffree-line-length-132',
                    'kernel.f90','-o','kernel'],cwd=root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                result=subprocess.run([str(root/'kernel')],cwd=root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                for mode in ('bad-length','bad-area'):
                    result=subprocess.run([str(root/'kernel'),mode],cwd=root,capture_output=True,text=True)
                    self.assertNotEqual(result.returncode,0)
                    self.assertIn('nonpositive',result.stderr)


    @unittest.skipUnless(os.environ.get('WAVETRISK_RUN_MPI_TESTS')=='1','Set WAVETRISK_RUN_MPI_TESTS=1')
    def test_actual_selector_defaults_controls_and_disagreement(self):
        text=(ROOT/'src/parallel_block_mpi.f90').read_text()
        declaration=next(line for line in text.splitlines() if 'logical, save :: stage194_ready' in line)
        source='''module actual
use mpi_f08
integer::rank
type(MPI_Comm)::comm=MPI_COMM_WORLD
'''+declaration+'\ncontains\n'+routine(text,'configure_stage194')+'''
subroutine fail(message)
character(*),intent(in)::message
integer::ierr
call MPI_Abort(comm,1,ierr)
error stop message
end subroutine
subroutine check_mpi(ierr,message)
integer,intent(in)::ierr
character(*),intent(in)::message
if(ierr/=MPI_SUCCESS)call fail(message)
end subroutine
end module
program check_options
use actual
integer::ierr
character(16)::arg
call MPI_Init(ierr)
call MPI_Comm_rank(comm,rank,ierr)
call get_command_argument(1,arg)
if(arg=='disagree') stage194_inner_geometry=rank==0
call configure_stage194
call configure_stage194
call MPI_Finalize(ierr)
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'options.f90').write_text(source)
            subprocess.run(['mpif90','-O2','-fcheck=all','options.f90','-o','options'],cwd=root,check=True,capture_output=True)
            clean={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
            for setting,arg,success in [(None,'',True),('0','',True),('1','',True),('invalid','',False),(None,'disagree',False)]:
                env=dict(clean)
                if setting is not None:env['WAVETRISK_STAGE194_INNER_GEOMETRY']=setting
                result=subprocess.run(['mpirun','-n','2',str(root/'options'),arg],cwd=root,env=env,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                if success:self.assertEqual(result.stdout.count('Stage 194 requested inner inverse geometry = '+(setting or '0')),1)



if __name__=='__main__':unittest.main()
