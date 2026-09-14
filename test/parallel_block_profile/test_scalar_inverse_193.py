"""Compare the actual scalar and column kernels, including mixed storage and soil."""
from pathlib import Path
import re
import os
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]


def routine(text,name,kind='subroutine'):
    return re.search(r'^  '+kind+' '+name+r'\b.*?^  end '+kind+' '+name+r'\b',text,re.M|re.S).group()


class ScalarInverse193(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_kernels_match_across_fields_and_storage(self):
        text=(ROOT/'src/parallel_block_mpi.f90').read_text()
        routines='\n'.join(routine(text,n) for n in ('reconstruct_block_scalar_level',
            'reconstruct_block_scalar_columns','set_block_scalar_sol_value'))
        family=re.search(r'^  real\(dp\) function block_scalar_family_value\b.*?^  end function block_scalar_family_value',text,re.M|re.S).group()
        header='''module actual
use iso_fortran_env, only: dp=>real64,int64,error_unit
use ieee_arithmetic
implicit none
integer,parameter::PATCH_SIZE=4,N_CHDRN=4,zlevels=4,FROZEN=-1
integer,parameter::STORE_PATCH=1,STORE_BDRY=2,STORE_GHOST=3,BLOCK_PAYLOAD_SOL=1,BLOCK_PAYLOAD_WAV_COEFF=2
integer,parameter::BLOCK_SCALAR_AREA_INDEX=7,BLOCK_SCALAR_OVERLAP_START=27,BLOCK_SCALAR_WAVELET_MASK_INDEX=36
integer,parameter::chd_offs(2,4)=reshape([0,0,2,0,0,2,2,2],[2,4])
integer::rank=0
logical::stage193_scalar_columns=.false.,block_profile=.true.
integer(int64)::scalar_inverse_work(6)=0,lookups=0,geometry=0
type Patch
 integer::level,children(4)=0,elts_start
end type
type Block_Data
 type(Patch)::patch(5)
 integer::n_scalar_variable=3,n_field_level=8,field_level=-2
 integer::node(80),bdry_node(64),ghost_node(48)
 real(dp)::scalar(1920),wavelet_scalar(1920),bdry_scalar(1536),bdry_wavelet_scalar(1536)
 real(dp)::ghost_scalar(1152),ghost_wavelet_scalar(1152)
end type
type Ready
 logical::ready=.true.
 integer::catalog_index=1
end type
type(Ready)::block_scalar_tendency(1)
type Block_Scalar_Inverse_Context
 integer::target_level=4,phase=1
 integer(int64)::parent_patch_count=0,candidate_count=0,produced_count=0,frozen_count=0
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
subroutine locate_block_scalar_record(block,local_index,p,i,j,storage,record,node,remap)
type(Block_Data),intent(in)::block
integer,intent(in)::local_index,p,i,j
integer,intent(out)::storage,record,node
logical,optional,intent(in)::remap
lookups=lookups+1
record=p-1
if(i>=0.and.i<4.and.j>=0.and.j<4) then
 storage=STORE_PATCH;node=j*4+i
else if(i<0.or.j<0) then
 storage=STORE_BDRY;node=modulo((j+2)*8+i+2,64)
else
 storage=STORE_GHOST;node=modulo((j+2)*8+i+2,48)
end if
end subroutine
real(dp) function block_scalar_record_value(block,local_index,p,s,l,i,j,slot) result(v)
type(Block_Data),intent(in)::block
integer,intent(in)::local_index,p,s,l,i,j,slot
geometry=geometry+1
! All physical fields share geometry; nonphysical storage must never be read.
if(block%field_level+l-1<1.or.block%field_level+l-1>zlevels) error stop 'read soil geometry'
if(slot==BLOCK_SCALAR_WAVELET_MASK_INDEX) then
 v=0
 if(p==2.and.i==0.and.j==0) v=FROZEN
else
 v=real(100+p*3+i+j*2+slot,dp)/10003.0_dp
end if
end function
'''
        driver='''
end module
program check
use actual
implicit none
type(Block_Data)::a,b,seed
type(Block_Scalar_Inverse_Context)::ca,cb
integer::n,p,phase,trial,samples,levels(4)=[-2,1,3,5]
integer(int64)::old_lookups,old_geometry,new_lookups,new_geometry
character(16)::mode
call get_command_argument(1,mode)
do p=1,5
 seed%patch(p)%level=5;seed%patch(p)%elts_start=(p-1)*16
end do
seed%patch(1)%level=4;seed%patch(1)%children=[1,2,3,4]
do n=1,size(seed%scalar)
 seed%scalar(n)=sin(real(n,dp));seed%wavelet_scalar(n)=cos(real(n,dp))/7.0_dp
end do
do n=1,size(seed%bdry_scalar)
 seed%bdry_scalar(n)=sin(real(n,dp)/3);seed%bdry_wavelet_scalar(n)=cos(real(n,dp))/11
end do
do n=1,size(seed%ghost_scalar)
 seed%ghost_scalar(n)=sin(real(n,dp)/5);seed%ghost_wavelet_scalar(n)=cos(real(n,dp))/13
end do
if(mode=='poison') seed%ghost_wavelet_scalar=ieee_value(0.0_dp,ieee_quiet_nan)
do trial=1,4
 seed%field_level=levels(trial)
 samples=max(0,min(seed%n_field_level,zlevels-seed%field_level+1)-max(1,2-seed%field_level)+1)*seed%n_scalar_variable
 a=seed;b=seed
 do phase=1,2
  ca=Block_Scalar_Inverse_Context(phase=phase);cb=ca
  lookups=0;geometry=0
  if(mode/='poison') call reconstruct_block_scalar_level(1,a,ca)
  old_lookups=lookups;old_geometry=geometry
  lookups=0;geometry=0;scalar_inverse_work=0
  call reconstruct_block_scalar_columns(1,b,cb)
  new_lookups=lookups;new_geometry=geometry
  if(any(a%scalar/=b%scalar)) error stop 'scalar and column result differ'
  if(any(a%bdry_scalar/=b%bdry_scalar).or.any(a%ghost_scalar/=b%ghost_scalar)) error stop 'alias overwritten'
  if(any(b%wavelet_scalar/=seed%wavelet_scalar)) error stop 'wavelet overwritten'
  if(ca%candidate_count/=cb%candidate_count.or.ca%produced_count/=cb%produced_count.or. &
     ca%frozen_count/=cb%frozen_count.or.ca%parent_patch_count/=cb%parent_patch_count) error stop 'coverage differs'
  if(samples>0) then
   if(new_lookups*samples/=old_lookups.or.new_geometry*samples/=old_geometry) error stop 'hoisting factor differs'
   if(any(scalar_inverse_work(3:4)/=[new_lookups,new_geometry])) error stop 'measured work counters differ'
   if(any(scalar_inverse_work(5:6)/=[old_lookups,old_geometry])) error stop 'reference-equivalent counters differ'
  else
   if(new_lookups/=0.or.new_geometry/=0) error stop 'work on empty physical column'
  end if
 end do
end do
print *, 'PASS exact kernels, frozen nodes, three stores, soil, empty columns and work counts'
end program
'''
        source=header+routines+'\n'+family+'\n'+driver
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'kernel.F90').write_text(source)
            for opt in ('-O0','-O2'):
                result=subprocess.run(['gfortran',opt,'-cpp','-fcheck=all','-ffp-contract=off',
                    '-ffree-line-length-132','kernel.F90','-o','kernel'],cwd=root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                result=subprocess.run([str(root/'kernel')],cwd=root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                result=subprocess.run([str(root/'kernel'),'poison'],cwd=root,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('source field is non-finite',result.stderr)


    @unittest.skipUnless(os.environ.get('WAVETRISK_RUN_MPI_TESTS')=='1','Set WAVETRISK_RUN_MPI_TESTS=1')
    def test_actual_selector_defaults_controls_and_disagreement(self):
        text=(ROOT/'src/parallel_block_mpi.f90').read_text()
        declaration=next(line for line in text.splitlines() if 'logical, save :: stage193_ready' in line)
        source='''module actual
use mpi_f08
integer::rank
type(MPI_Comm)::comm=MPI_COMM_WORLD
'''+declaration+'\ncontains\n'+routine(text,'configure_stage193')+'''
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
if(arg=='disagree') stage193_scalar_columns=rank==0
call configure_stage193
call configure_stage193
call MPI_Finalize(ierr)
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'options.f90').write_text(source)
            subprocess.run(['mpif90','-O2','-fcheck=all','options.f90','-o','options'],cwd=root,check=True,capture_output=True)
            clean={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
            for setting,arg,success in [(None,'',True),('0','',True),('1','',True),('invalid','',False),(None,'disagree',False)]:
                env=dict(clean)
                if setting is not None:env['WAVETRISK_STAGE193_SCALAR_COLUMNS']=setting
                result=subprocess.run(['mpirun','-n','2',str(root/'options'),arg],cwd=root,env=env,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                if success:self.assertEqual(result.stdout.count('Stage 193 requested scalar inverse columns = '+(setting or '0')),1)


if __name__=='__main__': unittest.main()
