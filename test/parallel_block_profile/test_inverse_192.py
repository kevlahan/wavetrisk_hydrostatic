"""Exercise actual sparse transfer selection and rank-consistent runtime control."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from test_install_overlap_191 import sub, ROOT


class Inverse192(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_plan_keeps_writes_and_order(self):
        text=(ROOT/'src/parallel_block_inverse.f90').read_text()
        plan=re.search(r'  type :: Transfer_Plan\n.*?  end type',text,re.S).group()
        source='module actual\ninteger,parameter :: n_process=3\n'+plan+'\ncontains\n'
        source+=sub(text,'displacements')+'\n'+sub(text,'select_transfer')+'''
subroutine compile_local_block_inverse_routes(key,address)
integer,intent(in)::key(:,:)
integer,allocatable,intent(out)::address(:,:)
address=key
end subroutine
end module
program check
use actual
type(Transfer_Plan)::base,full,sparse
integer::i,n,at,peer,level,flags,j
! Every combination of alias/scaffold/outer flags at four levels. Different
! send and receive peer segmentation exercises both independently derived ends.
n=4*64
allocate(base%sc(3),base%sd(3),base%rc(3),base%rd(3),base%slot(n),base%key(4,n))
allocate(base%source_info(2,n),base%target_info(2,n))
base%sc=[64,128,64];base%rc=[128,0,128]
call displacements(base%sc,base%sd)
call displacements(base%rc,base%rd)
at=0
do level=3,6
 do flags=0,63
  at=at+1
  base%source_info(:,at)=[level,flags]
  base%target_info(:,at)=[level,flags]
  base%slot(at)=at
  base%key(:,at)=at
 end do
end do
call select_transfer(full,base,2,5,.true.,.true.)
call select_transfer(sparse,base,2,5,.true.,.true.,.true.)
! Independently specify the writes eligible at the fine/coarse epoch.
at=0
do i=1,n
 level=base%source_info(1,i);flags=base%source_info(2,i)
 if(level/=4.and.level/=5.and..not.btest(flags,4).and..not.btest(flags,2)) cycle
 if(.not.(btest(flags,1).or.btest(flags,2).or.btest(flags,5))) cycle
 at=at+1
 if(sparse%slot(at)/=i.or.any(sparse%key(:,at)/=i)) error stop 'write/order lost'
end do
if(at/=size(sparse%slot).or.at/=size(sparse%key,2)) error stop 'extra records'
if(size(sparse%slot)>=size(full%slot)) error stop 'no pruning'
do peer=1,3
 i=base%sd(peer)+1;j=base%sd(peer)+base%sc(peer)
 if(count(sparse%slot>=i.and.sparse%slot<=j)/=sparse%sc(peer)) error stop 'send counts'
 i=base%rd(peer)+1;j=base%rd(peer)+base%rc(peer)
 if(count(sparse%key(1,:)>=i.and.sparse%key(1,:)<=j)/=sparse%rc(peer)) error stop 'recv counts'
end do
call select_transfer(sparse,base,2,5,.true.,.true.,.false.)
if(any(sparse%slot/=full%slot).or.any(sparse%key/=full%key)) error stop 'control differs'
! Boundary path remains unpruned: only the caller's interior outer plan opts in.
print *, 'PASS sparse writes and control'
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'plan.f90').write_text(source)
            for opt in ('-O0','-O2'):
                subprocess.run(['gfortran',opt,'-fcheck=all','-ffree-line-length-132','plan.f90','-o','plan'],cwd=root,check=True,capture_output=True)
                subprocess.run([str(root/'plan')],check=True,capture_output=True)

    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_edge_scatter_preserves_other_edges(self):
        text=(ROOT/'src/parallel_block.f90').read_text()
        source="""module actual
use iso_fortran_env, only: dp=>real64
integer,parameter::EDGE=3,zlevels=2,STORE_PATCH=1,BLOCK_PAYLOAD_SOL=1,BLOCK_PAYLOAD_WAV_COEFF=2
type Block_Data
 integer::n_scalar_variable=2
 real(dp)::scalar(24),wavelet_scalar(24),vector(36),wavelet_vector(36)
 real(dp)::bdry_scalar(24),bdry_wavelet_scalar(24),bdry_vector(36),bdry_wavelet_vector(36)
end type
type(Block_Data)::block_local(1)
contains
"""+re.search(r'^subroutine transfer_local_block_inverse_routes\b.*?^end subroutine transfer_local_block_inverse_routes',text,re.M|re.S).group()+"""
end module
program check
use actual
integer::addresses(7,2),i,k,store,family
real(dp)::buffer(4),expected(36)
addresses(1:6,1)=[1,STORE_PATCH,1,4,1,3]
addresses(1:6,2)=[1,STORE_PATCH,2,4,1,3]
addresses(7,:)=[2,3]
buffer=[10.0_dp,20.0_dp,30.0_dp,40.0_dp]
expected=-1.0_dp
expected(17)=10.0_dp;expected(29)=20.0_dp
expected(21)=30.0_dp;expected(33)=40.0_dp
do store=1,2
 do family=1,2
  block_local(1)%vector=-1.0_dp
  block_local(1)%wavelet_vector=-1.0_dp
  block_local(1)%bdry_vector=-1.0_dp
  block_local(1)%bdry_wavelet_vector=-1.0_dp
  addresses(2,:)=store
  call transfer_local_block_inverse_routes(addresses(1:6,:),family,2,.true.,buffer,addresses(7,:))
  if(store==1.and.family==1) then
   if(any(block_local(1)%vector/=expected)) error stop 'interior solution'
  else if(store==1) then
   if(any(block_local(1)%wavelet_vector/=expected)) error stop 'interior wavelet'
  else if(family==1) then
   if(any(block_local(1)%bdry_vector/=expected)) error stop 'boundary solution'
  else
   if(any(block_local(1)%bdry_wavelet_vector/=expected)) error stop 'boundary wavelet'
  end if
 end do
end do
print *, 'PASS selected edge writes with atmospheric offset'
end program
"""
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'edges.f90').write_text(source)
            for opt in ('-O0','-O2'):
                subprocess.run(['gfortran',opt,'-fcheck=all','-ffree-line-length-132','edges.f90','-o','edges'],cwd=root,check=True,capture_output=True)
                subprocess.run([str(root/'edges')],check=True,capture_output=True)

    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_combined_edge_schedule(self):
        text=(ROOT/'src/parallel_block_inverse.f90').read_text()
        plan=re.search(r'  type :: Transfer_Plan\n.*?  end type',text,re.S).group()
        source='module actual\ninteger,parameter :: n_process=3,EDGE=3\n'+plan+'\ncontains\n'
        source+=sub(text,'displacements')+'\n'+sub(text,'combine_outer_edges')+"""
subroutine compile_local_block_inverse_routes(key,address)
integer,intent(in)::key(:,:)
integer,allocatable,intent(out)::address(:,:)
allocate(address(6,size(key,2)))
address=0
address(1:4,:)=key
end subroutine
end module
program check
use actual
type(Transfer_Plan)::a,b,p
integer::i
a%sc=[2,0,1];a%rc=[1,2,0];a%sd=[0,2,2];a%rd=[0,1,3]
a%slot=[1,2,3]
allocate(a%source_info(2,3),a%target_info(2,3),a%key(4,3))
a%source_info(1,:)=0;a%target_info(1,:)=0
a%source_info(2,:)=[64,256,192];a%target_info(2,:)=[256,64,192]
do i=1,3
 a%key(:,i)=10+i
end do
b%sc=[1,0,1];b%rc=[0,1,1];b%sd=[0,1,1];b%rd=[0,0,1]
b%slot=[4,5]
allocate(b%key(4,2),b%source_info(2,2),b%target_info(2,2))
b%source_info=0;b%target_info=0
b%source_info(2,:)=[448,448];b%target_info(2,:)=[448,448]
b%key(:,1)=14;b%key(:,2)=15
call combine_outer_edges(p,a,b)
if(.not.p%edgewise) error stop 'not edgewise'
if(any(p%sc/=[5,0,5]).or.any(p%rc/=[1,6,3])) error stop 'peer counts'
if(any(p%slot/=[1,6,10,11,12,7,8,13,14,15])) error stop 'source order'
if(any(p%address(7,:)/=[3,1,1,2,1,2,3,1,2,3])) error stop 'target edge order'
if(any(p%address(1,:)/=[11,12,13,13,14,14,14,15,15,15])) error stop 'target node order'
b%source_info(2,:)=[0,64];b%target_info(2,:)=[0,64]
call combine_outer_edges(p,a,b)
if(any(p%sc/=[2,0,3]).or.any(p%rc/=[1,3,1])) error stop 'sparse boundary counts'
if(any(p%slot/=[1,6,7,8,13])) error stop 'sparse boundary source order'
if(any(p%address(7,:)/=[3,1,1,2,1])) error stop 'sparse boundary target order'
print *, 'PASS sparse interiors and boundaries'
end program
"""
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'schedule.f90').write_text(source)
            for opt in ('-O0','-O2'):
                subprocess.run(['gfortran',opt,'-fcheck=all','-ffree-line-length-132','schedule.f90','-o','schedule'],cwd=root,check=True,capture_output=True)
                subprocess.run([str(root/'schedule')],check=True,capture_output=True)

    @unittest.skipUnless(os.environ.get('WAVETRISK_RUN_MPI_TESTS')=='1','Set WAVETRISK_RUN_MPI_TESTS=1')
    def test_actual_selector_defaults_controls_and_disagreement(self):
        text=(ROOT/'src/parallel_block_inverse.f90').read_text()
        declaration=next(line for line in text.splitlines() if 'logical :: sparse_outer_install' in line)
        source='''module actual
use mpi_f08
integer::rank
type(MPI_Comm)::comm=MPI_COMM_WORLD
'''+declaration+'\ncontains\n'+sub(text,'configure_stage192')+'''
subroutine die(message)
character(*),intent(in)::message
integer::ierr
call MPI_Abort(comm,1,ierr)
error stop message
end subroutine
subroutine check(ierr,message)
integer,intent(in)::ierr
character(*),intent(in)::message
if(ierr/=MPI_SUCCESS)call die(message)
end subroutine
end module
program check_options
use actual
integer::ierr
character(16)::arg
call MPI_Init(ierr)
call MPI_Comm_rank(comm,rank,ierr)
call get_command_argument(1,arg)
if(arg=='disagree') sparse_outer_install=rank==0
call configure_stage192
call configure_stage192
call MPI_Finalize(ierr)
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'options.f90').write_text(source)
            subprocess.run(['mpif90','-O2','-fcheck=all','options.f90','-o','options'],cwd=root,check=True,capture_output=True)
            clean={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
            for setting,arg,success in [(None,'',True),('0','',True),('1','',True),('invalid','',False),(None,'disagree',False)]:
                env=dict(clean)
                if setting is not None:env['WAVETRISK_STAGE192_SPARSE_OUTER']=setting
                result=subprocess.run(['mpirun','-n','2',str(root/'options'),arg],cwd=root,env=env,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                if success:self.assertEqual(result.stdout.count('Stage 192 requested sparse outer install = '+(setting or '0')),1)

if __name__=='__main__':unittest.main()
