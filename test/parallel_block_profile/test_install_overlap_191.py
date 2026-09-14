"""Exercise actual ownership-transfer and deferred MPI routines with small stores."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]


def sub(text,name):
    return re.search(r'^\s*subroutine '+name+r'\b.*?^\s*end subroutine(?: '+name+r')?\s*$',text,re.M|re.S|re.I).group()


class InstallOverlap191(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('WAVETRISK_RUN_MPI_TESTS')=='1','Set WAVETRISK_RUN_MPI_TESTS=1 for selector MPI test')
    def test_default_and_explicit_selectors(self):
        text=(ROOT/'src/parallel_block_mpi.f90').read_text()
        declarations='\n'.join(line for line in text.splitlines() if 'logical, save :: stage191_options' in line)
        function=re.search(r'  logical function block_stage191_option\(option\).*?  end function block_stage191_option',text,re.S).group()
        source='''module actual
use mpi_f08
integer :: rank
type(MPI_Comm) :: comm=MPI_COMM_WORLD
'''+declarations+'\ncontains\n'+function+'''
subroutine fail(message)
character(*),intent(in)::message
error stop message
end subroutine
subroutine check_mpi(ierr,message)
integer,intent(in)::ierr
character(*),intent(in)::message
if(ierr/=MPI_SUCCESS)error stop message
end subroutine
end module
program check
use actual
integer::ierr,i,j
call MPI_Init(ierr)
call MPI_Comm_rank(comm,rank,ierr)
i=merge(1,0,block_stage191_option(1))
j=merge(1,0,block_stage191_option(2))
write(*,'(a,2i2)') 'SELECTORS',i,j
call MPI_Finalize(ierr)
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'options.f90').write_text(source)
            subprocess.run(['mpif90','-O2','-fcheck=all','options.f90','-o','options'],cwd=p,check=True,capture_output=True)
            clean={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
            for overrides,expected in [({},'SELECTORS 1 0'),
                ({'WAVETRISK_STAGE191_MASS_OVERLAP':'1'},'SELECTORS 1 1'),
                ({'WAVETRISK_STAGE191_INSTALL':'0'},'SELECTORS 0 0')]:
                result=subprocess.run(['mpirun','-n','2',str(p/'options')],cwd=p,env={**clean,**overrides},check=True,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.stdout.count(expected),2)

    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_move_every_component_and_replacement(self):
        text=(ROOT/'src/parallel_block.f90').read_text()
        body=re.search(r'  type, public :: Block_Data\n.*?  end type Block_Data',text,re.S).group()
        meta=[];arrays=[];types=set()
        for line in body.splitlines()[1:-1]:
            if '::' not in line:continue
            decl,var=line.split('::');name=var.split('=')[0].split('(')[0].strip()
            if 'allocatable' not in decl:meta.append(name);continue
            derived=re.search(r'type\((\w+)\)',decl)
            if derived:types.add(derived[1])
            arrays.append((name,2 if '(:,:)' in var else 1,derived is not None))
        # Minimal leaf payload types; Block_Data and its move routine are actual source.
        stubs='\n'.join(f'type :: {t}\n integer :: tag\nend type' for t in sorted(types))
        init=[];check=[]
        for i,name in enumerate(meta,1):
            init.append(f'a%{name}={i}')
            check.append(f'if(b%{name}/={i}) error stop "metadata {name}"')
        for i,(name,rank,derived) in enumerate(arrays,1):
            shape='-1:2' if rank==1 else '-1:2,0:2'
            dest='a%'+name+('%tag' if derived else '')
            value=str(i) if derived or name=='neigh_class' else str(i)+'.125_dp'
            init.extend([f'allocate(a%{name}({shape}))',f'{dest}={value}',f'allocate(b%{name}(1'+(',1' if rank==2 else '')+'))'])
            dest='b%'+name+('%tag' if derived else '')
            check.extend([f'if(allocated(a%{name})) error stop "source still owns {name}"',
                          f'if(any({dest}/={value})) error stop "payload {name}"',
                          f'if(lbound(b%{name},1)/=-1) error stop "lower bound {name}"'])
        source='module actual\nuse iso_fortran_env, only: dp=>real64\n'+stubs+'\n'+body+'\ncontains\n'+sub(text,'move_block_storage')+'\nend module\n'
        source+='program check\nuse actual\ntype(Block_Data)::a,b\n'+'\n'.join(init)+'\ncall move_block_storage(a,b)\n'+'\n'.join(check)+'\nprint *, "PASS moved all components"\nend program\n'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'move.f90').write_text(source)
            for opt in ['-O0','-O2']:
                subprocess.run(['gfortran',opt,'-fcheck=all','-ffp-contract=off','move.f90','-o','move'],cwd=p,check=True,capture_output=True)
                result=subprocess.run([str(p/'move')],check=True,text=True,capture_output=True)
                self.assertIn('PASS',result.stdout)

    @unittest.skipUnless(os.environ.get('WAVETRISK_RUN_MPI_TESTS')=='1','Set WAVETRISK_RUN_MPI_TESTS=1 for two-rank MPI test')
    def test_deferred_mass_buffers_and_guards(self):
        text=(ROOT/'src/multi_level.f90').read_text()
        plan=re.search(r'  type :: Mass_Boundary_Plan\n.*?  end type',text,re.S).group()
        header='''module actual
use mpi_f08
use iso_fortran_env, only: dp=>real64,int64
integer,parameter :: AT_NODE=1,AT_EDGE=2,DP_MASS_WAIT=1
integer :: n_process
integer(int64) :: mass_work(6)=0
type(MPI_Comm) :: comm=MPI_COMM_WORLD
type(MPI_Datatype) :: MPI_DP=MPI_DOUBLE_PRECISION
type :: Work
 real(dp),allocatable :: flux(:),tendency(:,:)
end type
type(Work) :: native_mass(1)
'''+plan+'''
type(Mass_Boundary_Plan),allocatable,asynchronous :: mass_boundary(:,:)
contains
subroutine detail_enter(p)
integer,intent(in)::p
end subroutine
subroutine detail_leave(p)
integer,intent(in)::p
end subroutine
'''
        routines='\n'.join(sub(text,n) for n in ['begin_native_mass','finish_native_mass','exchange_native_mass'])
        driver='''
end module
program check
use actual
integer::i,k,r,peer,ierr,pos
real(dp)::started
character(16)::mode
call MPI_Init(ierr)
call MPI_Comm_rank(comm,r,ierr)
call MPI_Comm_size(comm,n_process,ierr)
if(n_process/=2) error stop "two ranks required"
call get_command_argument(1,mode)
peer=1-r
allocate(mass_boundary(2,1),native_mass(1)%flux(3),native_mass(1)%tendency(3,2))
do pos=1,2
 associate(p=>mass_boundary(pos,1))
 allocate(p%send_count(2),p%recv_count(2),p%send_displ(2),p%recv_displ(2))
 p%send_count=0;p%recv_count=0;p%send_displ=0;p%recv_displ=0
 p%send_count(peer+1)=1;p%recv_count(peer+1)=1
 allocate(p%source(2,1),p%destination(3,1),p%local(5,1),p%send_value(1),p%recv_value(1))
 p%source(:,1)=[1,1];p%destination(:,1)=[1,3,1];p%local(:,1)=[1,1,1,2,1]
 end associate
end do
do i=1,20
 k=1+modulo(i,2)
 native_mass(1)%tendency=-999.0_dp
 native_mass(1)%tendency(1,k)=100*r+i
 if(r==1) then
  started=MPI_Wtime()
  do while(MPI_Wtime()-started<0.002_dp)
  end do
 end if
 call begin_native_mass(AT_NODE,1,k)
 if(mode=='reuse')call begin_native_mass(AT_NODE,1,k)
 if(mode=='layer')call finish_native_mass(AT_NODE,1,3-k)
 ! Independent velocity-source-like work while requests own the buffers.
 native_mass(1)%flux=real(i,dp)
 call finish_native_mass(AT_NODE,1,k)
 if(mass_boundary(AT_NODE,1)%pending)error stop "pending after finish"
 if(native_mass(1)%tendency(2,k)/=100*r+i)error stop "local copy"
 if(native_mass(1)%tendency(3,k)/=100*peer+i)error stop "remote copy"
end do
call MPI_Finalize(ierr)
print *, 'PASS deferred mass',r
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'mass.f90').write_text(header+routines+driver)
            for opt in ['-O0','-O2']:
                subprocess.run(['mpif90',opt,'-fcheck=all','-ffp-contract=off','mass.f90','-o','mass'],cwd=p,check=True,capture_output=True)
                for mode in ['valid','reuse','layer']:
                    run=subprocess.run(['mpirun','-n','2',str(p/'mass'),mode],cwd=p,capture_output=True,text=True,timeout=30)
                    if mode=='valid':self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                    else:self.assertNotEqual(run.returncode,0,mode)

if __name__=='__main__': unittest.main()
