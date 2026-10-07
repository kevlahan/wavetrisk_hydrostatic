"""Exercise the actual old/new payload routines with four real MPI ranks.

Checks uneven signed-byte traffic, an idle rank, all-zero traffic, the exact
three global sums, collective count, and deliberate post-transfer corruption.
Builds and runs serially at O0/check and O2. No performance claim.
"""
import json,os,re,subprocess,sys,tarfile,tempfile
from pathlib import Path

PREFIX='''module proof_test
 use iso_fortran_env, only: int8,int64
 use mpi_f08, MPI_Allreduce_real => MPI_Allreduce
 implicit none
 integer :: rank,n_process,reductions
 type(MPI_Comm) :: comm
 integer(int64) :: observed(3)
 type Block_Migration_Manifest
  logical :: validated=.true.,sizes_validated=.true.,payload_validated=.false.
  integer(int64) :: total_send_nbyte=0,total_recv_nbyte=0
  integer,allocatable :: send_count(:),recv_count(:),send_displ(:),recv_displ(:),send_nbyte(:),recv_nbyte(:)
  integer,allocatable :: send_byte_count(:),recv_byte_count(:),send_byte_displ(:),recv_byte_displ(:)
  integer(int8),allocatable :: recv_payload(:)
 end type
 interface MPI_Allreduce
  module procedure reduce_scalar,reduce_vector
 end interface
contains
 subroutine reduce_scalar(s,r,n,t,op,c,ierr)
  integer(int64),intent(in) :: s
  integer(int64),intent(out) :: r
  integer,intent(in) :: n
  type(MPI_Datatype),intent(in) :: t
  type(MPI_Op),intent(in) :: op
  type(MPI_Comm),intent(in) :: c
  integer,intent(out) :: ierr
  call MPI_Allreduce_real(s,r,n,t,op,c,ierr)
  reductions=reductions+1
  observed(reductions)=r
 end subroutine
 subroutine reduce_vector(s,r,n,t,op,c,ierr)
  integer(int64),intent(in) :: s(:)
  integer(int64),intent(out) :: r(:)
  integer,intent(in) :: n
  type(MPI_Datatype),intent(in) :: t
  type(MPI_Op),intent(in) :: op
  type(MPI_Comm),intent(in) :: c
  integer,intent(out) :: ierr
  call MPI_Allreduce_real(s,r,n,t,op,c,ierr)
  reductions=reductions+1
  if(n/=3) error stop 'proof extent'
  observed=r
 end subroutine
 subroutine check_mpi(ierr,label)
  integer,intent(in) :: ierr
  character(*),intent(in) :: label
  if(ierr/=MPI_SUCCESS) call fail(label)
 end subroutine
 subroutine fail(label)
  character(*),intent(in) :: label
  integer :: ierr
  print *,label
  call MPI_Abort(comm,17,ierr)
 end subroutine
 integer function amount(src,dst,empty)
  integer,intent(in) :: src,dst
  logical,intent(in) :: empty
  amount=0
  if(empty.or.src==3.or.dst==3)return
  amount=modulo(src+2*dst+1,5)
 end function
 integer(int8) function value(src,dst,k)
  integer,intent(in) :: src,dst,k
  value=int(modulo(17*src+dst+125*k,256)-128,int8)
 end function
'''
SUFFIX='''end module
program main
 use proof_test
 implicit none
 type(Block_Migration_Manifest) :: m
 integer :: ierr,s,d,k,j,pos,expected_calls
 integer(int8),allocatable :: payload(:),expected(:)
 integer(int64) :: expected_sum,expected_bytes
 character(32) :: mode,arg
 logical :: empty
 call MPI_Init(ierr)
 comm=MPI_COMM_WORLD
 call MPI_Comm_rank(comm,rank,ierr)
 call MPI_Comm_size(comm,n_process,ierr)
 if(n_process/=4) error stop 'requires four ranks'
 call get_command_argument(1,mode)
 call get_command_argument(2,arg)
 read(arg,*) expected_calls
 empty=trim(mode)=='empty'
 allocate(m%send_count(4),m%recv_count(4),m%send_displ(4),m%recv_displ(4))
 m%send_displ=0;m%recv_displ=0
 do d=0,3
  m%send_count(d+1)=merge(1,0,amount(rank,d,empty)>0)
  m%recv_count(d+1)=merge(1,0,amount(d,rank,empty)>0)
 end do
 do d=2,4
  m%send_displ(d)=sum(m%send_count(:d-1))
  m%recv_displ(d)=sum(m%recv_count(:d-1))
 end do
 allocate(m%send_nbyte(sum(m%send_count)),m%recv_nbyte(sum(m%recv_count)))
 j=0
 do d=0,3
  k=amount(rank,d,empty)
  if(k==0)cycle
  j=j+1;m%send_nbyte(j)=k
 end do
 j=0
 do s=0,3
  k=amount(s,rank,empty)
  if(k==0)cycle
  j=j+1;m%recv_nbyte(j)=k
 end do
 m%total_send_nbyte=sum(m%send_nbyte);m%total_recv_nbyte=sum(m%recv_nbyte)
 allocate(payload(max(1,int(m%total_send_nbyte))),expected(max(1,int(m%total_recv_nbyte))))
 payload=0;expected=0;pos=0
 do d=0,3
  do k=1,amount(rank,d,empty)
   pos=pos+1;payload(pos)=value(rank,d,k)
  end do
 end do
 pos=0
 do s=0,3
  do k=1,amount(s,rank,empty)
   pos=pos+1;expected(pos)=value(s,rank,k)
  end do
 end do
 expected_sum=0;expected_bytes=0
 do s=0,3
  do d=0,3
   do k=1,amount(s,d,empty)
    expected_sum=expected_sum+int(value(s,d,k),int64)
    expected_bytes=expected_bytes+1
   end do
  end do
 end do
 reductions=0;observed=0
 call exchange_block_migration_payloads(m,payload,.false.)
 if(reductions/=expected_calls)error stop 'collective count'
 if(any(observed/=[expected_sum,expected_sum,expected_bytes]))error stop 'proof changed'
 if(.not.m%payload_validated)error stop 'unvalidated payload'
 if(any(m%recv_payload/=expected))error stop 'payload changed'
 if(rank==0)print *,'PASS actual payload routine: ',trim(mode),reductions,observed
 call MPI_Finalize(ierr)
end program
'''

def check(kit,output):
    kit=Path(kit).resolve(strict=True)
    records=[]
    with tarfile.open(kit/'sources.tar.gz') as stream,tempfile.TemporaryDirectory(prefix='rk271-proof-') as temp:
        work=Path(temp)
        for family,calls in (('stage266',3),('candidate',1)):
            text=stream.extractfile(family+'/src/parallel_block_mpi.f90').read().decode()
            routine=re.search(r'  subroutine exchange_block_migration_payloads \(.*?  end subroutine exchange_block_migration_payloads',text,re.S).group()
            for opt in ('-O0','-O2'):
                for corrupt in (False,True):
                    body=routine
                    if corrupt:
                        anchor='    local_send_checksum = 0_int64'
                        body=body.replace(anchor,"    if(rank==1.and.n_recv_byte>0) manifest%recv_payload(1)=manifest%recv_payload(1)+1_int8\n"+anchor,1)
                    source=work/'proof.f90';source.write_text(PREFIX+body+'\n'+SUFFIX)
                    binary=work/'proof'
                    subprocess.run(['mpif90',opt,'-fcheck=all','-ffree-line-length-none',str(source),'-o',str(binary)],cwd=work,check=True,capture_output=True,text=True)
                    for mode in (('corrupt',) if corrupt else ('normal','empty')):
                        env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_MCA_mpi_yield_when_idle='1')
                        command=['mpirun','--oversubscribe','--bind-to','none','-n','4',str(binary),mode,str(calls)]
                        run=subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
                        log=run.stdout+run.stderr
                        if corrupt:
                            if run.returncode==0 or 'global payload byte checksums differ' not in log:raise ValueError('Corruption was not caught: '+log)
                        elif run.returncode or 'PASS actual payload routine' not in log:raise ValueError(log)
                        records.append(dict(family=family,optimization=opt,mode=mode,expected_collectives=calls,passed=True,returncode=run.returncode))
                        print('PASS',family,opt,mode,flush=True)
    report=dict(status='passed',cases=records,source_manifest_sha256=__import__('hashlib').sha256((kit/'source-files.json').read_bytes()).hexdigest(),
                scope='Twelve actual-routine four-rank MPI cases; exact three integer sums/payloads and retained corruption rejection, with three old or one new proof collective.')
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':check(*sys.argv[1:])
