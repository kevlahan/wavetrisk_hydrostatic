"""Apply the Stage254 timing probe to a verified, unmodified Stage251 source tree.

The probe changes only owner-program timers and an end-of-window report hook.
It does not alter the numerical execution order or the Stage251 feature set.
"""
import hashlib
from pathlib import Path
import sys

PROGRAM_SHA = 'eff68d917b06eebea30aa81f12ec7bb1c8163de3c502e27b087025da4a66a8e8'
MPI_SHA = '5a629218606a164e38cba80bc5101dcd17c6f1f492db8374e0a210b9869db505'
INSTRUMENTED_PROGRAM_SHA = '6d12e52f44fc0ded0a65326f4a48c08c743029c6c8384916d9580257198d22d7'
INSTRUMENTED_MPI_SHA = '119a9a055dcf7f86185ee82adaad9b099daf6b2d94ca53ce5bb835ea45e77f12'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace(text, before, after):
    if text.count(before) != 1:
        raise ValueError('Instrumentation anchor absent or ambiguous: '+before[:75])
    return text.replace(before, after, 1)


def instrument(source):
    source = Path(source)
    program = source/'src/parallel_block_program.f90'
    mpi = source/'src/parallel_block_mpi.f90'
    if sha(program) != PROGRAM_SHA or sha(mpi) != MPI_SHA:
        raise ValueError('Stage251 source differs before instrumentation')
    p = program.read_text()
    p = replace(p, '     logical :: compiled=.false.\n  end type\ncontains\n',
                '     logical :: compiled=.false.,profile_velocity=.false.\n  end type\n'
                '  integer,parameter :: PC_EXEC=1,PC_FETCH=2,PC_READ=3,PC_READ_MPI=4, &\n'
                '       PC_EVAL=5,PC_DELIVERY=6,PC_DELIVERY_MPI=7,PROFILE_DEPTH_MAX=16\n'
                '  real(dp),save :: owner_seconds(7)=0.0_dp,depth_seconds(0:PROFILE_DEPTH_MAX)=0.0_dp\n'
                '  integer(int64),save :: owner_calls(3)=0_int64,depth_calls(0:PROFILE_DEPTH_MAX)=0_int64\n'
                '  logical,save :: owner_profile_initialized=.false.,owner_profile_enabled=.false.\n'
                'contains\n'
                '  logical function owner_profile_active() result(enabled)\n'
                '    character(8)::value\n'
                '    integer::status\n'
                '    if(.not.owner_profile_initialized) then\n'
                '       value=""\n'
                '       call get_environment_variable("WAVETRISK_PROFILE_OWNER_PROGRAM",value,status=status)\n'
                '       owner_profile_enabled=status==0.and.trim(value)=="1"\n'
                '       owner_profile_initialized=.true.\n'
                '    end if\n'
                '    enabled=owner_profile_enabled\n'
                '  end function\n')
    p = replace(p, '  subroutine program_read(p,q)\n',
                '  subroutine program_read(p,q,profile_depth)\n')
    p = replace(p, '    type(Owner_Read),intent(inout)::q\n    integer::i,s,j,ierr\n    if(p%layer<=0)',
                '    type(Owner_Read),intent(inout)::q\n'
                '    integer,intent(in),optional::profile_depth\n'
                '    integer::i,s,j,ierr\n'
                '    real(dp)::start_time,mpi_start,elapsed\n'
                '    if(p%profile_velocity) start_time=MPI_Wtime()\n'
                '    if(p%layer<=0)')
    p = replace(p, '    call MPI_Alltoallv(q%send,q%rc,q%rd,MPI_DOUBLE_PRECISION, &\n'
                '         q%received,q%sc,q%sd,MPI_DOUBLE_PRECISION,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die(\'result delivery\')\n',
                '    if(p%profile_velocity) mpi_start=MPI_Wtime()\n'
                '    call MPI_Alltoallv(q%send,q%rc,q%rd,MPI_DOUBLE_PRECISION, &\n'
                '         q%received,q%sc,q%sd,MPI_DOUBLE_PRECISION,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die(\'result delivery\')\n'
                '    if(p%profile_velocity) then\n'
                '       elapsed=MPI_Wtime()-mpi_start\n'
                '       if(present(profile_depth)) then\n'
                '          owner_seconds(PC_READ_MPI)=owner_seconds(PC_READ_MPI)+elapsed\n'
                '       else\n'
                '          owner_seconds(PC_DELIVERY_MPI)=owner_seconds(PC_DELIVERY_MPI)+elapsed\n'
                '       end if\n'
                '    end if\n')
    p = replace(p, '    do i=1,size(q%value)\n'
                '       q%value(i)=q%orientation(i)*q%received(q%address(i))\n'
                '    end do\n  end subroutine\n  subroutine program_execute',
                '    do i=1,size(q%value)\n'
                '       q%value(i)=q%orientation(i)*q%received(q%address(i))\n'
                '    end do\n'
                '    if(p%profile_velocity) then\n'
                '       elapsed=MPI_Wtime()-start_time\n'
                '       if(present(profile_depth)) then\n'
                '          if(profile_depth<0.or.profile_depth>PROFILE_DEPTH_MAX) call die("profile depth")\n'
                '          owner_seconds(PC_READ)=owner_seconds(PC_READ)+elapsed\n'
                '          depth_seconds(profile_depth)=depth_seconds(profile_depth)+elapsed\n'
                '          depth_calls(profile_depth)=depth_calls(profile_depth)+1_int64\n'
                '          owner_calls(2)=owner_calls(2)+1_int64\n'
                '       else\n'
                '          owner_seconds(PC_DELIVERY)=owner_seconds(PC_DELIVERY)+elapsed\n'
                '          owner_calls(3)=owner_calls(3)+1_int64\n'
                '       end if\n'
                '    end if\n'
                '  end subroutine\n  subroutine program_execute')
    p = replace(p, '    real(dp),allocatable::coeff(:)\n    integer::depth,i,a,n\n'
                '    if(.not.p%compiled) call die(\'execute before compile\')\n',
                '    real(dp),allocatable::coeff(:)\n'
                '    real(dp)::execute_start,section_start\n'
                '    integer::depth,i,a,n\n'
                '    if(.not.p%compiled) call die(\'execute before compile\')\n'
                '    p%profile_velocity=present(external).and.owner_profile_active()\n'
                '    if(p%profile_velocity) execute_start=MPI_Wtime()\n')
    p = replace(p, '    call fetch_resident_operands(p%resident,k)\n'
                '    do depth=0,p%maxdepth\n'
                '       call program_read(p,p%read(depth))\n'
                '       do i=1,size(p%owned)\n',
                '    if(p%profile_velocity) section_start=MPI_Wtime()\n'
                '    call fetch_resident_operands(p%resident,k)\n'
                '    if(p%profile_velocity) owner_seconds(PC_FETCH)= &\n'
                '         owner_seconds(PC_FETCH)+MPI_Wtime()-section_start\n'
                '    do depth=0,p%maxdepth\n'
                '       call program_read(p,p%read(depth),depth)\n'
                '       if(p%profile_velocity) section_start=MPI_Wtime()\n'
                '       do i=1,size(p%owned)\n')
    p = replace(p, '          end associate\n'
                '       end do\n'
                '    end do\n'
                '  end subroutine\n'
                '  function program_value',
                '          end associate\n'
                '       end do\n'
                '       if(p%profile_velocity) owner_seconds(PC_EVAL)= &\n'
                '            owner_seconds(PC_EVAL)+MPI_Wtime()-section_start\n'
                '    end do\n'
                '    if(p%profile_velocity) then\n'
                '       owner_seconds(PC_EXEC)=owner_seconds(PC_EXEC)+MPI_Wtime()-execute_start\n'
                '       owner_calls(1)=owner_calls(1)+1_int64\n'
                '    end if\n'
                '  end subroutine\n'
                '  function program_value')
    p = replace(p, '  subroutine die(message)\n',
                '  subroutine report_owner_program_profile(reset)\n'
                '    logical,intent(in)::reset\n'
                '    real(dp)::seconds_sum(7),seconds_min(7),seconds_max(7),depth_sum(0:PROFILE_DEPTH_MAX)\n'
                '    integer(int64)::calls_sum(3),depth_count_sum(0:PROFILE_DEPTH_MAX)\n'
                '    integer::i,d,ierr\n'
                '    character(len=18),parameter::names(7)=[character(len=18):: &\n'
                '         "execute","resident fetch","internal read","read MPI", &\n'
                '         "producer loop","delivery read","delivery MPI"]\n'
                '    if(.not.owner_profile_active()) return\n'
                '    call MPI_Allreduce(owner_seconds,seconds_sum,7,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile seconds sum")\n'
                '    call MPI_Allreduce(owner_seconds,seconds_min,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile seconds min")\n'
                '    call MPI_Allreduce(owner_seconds,seconds_max,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile seconds max")\n'
                '    call MPI_Allreduce(owner_calls,calls_sum,3,MPI_INTEGER8,MPI_SUM,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile calls")\n'
                '    call MPI_Allreduce(depth_seconds,depth_sum,PROFILE_DEPTH_MAX+1, &\n'
                '         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile depth seconds")\n'
                '    call MPI_Allreduce(depth_calls,depth_count_sum,PROFILE_DEPTH_MAX+1, &\n'
                '         MPI_INTEGER8,MPI_SUM,comm,ierr)\n'
                '    if(ierr/=MPI_SUCCESS) call die("profile depth calls")\n'
                '    if(rank==0) then\n'
                '       write(6,\'(a,3(i0,1x))\') "Owner program calls execute/read/delivery = ", &\n'
                '            calls_sum/int(n_process,int64)\n'
                '       do i=1,7\n'
                '          write(6,\'(a,a,a,3(es13.5,1x))\') "Owner program ",trim(names(i)), &\n'
                '               " seconds min/avg/max = ",seconds_min(i),seconds_sum(i)/real(n_process,dp),seconds_max(i)\n'
                '       end do\n'
                '       do d=0,PROFILE_DEPTH_MAX\n'
                '          if(depth_count_sum(d)==0_int64) cycle\n'
                '          write(6,\'(a,i0,a,i0,a,es13.5)\') "Owner program depth ",d, &\n'
                '               " calls = ",depth_count_sum(d)/int(n_process,int64), &\n'
                '               " read seconds avg = ",depth_sum(d)/real(n_process,dp)\n'
                '       end do\n'
                '    end if\n'
                '    if(reset) then\n'
                '       owner_seconds=0.0_dp;depth_seconds=0.0_dp\n'
                '       owner_calls=0_int64;depth_calls=0_int64\n'
                '    end if\n'
                '  end subroutine\n'
                '  subroutine die(message)\n')
    p = replace(p, '  public :: program_value\n',
                '  public :: program_value,report_owner_program_profile\n')
    program.write_text(p)

    m = mpi.read_text()
    m = replace(m, '  use parallel_block_diffusion_mod, only: owned_mass_node\n',
                '  use parallel_block_diffusion_mod, only: owned_mass_node\n'
                '  use parallel_block_program_mod, only: report_owner_program_profile\n')
    m = replace(m, '    call report_block_detail(clear_after)\n',
                '    call report_owner_program_profile(clear_after)\n'
                '    call report_block_detail(clear_after)\n')
    mpi.write_text(m)
    result = {'program_sha256': sha(program), 'mpi_sha256': sha(mpi)}
    if result != {'program_sha256': INSTRUMENTED_PROGRAM_SHA,
                  'mpi_sha256': INSTRUMENTED_MPI_SHA}:
        raise ValueError('Timing-probe output hash differs')
    return result


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: instrument.py SOURCE_TREE')
    print(instrument(sys.argv[1]))
