"""Production patch and separate test-only ownership/lifecycle observer.

The production candidate contains no Stage223 bypass, layout override or event
files. Only the checked observer adds lifecycle evidence and empty-root layout.
"""
import subprocess
from pathlib import Path
NAMES = ('src/arch.f90', 'src/main.f90', 'src/checkpoint.f90', 'src/parallel_block_mpi.f90')


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one source anchor: '+old[:100])
    return text.replace(old, new)


def audit(tree):
    changed = {}
    p = tree/'src/arch.f90'
    s = once(p.read_text(), '\n  implicit none\n',
             '\n  use shared_mod, only : time, itime, iwrite\n\n  implicit none\n')
    s = once(s, '  public :: MPI_IN, MPI_DP, MPI_SP',
             '  public :: stage223_event\n  public :: MPI_IN, MPI_DP, MPI_SP')
    s = once(s, '    integer                            :: d, r, n_domain_floor',
             '''    integer                            :: d, r, n_domain_floor
    integer :: status223
    character(16) :: layout223''')
    s = once(s, '    call MPI_Bcast ( &\n         owner, n_glo_block, MPI_IN, 0, comm)', '''    call get_environment_variable('WAVETRISK_STAGE223_LAYOUT',layout223,status=status223)
    if(status223/=0) error stop 'Missing Stage223 ownership layout'
    select case(trim(layout223))
    case('native')
       ! Preserve the original Domain distributor exactly.
    case('empty-root')
       if(n_process<2) error stop 'Empty-root probe needs at least two ranks'
       if(rank==0) then
          do d=1,n_glo_block
             owner(d)=1+mod(d-1,n_process-1)
          end do
       end if
    case default
       error stop 'Unknown Stage223 ownership layout'
    end select
    call MPI_Bcast ( &
         owner, n_glo_block, MPI_IN, 0, comm)''')
    s = once(s, '\ncontains\n', '''
contains

  subroutine stage223_event(phase,domains,blocks)
    ! Separate files prevent concurrent MPI stdout from corrupting evidence.
    character(*),intent(in) :: phase
    integer,intent(in) :: domains,blocks
    integer,save :: sequence=0
    integer :: unit
    character(128) :: filename
    sequence=sequence+1
    write(filename,'(a,i0,a)') 'capability223-rank-',rank,'.txt'
    open(newunit=unit,file=trim(filename),status='unknown',position='append',action='write')
    write(unit,'(a,7(1x,i0),1x,es25.17)') trim(phase),sequence,rank,domains,blocks, &
         owner(1),iwrite,itime,time
    close(unit)
  end subroutine stage223_event
''')
    changed[p] = s
    p = tree/'src/parallel_block_mpi.f90'
    s = once(p.read_text(), 'module parallel_block_mpi_mod\n  use arch_mod, only : check_parallel_task_capacity',
             'module parallel_block_mpi_mod\n  use arch_mod, only : stage223_event\n  use arch_mod, only : check_parallel_task_capacity')
    s = once(s, '    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)',
             '''    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)
    call stage223_event('install',size(grid),n_local_blocks())''')
    changed[p] = s
    p = tree/'src/main.f90'
    s = once(p.read_text(), '  use arch_mod, only : abort_run, barrier, distribute_grid, glo_id, init_arch_mod, &',
             '''  use parallel_block_mod, only : n_local_blocks
  use arch_mod, only : stage223_event
  use arch_mod, only : abort_run, barrier, distribute_grid, glo_id, init_arch_mod, &''')
    anchor = '               "multistage block candidate did not retain synchronized state"\n       end if\n'
    s = once(s, anchor, anchor+"       call stage223_event('dynamics',size(grid),n_local_blocks())\n")
    s = once(s, '  end subroutine time_step',
             '''    if(.not.block_multistage_candidate.or..not.block_state_ready) &
         error stop 'Stage223 requires completed resident multistage block steps'
    call stage223_event('step',size(grid),n_local_blocks())
  end subroutine time_step''')
    changed[p] = s
    p = tree/'src/checkpoint.f90'
    s = once(p.read_text(), '  use mpi_f08', '  use mpi_f08\n  use arch_mod, only : stage223_event')
    # Accept the frozen Stage223 source and the subsequent immutable-position fix.
    position = 'global_prefix_position' if 'integer(int64) :: global_prefix_position' in s else 'domain_pos'
    anchor = f'         call distribute_checkpoint_global_state(fid,{position},global_prefix_bytes)'
    s = once(s, anchor, anchor+"\n    call stage223_event('restart',size(grid),-1)")
    changed[p] = s

    for p, text in changed.items():
        p.write_text(text)


if __name__ == '__main__':
    candidate(Path(sys.argv[1]))


def candidate(tree):
    subprocess.run(['git','apply',str(Path(__file__).parent/'stage226-production.patch')],cwd=tree,check=True)
