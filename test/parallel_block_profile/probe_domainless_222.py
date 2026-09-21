"""Test-only four-rank geometry layout: rank 3 owns no Domains.

Adaptive block ownership is left to the real block distributor. This does not
relax the startup task-count guard or claim more-tasks-than-Domains support.
Apply only to a disposable checked source tree.
"""
from pathlib import Path
import sys
from audit_cleanup_222 import once, routine


def instrument(tree):
    p = tree/'src/arch.f90'
    p.write_text(once(p.read_text(), '    call MPI_Bcast ( &\n         owner, n_glo_block, MPI_IN, 0, comm)', '''    if(n_process/=4) error stop 'Stage222 domainless probe requires four tasks'
    if(rank==0) then
       do d=1,n_glo_block
          owner(d)=mod(d-1,n_process-1)
       end do
       write(*,'(a)') 'Stage222 domainless probe: geometry owners 0,1,2; block ranks 0,1,2,3'
    end if
    call MPI_Bcast ( &
         owner, n_glo_block, MPI_IN, 0, comm)'''))
    p = tree/'src/parallel_block_mpi.f90'
    p.write_text(once(p.read_text(),
                     '    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)',
                     '''    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)
    if(size(grid)==0) then
       write(*,'(a,i0,a,i0)') 'Stage222 domainless owner rank = ',rank,' blocks = ',n_local_blocks()
       if(n_local_blocks()==0) call fail('domainless probe has no owned block work')
    end if'''))
    s = p.read_text()
    if 'logical :: native_publication_ready' in s:
        if 'ieee_quiet_nan' not in s[:s.index('  implicit none')]:
            s = once(s, '  use ieee_arithmetic, only : ieee_is_finite',
                     '  use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan')
        s = routine(s, 'prepare_parallel_block_grid_change', lambda t: once(t,
                    '    if (.not.native_publication_ready) call synchronize_parallel_block_checkpoint(.true.)',
                    '''    if(native_publication_ready) then
       ! Do not recreate the omitted writebacks in this probe. Poison their
       ! scratch so a later reader cannot benefit from old staging/transport.
       block_writeback_plan%scalar_domain_stage=ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%vector_domain_stage=ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%scalar_send_buffer=ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%vector_send_buffer=ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%scalar_recv_buffer=ieee_value(0.0_dp,ieee_quiet_nan)
       block_writeback_plan%vector_recv_buffer=ieee_value(0.0_dp,ieee_quiet_nan)
       write(*,'(a,i0)') 'Stage222 discarded grid buffers rank = ',rank
    end if
    if (.not.native_publication_ready) call synchronize_parallel_block_checkpoint(.true.)'''))
        p.write_text(s)


if __name__ == '__main__':
    instrument(Path(sys.argv[1]))
