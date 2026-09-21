"""Check existing Domain interiors before final grid-change republication.

Run with production oracles OFF: the existing comparison then requires exact
equality. The actual writeback and subsequent saved-alias publication remain.
This probes liveness/version equivalence; it is not a performance experiment.
"""
from pathlib import Path
import sys
from audit_cleanup_222 import once, routine


def instrument(tree):
    p = tree/'src/parallel_block_mpi.f90'
    p.write_text(routine(p.read_text(), 'prepare_parallel_block_grid_change', lambda t: once(t,
                        '    call synchronize_parallel_block_checkpoint(.true.)', '''    if(production_inverse_consumer_aliases_pending) then
       if(production_inverse_completed_mode/=0) call fail('grid overlap probe before final refresh')
       call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
       call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)
       write(*,'(a,i0)') 'Stage222 pre-grid publication matches rank = ',rank
    end if
    call synchronize_parallel_block_checkpoint(.true.)''')))


if __name__ == '__main__':
    instrument(Path(sys.argv[1]))
