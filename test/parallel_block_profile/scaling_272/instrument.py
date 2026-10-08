"""Diagnostic copy only: bounded timers and actual ownership/work accounting."""
from pathlib import Path
import instrument_base as base
from instrument_base import inject, scoped_inject, wrap, mirror_scopes


def instrument(root):
    root=Path(root)
    changed=base.instrument(root,True)
    path=root/'src/multi_level.f90';text=path.read_text()
    text=wrap(text,'basic_operators','BASIC')
    text=scoped_inject(text,'basic_operators','    do d = 1, size(grid)\n',
        '    do d = 1, size(grid)\n       call cost_work(C_BASIC,int(grid(d)%lev(l)%length,int64),1_int64,0_int64)\n')
    text=mirror_scopes(text,'block_tendency_compatibility_ml','detail',{'DP_PRIMITIVE':'PRIMITIVE'})
    text=mirror_scopes(text,'finish_native_mass','detail',{'DP_MASS_WAIT':'MASS_WAIT'})
    text=scoped_inject(text,'basic_operators','    do d = 1, size(grid)\n',
        '    call cost_enter(C_PRODUCER_PATCH)\n    do d = 1, size(grid)\n')
    text=scoped_inject(text,'basic_operators','       nullify (bernoulli, exner)\n    end do\n',
        '       nullify (bernoulli, exner)\n    end do\n    call cost_leave(C_PRODUCER_PATCH)\n')
    path.write_text(text)
    path=root/'src/comm_mpi.f90';text=path.read_text()
    assert text.startswith('module comm_mpi_mod\n')
    text=text.replace('module comm_mpi_mod\n','module comm_mpi_mod\n  use cost_observer_mod\n',1)
    for routine in ('update_bdry__finish1_0','update_bdry__finish1_1','update_bdry__finish1_2'):
        call='       call MPI_Waitall (nreq, req(1:nreq), MPI_STATUSES_IGNORE)\n'
        text=scoped_inject(text,routine,call,
            '       call cost_enter(C_DOMAIN_WAIT)\n'+call+'       call cost_leave(C_DOMAIN_WAIT)\n')
    path.write_text(text);changed.append('src/comm_mpi.f90')
    make=root/'Makefile';make.write_text(make.read_text()+'\n$(BUILD_DIR)/comm_mpi.o: $(BUILD_DIR)/cost_observer.o\n')
    path=root/'src/parallel_block_mpi.f90';text=path.read_text()
    text=scoped_inject(text,'begin_block_scalar_divergence_capture','    call cost_enter(C_CAPTURE_BEGIN)\n',
        '    call cost_enter(C_CAPTURE_BEGIN)\n    call cost_census272\n')
    text=scoped_inject(text,'evaluate_candidate_block_scalar_restriction','    count_local = [statistics%parent_patch_count, &\n',
        '    call cost_work(C_SCALAR_LOCAL,statistics%parent_patch_count, &\n         statistics%restricted_edge_count,statistics%divergence_replay_count)\n    count_local = [statistics%parent_patch_count, &\n')
    # Each actual transport branch has one matching count-array pair. Logical
    # remote doubles exclude the self copy; names run opposite to MPI direction.
    anchor='    call check_mpi(ierr,"MPI_Alltoallv scalar-divergence input")\n'
    counters='''    if (deduplicated_transport) then
       call cost_capture_volume272(block_scalar_divergence_plan%deduplicated_recv_count, &
            block_scalar_divergence_plan%deduplicated_send_count)
    else if (block_scalar_divergence_plan%full_transport) then
       call cost_capture_volume272(block_scalar_divergence_plan%recv_count, &
            block_scalar_divergence_plan%send_count)
    else
       call cost_capture_volume272(block_scalar_divergence_plan%production_recv_count, &
            block_scalar_divergence_plan%production_send_count)
    end if
'''
    text=scoped_inject(text,'finalize_block_scalar_divergence_capture',anchor,anchor+counters)
    extra='''
  subroutine cost_capture_volume272(sent,received)
    integer, intent(in) :: sent(:),received(:)
    call cost_work(C_CAPTURE_MPI,1_int64, &
         sum(int(sent,int64))-int(sent(rank+1),int64), &
         sum(int(received,int64))-int(received(rank+1),int64))
  end subroutine

  subroutine cost_census272
    integer :: d,l,i
    integer(int64) :: domain_patches,block_patches
    call cost_enter(C_OWNERS)
    domain_patches=0_int64
    block_patches=0_int64
    do d=1,size(grid)
       do l=level_start,level_end
          domain_patches=domain_patches+int(grid(d)%lev(l)%length,int64)
       end do
    end do
    do i=1,n_local_blocks()
       block_patches=block_patches+int(local_block_patch_count(local_block_catalog(i)),int64)
    end do
    call cost_work(C_OWNERS,int(size(grid),int64),int(n_local_blocks(),int64),1_int64)
    call cost_work(C_PATCHES,domain_patches,block_patches,1_int64)
    call cost_leave(C_OWNERS)
  end subroutine

'''
    text=inject(text,'end module parallel_block_mpi_mod',extra+'end module parallel_block_mpi_mod')
    path.write_text(text)
    return sorted(changed)
