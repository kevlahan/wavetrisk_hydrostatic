#!/usr/bin/env python3
"""Add bounded integer work accounting to a disposable Stage325 source copy.

Every replacement is unique and fail-closed. No source tree is assumed writable.
This must run after the owner-map patch and before the independent semantic sink.
"""
from pathlib import Path
import argparse, hashlib, json
from boundary_hook325 import patch as boundary_patch,SCHEMA as BOUNDARY_SCHEMA

HERE=Path(__file__).resolve().parent
SCHEMA={}
OPERATIONS=[]

def apply(root):
    root=Path(root); src=root/'src'; changed={}
    def read(name):
        if name not in changed: changed[name]=(src/name).read_text()
        return changed[name]
    def replace(name,old,new):
        text=read(name)
        if text.count(old)!=1: raise ValueError(f'account325 anchor {name}: {old[:100]!r} count={text.count(old)}')
        changed[name]=text.replace(old,new)
        OPERATIONS.append({'file':name,'anchor_sha256':hashlib.sha256(old.encode()).hexdigest()})
    def before(name,anchor,code): replace(name,anchor,code+'\n'+anchor)
    def after(name,anchor,code): replace(name,anchor,anchor+'\n'+code)
    def use(name):
        line,second=read(name).splitlines()[:2]
        replace(name,line+'\n'+second,line+'\n  use audit316_mod, only: work325_enabled,work325_record\n'+second)
    def schema(tag,fields,peer_unit=None):
        SCHEMA[tag]={'fields':fields.split(), 'peer_unit':peer_unit}
    def record(tag,phase,expr,fields,send=None,recv=None):
        schema(tag,fields, 'bytes' if send is not None else None)
        suffix=''
        if send is not None: suffix=', &\n         '+send+', &\n         '+recv
        return f"    if(work325_enabled()) call work325_record('{tag}',{phase}, &\n         {expr}{suffix})"

    name='audit316.f90'
    before(name,'  interface', (HERE/'work325_decl.f90').read_text().rstrip())
    after(name,'contains',(HERE/'work325_impl.f90').read_text().rstrip())
    after(name,'    initialized316=.true.','    call work325_init(rank)')
    before(name,'    finished316=.true.','    call work325_finish()')

    name='parallel_block_mpi.f90'; use(name)
    after(name,'    mode=owner_bridge322_mode()\n    if(.not.mode_agreed) then',
          record('selector.consensus','0','[1_int64,2_int64*4_int64]','collectives payload_bytes'))
    if '    if(mode==4) call check_locality325_permutation()' in read(name):
        before(name,'    if(mode==4) call check_locality325_permutation()',
               "    if(mode==4) then\n"+record('selector.locality.guard','4','[640_int64,256_int64]',
               'catalogue_iterations permutation_iterations')+'\n    end if')
    before(name,'    if(mode==0) return\n    call check_owner_bridge322_catalog(.false.)',
           record('selector.call','mode','[int(size(selected),int64),int(merge(0,640,mode==0),int64)]',
                  'native_catalogue_entries overridden_entries'))
    before(name,'    if(audit316_enabled()) call audit_migration316(manifest)',
           record('migration','0','[int(size(block_source),int64), &\n         int(size(block_source)-manifest%n_send,int64),int(manifest%n_send,int64),int(manifest%n_recv,int64), &\n         manifest%total_send_nbyte,manifest%total_recv_nbyte, &\n         int(merge(n_local_blocks(),0,transfer_storage),int64),17_int64, &\n         int(size(send_payload),int64)+int(size(manifest%recv_payload),int64)]',
                  'source_blocks retained_blocks remote_send_blocks remote_receive_blocks packed_copy_checksum_send_bytes unpack_checksum_receive_bytes storage_moves direct_collectives live_aggregate_payload_bytes',
                  'int(manifest%send_byte_count,int64)','int(manifest%recv_byte_count,int64)'))
    after(name,'    send_nbyte = 0\n    seen_source = 0',
          "    ! Diagnostic extent census includes retained blocks; these are not bytes copied.\n"
          "    if(work325_enabled()) then\n"
          "       do ib=1,size(block_source)\n"
          "          b=block_source_catalog_index(ib)\n"
          "          call work325_record('migration.source',b,[int(packed_block_nbyte(block_source(ib)),int64), &\n"
          "               int(merge(1,0,block_catalog(b)%owner/=rank),int64)])\n"
          "       end do\n"
          "    end if")
    schema('migration.source','serialized_equivalent_bytes remote_flag')
    # The live component-family and payload-family branches are retained exactly.
    before(name,'    call detail_enter(DP_WRITEBACK_MPI)\n    if (component /= BLOCK_WRITEBACK_VECTOR) call MPI_Alltoallv(',
           "    if(component/=BLOCK_WRITEBACK_VECTOR) then\n"+
           record('writeback.scalar','family','[int(sum(block_writeback_plan%scalar_send_count),int64), &\n         int(sum(block_writeback_plan%scalar_recv_count),int64)]','gather_values receive_values',
                  '8_int64*int(block_writeback_plan%scalar_send_count,int64)',
                  '8_int64*int(block_writeback_plan%scalar_recv_count,int64)')+'\n    end if\n'+
           "    if(component/=BLOCK_WRITEBACK_SCALAR) then\n"+
           record('writeback.vector','family','[int(sum(block_writeback_plan%vector_send_count),int64), &\n         int(sum(block_writeback_plan%vector_recv_count),int64)]','gather_values receive_values',
                  '8_int64*int(block_writeback_plan%vector_send_count,int64)',
                  '8_int64*int(block_writeback_plan%vector_recv_count,int64)')+'\n    end if')
    before(name,'    call detail_enter(DP_WRITEBACK_MPI)\n    call MPI_Alltoallv(block_writeback_plan%scalar_send_buffer,sc,sd,MPI_DOUBLE_PRECISION,',
           record('writeback.temperature','1','[int(sum(sc),int64),int(sum(rc),int64),int(sum(rc),int64)]',
                  'gather_values received_values expansion_copy_values',
                  '8_int64*int(sc,int64)','8_int64*int(rc,int64)'))
    before(name,'    if (.not. transfer_vector_only) then\n       call MPI_Alltoallv(',
           "    if(.not.transfer_vector_only) then\n"+
           record('reverse.scalar','payload_family','[int(sum(block_writeback_plan%scalar_recv_count),int64), &\n         int(sum(block_writeback_plan%scalar_send_count),int64)]','gather_values receive_values',
                  '8_int64*int(block_writeback_plan%scalar_recv_count,int64)',
                  '8_int64*int(block_writeback_plan%scalar_send_count,int64)')+'\n    end if\n'+
           record('reverse.vector','payload_family','[int(sum(block_writeback_plan%vector_recv_count),int64), &\n         int(sum(block_writeback_plan%vector_send_count),int64)]','gather_values receive_values',
                  '8_int64*int(block_writeback_plan%vector_recv_count,int64)',
                  '8_int64*int(block_writeback_plan%vector_send_count,int64)'))
    after(name,'       imported_patch_count = imported_patch_count + local_patch\n    end do',
          record('reverse.local','payload_family','[int(imported_patch_count,int64), &\n         int(imported_patch_count,int64)*int(merge(0,block_writeback_plan%scalar_patch_nvalue, &\n         transfer_vector_only),int64), &\n         int(imported_patch_count,int64)*int(block_writeback_plan%vector_patch_nvalue,int64)]',
                 'patches scalar_direct_install_values vector_direct_install_values'))
    anchor='       reconstructed_patch_count = reconstructed_patch_count + n_patch\n    end do\n\n    do r = 1,n_process\n       pos_scalar = block_writeback_plan%scalar_recv_displ(r) + 1'
    replace(name,anchor,anchor.replace('\n    do r = 1,n_process', '\n'+record('writeback.local','payload_family',
            '[int(reconstructed_patch_count,int64)]','direct_staged_patches')+'\n    do r = 1,n_process'))

    changed['parallel_block_mpi.f90']=boundary_patch(changed['parallel_block_mpi.f90'])
    SCHEMA.update(BOUNDARY_SCHEMA)

    name='parallel_block.f90'; use(name)
    after(name,'  type(Block_Data), intent(out) :: destination',
          record('storage.move','0','[36_int64]','move_alloc_calls'))
    # Literal is certified against exactly the move_alloc calls in this procedure.
    part=read(name).split('subroutine move_block_storage(source,destination)',1)[1].split('end subroutine move_block_storage',1)[0]
    if part.count('call move_alloc(')!=36: raise ValueError('storage move count changed')

    name='parallel_block_rhs_distribute.f90'; use(name)
    before(name,'    deallocate(packet_list,destination)',
           record('rhs.geometry','0','[int(used,int64),int(sum(sc),int64),int(sum(rc),int64), &\n         int(sc(rank+1),int64),int(rc(rank+1),int64)]',
                  'selection_pack_packets all_packed_bytes all_incoming_bytes self_packed_bytes self_incoming_bytes',
                  'int(sc,int64)','int(rc,int64)'))

    name='parallel_block_transform_placement.f90'; use(name)
    after(name,"    if(ierr/=MPI_SUCCESS) error stop 'transform metadata exchange'",
          record('transform.geometry','0','[int(sum(sc),int64),int(sum(rc),int64), &\n         int(sc(rank+1),int64),int(rc(rank+1),int64),int(size(input%scalar),int64), &\n         int(size(input%inner),int64),int(size(input%outer),int64),2_int64]',
                 'all_packed_bytes all_incoming_bytes self_packed_bytes self_incoming_bytes source_scalar_ops source_inner_ops source_outer_ops emission_passes',
                 'int(sc,int64)','int(rc,int64)'))

    name='parallel_block_rhs_exchange.f90'; use(name)
    after(name,"    call mpi_check(comm,ierr,'metadata delivery')",
          record('directory.deliver','width','[int(size(record,2),int64),int(size(received,2),int64)]',
                 'sent_records received_records','4_int64*int(width,int64)*int(sc,int64)',
                 '4_int64*int(width,int64)*int(rc,int64)'))
    after(name,"    call mpi_check(comm,ierr,'indexed metadata reply')",
          record('metadata.indexed','phase','[int(nq,int64),int(u,int64),int(nc,int64), &\n         int(size(reply),int64),int(size(received),int64),int(size(values),int64)]',
                 'requested_rows deduplicated_rows columns source_gather_values reply_buffer_values expanded_values',
                 '8_int64*int(nc,int64)*int(rc,int64)','8_int64*int(nc,int64)*int(sc,int64)'))
    before(name,'    call detail_enter(plan%profile_region)\n    call rhs_detail_enter(DP_RHS_PACK)',
           record('route.exchange','phase','[int(columns,int64), &\n         int(columns,int64)*int(sum(plan%sc),int64),int(columns,int64)*int(sum(plan%rc),int64), &\n         int(columns,int64)*int(plan%sc(plan%me+1),int64),int(columns,int64)*int(size(plan%target),int64), &\n         int(ns,int64)+int(nr,int64),int(size(plan%requests),int64), &\n         int(count(plan%sc>0)-merge(1,0,plan%sc(plan%me+1)>0),int64), &\n         int(count(plan%rc>0)-merge(1,0,plan%rc(plan%me+1)>0),int64)]',
                  'columns snapshot_values received_values self_copy_values installed_values active_scratch_values request_capacity remote_send_requests remote_receive_requests',
                  '8_int64*(int(columns,int64)*int(plan%sc,int64)+int(merge(1,0,plan%sc>0),int64))',
                  '8_int64*(int(columns,int64)*int(plan%rc,int64)+int(merge(1,0,plan%rc>0),int64))'))

    name='parallel_block_transform_adapter.f90'; use(name)
    after(name,'    integer :: i,j,k,d,p,e,r,q,nreq,ierr',
          record('transform.adapter','plan%tag','[int(size(plan%source,2),int64),int(size(plan%target),int64), &\n         int(plan%columns,int64)*int(sum(plan%rc),int64),int(plan%columns,int64)*int(plan%sc(rank+1),int64), &\n         int(plan%columns,int64)*int(size(plan%target),int64),int(size(plan%sent),int64)+int(size(plan%received),int64), &\n         int(size(plan%request),int64)]',
                 'source_rows target_rows gather_values self_copy_values installed_values buffer_values request_capacity',
                 '8_int64*int(plan%columns,int64)*int(plan%rc,int64)',
                 '8_int64*int(plan%columns,int64)*int(plan%sc,int64)'))

    name='parallel_block_transform_backend.f90'; use(name)
    after(name,'    call prepare_exchange_scratch\n    call compile_forward_order',
          record('transform.workspace','0','[int(nnode,int64),int(nrow,int64),int(nrow-nnode,int64), &\n         int(npublication,int64),int(np,int64),int(nboundary,int64),int(size(fixed_node),int64), &\n         int(nrow-nnode-size(save_node)-size(fixed_node),int64), &\n         int(size(scalar),int64)+int(size(vector),int64), &\n         int(size(scalar_wavelet),int64)+int(size(vector_wavelet),int64),int(size(exchange_scratch),int64), &\n         int(size(wire),int64)]',
                 'canonical_rows total_rows noncanonical_rows publication_requests domain_requests block_boundary_requests fixed_rows imported_publication_rows live_value_slots wavelet_value_slots exchange_scratch_values wire_values'))
    for name,text in changed.items():
        (src/name).write_text(text)
    report={'schema':'stage325-work-patch-v1','files':{name:hashlib.sha256(text.encode()).hexdigest() for name,text in changed.items()},
            'hooks':OPERATIONS,'records':SCHEMA,
            'limitations':['Integer counts and serialized-equivalent sizes are not elapsed attribution.',
                'Counts aggregate by tag and phase across the complete qualified lifecycle; no phase costs are added.',
                'Payload-only scratch excludes compiler descriptors, MPI internals, and transient heap peak.',
                'Diagnostic retained packed-size queries add qualification overhead; performance must disable ledger.',
                'The 17 direct migration collective count is source-certified; nested validation collectives are excluded.']}
    (root/'work325-schema.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');a=p.parse_args()
    print(json.dumps(apply(a.source),indent=2))
