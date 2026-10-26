"""Single exact live boundary-import routine hook, shared by final patch/delta."""
SCHEMA={
 'boundary.scalar':{'fields':['remote_gather_values','remote_installed_values','verification_values'],'peer_unit':'bytes'},
 'boundary.vector':{'fields':['remote_gather_values','remote_installed_values','verification_values'],'peer_unit':'bytes'},
 'boundary.local':{'fields':['local_boundary_records','local_scalar_gather_install_values','local_vector_gather_install_values',
                            'verification_values'],'peer_unit':None},
}

def patch(text):
    start=text.index('  subroutine import_domain_boundary_field_family_to_blocks ( &')
    end=text.index('  end subroutine import_domain_boundary_field_family_to_blocks',start)
    old=text[start:end];body=old
    def change(a,b):
        nonlocal body
        if body.count(a)!=1:raise ValueError('boundary325.hook_anchor '+repr(a[:60]))
        body=body.replace(a,b)
    a='    logical :: transfer_vector\n'
    change(a,a+'    integer(int64) :: work325_boundary_local(3)\n')
    a='    call block_profile_enter(BLOCK_PROFILE_BOUNDARY_GHOST)'
    change(a,a+'\n    if(work325_enabled()) work325_boundary_local=0_int64')
    a='          call pack_domain_boundary_prognostic( &\n               d,source_bdry,n_node,payload_family, &\n               ghost_exchange_plan%scalar_patch_buffer(1:n_scalar), &'
    change(a,"          if(work325_enabled()) work325_boundary_local=work325_boundary_local+ &\n               [1_int64,int(merge(n_scalar,0,transfer_scalar),int64), &\n               int(merge(n_vector,0,transfer_vector),int64)]\n"+a)
    a='    call record_parallel_block_profile_volume( &\n         BLOCK_PROFILE_BOUNDARY_GHOST, &'
    code=''
    for component in ('scalar','vector'):
        prefix='block_writeback_plan%boundary_'+component
        sent=prefix+'_domain_send_count';recv=prefix+'_block_recv_count'
        code+=f"""    if(work325_enabled().and.transfer_{component}) call work325_record('boundary.{component}',payload_family, &
         [int(sum({sent}),int64), &
         int(sum({recv}),int64), &
         int(merge(sum({recv}),0,check_current),int64)], &
         8_int64*int({sent},int64), &
         8_int64*int({recv},int64))
"""
    code+="""    if(work325_enabled()) call work325_record('boundary.local',payload_family, &
         [work325_boundary_local,merge(sum(work325_boundary_local(2:)),0_int64,check_current)])
"""
    change(a,code+a)
    return text[:start]+body+text[end:]
