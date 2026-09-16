"""Observer-only batch/call scopes. Every original numerical byte is recoverable."""
from pathlib import Path
import re
# Shared labels mean related responsibilities, not identical work or removable cost.
GROUPS = [
('inverse entry and seed', ['prepare_native_inverse','prepare_block_native_inverse_boundaries']),
('inverse driver remainder', ['activate_native_inverse','activate_block_native_inverse_transform','inverse_wavelet_transform','inverse_scalar_transform','inverse_velo_transform']),
('inverse scalar kernel', ['scalar_inverse203','prolong_block_scalar_level','reconstruct_block_scalar_level','prolong_scalar','reconstruct_scalar']),
('inverse outer kernel', ['native_inverse_outer','reconstruct_native_inverse_outer_vectors','reconstruct_outer_velo','reconstruct_velo_penta']),
('inverse inner kernel', ['inner_inverse203','reconstruct_block_inner_vector_level','reconstruct_inner_velo']),
('inverse gather and restore', ['native_inverse_gather']),
('inverse scatter and cache', ['native_inverse_scatter']),
('inverse final block install', ['install_inverse203']),
('inverse transfer local preparation/install', ['transfer_nodes']),
('inverse alias local preparation/install', ['native_inverse_boundary']),
('inverse consumer alias snapshot', ['capture_native_inverse_consumer_aliases']),
('inverse Domain publication', ['publish_native_inverse']),
('block ghost refresh', ['refresh_block_wav_coeff_ghosts','refresh_block_sol_ghosts','refresh_block_ghosts']),
('forward scalar wavelet kernel', ['compute_scalar_wavelets','validate_scalar_wavelets']),
('forward vector wavelet kernel', ['compute_velo_wavelets','compute_velo_wavelets_penta','validate_outer_vector_wavelets','validate_complete_block_vector_wavelets']),
('forward restriction', ['restrict_scalar','restrict_velo','validate_velocity_restriction']),
('compression', ['compress_wavelets','compress_wavelets_scalar','compress_wavelets_velo','prepare_native_compression','activate_native_compression']),
('Domain boundary local preparation/install', ['update_bdry','update_bdry1','update_bdry__start','update_bdry__finish','refresh_scalar_boundary','refresh_vector_boundary']),
('Domain to block import', ['import_domain_field_family_to_blocks','exchange_domain_to_block_payloads','import_domain_boundary_field_family_to_blocks']),
('block to Domain publication', ['write_block_field_family_to_domains','refresh_parallel_block_domain_prognostic_state','reconstruct_block_writeback_domain_stage','publish_native_mass_rhs']),
('scalar capture preparation', ['prepare_block_scalar_capture_addresses','seed_scalar_producer_geometry','begin_block_scalar_divergence_capture']),
('scalar capture transport/install', ['finalize_block_scalar_divergence_capture','capture_block_scalar_divergence_level']),
('scalar flux arithmetic/access', ['compute_block_scalar_direct_flux','cpt_or_restr_flux']),
('scalar restriction arithmetic/access', ['restrict_block_scalar_flux_level','cpt_or_restr_scalar']),
('scalar divergence arithmetic/access', ['recompute_block_scalar_divergence_level','cal_scalar_trend']),
('scalar ghost transport local work', ['exchange_block_scalar_restriction_ghosts','exchange_temperature_edges']),
('scalar pipeline driver remainder', ['evaluate_candidate_block_scalar_restriction']),
('mass plan preparation', ['prepare_native_mass']),
('mass transport local work', ['begin_native_mass','finish_native_mass']),
('velocity plan preparation', ['prepare_native_velocity_programs']),
('velocity tendency arithmetic/access', ['velocity_trend_source','velocity_trend_grad','compute_native_velocity_source','compute_native_velocity_gradient','activate_block_native_velocity_tendency']),
('basic operators remainder', ['basic_operators']),
('surface pressure', ['cal_surf_press']),
('diffusion preparation', ['cal_divu_ml','cal_Laplacian_scalars','cal_Laplacian_divu','cal_Laplacian_vector_rot']),
('RK retention/assembly', ['retain_block_native_multistage_candidate','RK_sub_step_compatibility','RK_sub_step','manage_q1_mem']),
('adaptation masks', ['init_masks_zero','mask_active','mask_adj_same_scale','mask_restrict_same_scale','mask_adj_finer_scale','complete_masks','mask_second_neighbours']),
('adaptation topology update', ['update','post_refine','fill_up_level']),
('block catalog/ownership', ['build_parallel_block_catalog','build_block_catalog']),
('block source construction', ['build_source_blocks','build_one_source_block']),
('block migration remainder', ['migrate_blocks']),
('migration routes and storage', ['build_block_migration_manifest','exchange_block_migration_sizes','clear_block_migration_manifest','clear_block_staging']),
('migration install/retained state', ['install_local_blocks','prepare_parallel_block_production_state']),
('adaptation reconstruction retention', ['retain_post_grid_change_block_reconstruction']),
('adaptation lifecycle remainder', ['complete_block_adaptation_lifecycle']),
('Domain storage initialization', ['zero_float']),
('MPI post', ['MPI_Isend','MPI_Irecv']),
('MPI completion', ['MPI_Waitall','MPI_Wait','MPI_Test','MPI_Testall']),
('MPI collective/blocking', ['MPI_Alltoallv','MPI_Alltoall','MPI_Allreduce','MPI_Reduce','MPI_Bcast','MPI_Barrier','MPI_Allgather','MPI_Allgatherv','MPI_Exscan','MPI_Send','MPI_Recv','MPI_Sendrecv']),
]
GROUPS.extend([
('inverse transfer/alias pack', []),('inverse ordered local aliases', []),('inverse transfer/alias install', []),
('Domain boundary pack', []),('Domain boundary route construction', []),('Domain boundary buffer management', []),
('Domain boundary local copies', []),('Domain boundary receive installation', [])])
DETAIL={'parallel_block_inverse.f90':{'DP_INVERSE_PACK':49,'DP_INVERSE_LOCAL':50,'DP_INVERSE_INSTALL':51},
        'comm_mpi.f90':{'DP_BDRY_PACK':52,'DP_BDRY_ROUTE':53,'DP_BDRY_STORAGE':54,'DP_BDRY_LOCAL':55,'DP_BDRY_INSTALL':56}}
GROUPS.extend([
('native mass arithmetic', ['execute_mass_divergence','execute_mass_restriction']),
('native velocity publication', ['prepare_block_native_velocity_remainder']),
('RK candidate capture/staging', ['capture_block_domain_multistage_candidate_tendency','begin_block_domain_multistage_candidate_stage',
 'prepare_block_native_multistage_wavelet_stage','prepare_block_native_multistage_wavelet_acceptance','refresh_parallel_block_candidate_boundary_state']),
('temperature boundary contract', ['compute_temperature_rk_boundary','initialize_scalar_restriction_boundary_flux','prepare_temperature_boundary_stage','apply_temperature_boundary_stage']),
('physics consumer', ['apply_domain_physics_consumers','physics_simple_step']),
('vertical remap', ['remap_vertical_coordinates']),
('output/diagnostic consumer', ['prepare_end_step_domain_consumers','cal_total_mass','write_output_consumer','write_checkpoint']),
('grid-change compatibility preparation', ['prepare_parallel_block_grid_change','update_domain_soil_wavelet_consumer']),
('tendency driver remainder', ['block_tendency_compatibility_ml','trend_ml'])])
REGIONS={i+1:label for i,(label,_) in enumerate(GROUPS)}
TARGETS={n.lower():i+1 for i,(_,names) in enumerate(GROUPS) for n in names}

def instrument(root):
    root=Path(root);inventory={}
    for p in sorted((root/'src').glob('*.f90')):
        if p.is_symlink() or p.name.startswith(('workflow_','parallel_block_profile')):continue
        original=p.read_text();lines=original.splitlines(True);out=[];edits=[];sites=[];i=0
        while i<len(lines):
            start=i;i+=1
            while lines[i-1].split('!')[0].rstrip().endswith('&') and i<len(lines):i+=1
            old=''.join(lines[start:i]);logical=' '.join(l.split('!')[0].strip().strip('&') for l in lines[start:i])
            m=re.match(r'\s*(?:if\s*\(.*\)\s*)?call\s+(\w+)\b\s*(?:\(\s*(\w+))?',logical,re.I)
            target=None
            if m:
                name=m[1].lower();arg=(m[2] or '').lower()
                if name in TARGETS:target=name
                # Match domain/layer or all-block batches, never point accessors.
                elif name in ('apply_interscale_d','apply_interscale_d2','apply_to_penta_d','apply_onescale_d',
                              'apply_local_block_field_producer','apply_local_block_field_consumer','call_adaptation_domain_consumer') and arg in TARGETS:target=arg
            if target:
                tag=TARGETS[target];indent=re.match(r'[ \t]*',old)[0]
                new=indent+f'call measure_enter204({tag})\n'+old+indent+f'call measure_leave204({tag})\n'
                edits.append((new,old));out.append(new);sites.append({'line':start+1,'target':target,'operation':tag})
            else:out.append(old)
        changed=''.join(out)
        def hook(m):
            key=m[3].upper();tag=DETAIL.get(p.name,{}).get(key)
            if tag is None:return m.group()
            indent=m[1];name=m[2].lower();old=m.group()
            new=(indent+f'call measure_enter204({tag})\n'+old if name=='enter' else old+'\n'+indent+f'call measure_leave204({tag})')
            edits.append((new,old));sites.append({'target':key,'operation':tag,'scope':'batch detail boundary'})
            return new
        changed=re.sub(r'^([ \t]*)call detail_(enter|leave)\((\w+)\)[ \t]*$',hook,changed,flags=re.M|re.I)
        if not sites:continue
        m=re.search(r'^module \w+[^\n]*\n',changed,re.M|re.I)
        if not m:raise ValueError('Missing module import site '+p.name)
        new=m.group()+'  use workflow_profile_195_mod, only: measure_enter204,measure_leave204\n'
        edits.append((new,m.group()));changed=changed[:m.start()]+new+changed[m.end():]
        restored=changed
        for new,old in reversed(edits):restored=restored.replace(new,old)
        if restored!=original:raise ValueError('Cannot restore original numerical bytes '+p.name)
        p.write_text(changed);inventory[str(p.relative_to(root))]={'sites':sites,'original_bytes_recoverable':True}
    seen={x['operation'] for f in inventory.values() for x in f['sites']}
    if not {2,18,36,46,47}.issubset(seen):raise ValueError('Incomplete core instrument sites '+str(seen))
    return inventory
