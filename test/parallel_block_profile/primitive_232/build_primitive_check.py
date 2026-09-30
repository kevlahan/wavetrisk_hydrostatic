"""Build a test-only comparison of explicit primitive kernels against live arithmetic.

This is not a block-owner production implementation or a performance campaign.
"""
from pathlib import Path
import argparse,re,shutil,subprocess,sys
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import save,digest

def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Expected one source anchor: '+old[:100])
    return text.replace(old,new)


def routine(text,name,edit):
    start=re.search(r'^\s*(?:subroutine|function)\s+'+re.escape(name)+r'\b',text,re.M|re.I)
    if start is None:raise ValueError('Missing routine '+name)
    end=re.search(r'^\s*end\s+(?:subroutine|function)\s+'+re.escape(name)+r'\b',
                  text[start.end():],re.M|re.I)
    if end is None:raise ValueError('Missing end '+name)
    stop=start.end()+end.end()
    return text[:start.start()]+edit(text[start.start():stop])+text[stop:]


def build(reference,out):
    reference,out=reference.resolve(),out.resolve()
    shutil.copytree(reference,out,symlinks=True)
    shutil.copy2(BASE/'src/parallel_block_primitive_kernels.f90',out/'src')
    shutil.copy2(Path(__file__).with_name('primitive_oracle_232.f90'),out/'src')
    p=out/'Makefile';s=p.read_text()
    s=once(s,'      parallel_block_profile.f90 \\\n', '      parallel_block_primitive_kernels.f90 \\\n      primitive_oracle_232.f90 \\\n      parallel_block_profile.f90 \\\n')
    s+='''
$(BUILD_DIR)/parallel_block_primitive_kernels.o: $(BUILD_DIR)/kind.o
$(BUILD_DIR)/primitive_oracle_232.o: $(BUILD_DIR)/domain.o $(BUILD_DIR)/parallel_block_primitive_kernels.o
$(BUILD_DIR)/ops.o $(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/test.o $(BUILD_DIR)/test_case_module.o: $(BUILD_DIR)/primitive_oracle_232.o
'''
    p.write_text(s)
    p=out/'src/ops.f90';s=p.read_text().replace('module ops_mod','module ops_mod\n  use primitive_oracle_232_mod',1)
    # Check results only after the original expression has executed. q/mean and
    # incoming context are immutable through these selected operations.
    def step(t):
        a='''      end if
    end subroutine comput'''
        b='''      end if
      if (itype==0.and.compressible) call check_regular232(dom, &
           [id,idN,idE,idS,idW,idNE,id,idSW],q(S_MASS,zlev)%data(d)%elts,velo,zlev, &
           [pv_W_LORT,pv_S_UPLT],[u_dual_RT,u_dual_DG,u_dual_UP])
    end subroutine comput'''
        t=once(t,a,b)
        a='''      end if
    end subroutine comp_SW'''
        b='''      end if
      if (itype==0.and.compressible) call check_southwest232(dom, &
           [id,id,id,idS,idW,id,id,idSW],q(S_MASS,zlev)%data(d)%elts,velo,zlev,[pv_W_LORT,pv_S_UPLT])
    end subroutine comp_SW'''
        return once(t,a,b)
    s=routine(s,'step1',step)
    def pentagon(t):
        import re
        pattern=r'(       rho_dz\([^\n]+ = mass\([^\n]+\n)'
        def hook(m):
            line=m[0]
            if '1:5' in line:
                ids='[id,idN,idE,idS,idW,id,id,id]'
            else:ids='id_rhodz'
            return line+f'       if(compressible) call begin_pentagon232(c,dom,{ids},mass,velo,mean_m)\n'
        t,n=re.subn(pattern,hook,t)
        assert n==4,n
        # Exactly four outer corner branches; no nested branches in this routine.
        assert t.count('    end if')==4
        return t.replace('    end if','       if(compressible) call end_pentagon232(dom)\n    end if')
    s=routine(s,'post_step1',pentagon);p.write_text(s)
    p=out/'src/diagnostics.f90';s=p.read_text().replace('module diagnostics_mod','module diagnostics_mod\n  use primitive_oracle_232_mod',1)
    def pressure(t):
        t=once(t,'    dom%geopot_lower%elts(id) = dom%geopot%elts(id)', '''    if(compressible) call begin_pressure232(mass(id),mean_m(id),temp(id),mean_t(id), &
         grav_accel,c_p,p_0,kappa,dom%press_lower%elts(id),dom%geopot%elts(id))
    dom%geopot_lower%elts(id) = dom%geopot%elts(id)''')
        return once(t,'    dom%press_lower%elts(id) = p_upper','''    dom%press_lower%elts(id) = p_upper
    if(compressible) call end_pressure232(dom%press_lower%elts(id),dom%geopot%elts(id), &
         dom%press%elts(id),exner(id),dom%geopot_lower%elts(id))''')
    s=routine(s,'integrate_pressure_up',pressure)
    s=routine(s,'cal_surf_press',lambda t:once(t,'  end subroutine cal_surf_press',
           '    if(compressible) call check_surface232(q)\n  end subroutine cal_surf_press'))
    p.write_text(s)
    p=out/'test/climate/test_case_module.f90';s=p.read_text()
    s=s.replace('module test_case_mod','module test_case_mod\n  use primitive_oracle_232_mod',1)
    def scalar(t):
        a='       flux = (-1)**Laplace_sclr * nu_scale (v, zlev) * grad * l_e'
        b=a+'''
       if(Laplace_sclr==1) then
          call check_scalar_diffusion232(q(v,zlev)%data(d)%elts([id,idE,idNE,idN]+1), &
               d_e,l_e,Laplace_sclr,nu_scale(v,zlev),flux)
       else if(Laplace_sclr==2) then
          call check_scalar_diffusion232(Laplacian_scalar(v)%data(d)%elts([id,idE,idNE,idN]+1), &
               d_e,l_e,Laplace_sclr,nu_scale(v,zlev),flux)
       end if'''
        return once(t,a,b)
    s=routine(s,'physics_scalar_flux_case',scalar)
    def velocity(t):
        # Restrict the test adapter to the active 2/2 climate workload. Other
        # selectors may leave the unused global pointers unassociated.
        a='    if (penta) C_visc(S_DIVU,1:zlevels) = C_temp'
        b='''    if(Laplace_divu==2.and.Laplace_rotu==2) then
       call check_velocity_diffusion232(divu([id,idE,idNE,idN]+1), &
            [vort(TRIAG*id+LORT+1),vort(TRIAG*id+UPLT+1), &
             vort(TRIAG*idx(i,j-1,offs,dims)+UPLT+1),vort(TRIAG*idx(i-1,j,offs,dims)+LORT+1)], &
            dom%len%elts(EDGE*id+1:EDGE*id+3),dom%pedlen%elts(EDGE*id+1:EDGE*id+3), &
            Laplace_divu,Laplace_rotu,nu_scale(S_DIVU,zlev),nu_scale(S_ROTU,zlev),source)
    end if
'''+a
        return once(t,a,b)
    s=routine(s,'physics_velo_source_case',velocity);p.write_text(s)
    p=out/'test/climate/climate.f90';s=p.read_text().replace('program climate','program climate\n  use primitive_oracle_232_mod, only: report232',1).replace('  call finalize','  call report232\n  call finalize');p.write_text(s)
    command=['make','-j1','TEST_CASE=climate','PARAM=param_J4','DEBUG=check','FP_CONTRACT=off',
             'STAGE202=on','STAGE203=on','BUILD_DIR=build-230-checked','BIN_DIR=bin-230-checked']
    with (out/'build-232.log').open('w') as f:
        code=subprocess.run(command,cwd=out,stdout=f,stderr=subprocess.STDOUT).returncode
    if code:raise RuntimeError(f'build failed: {out}/build-232.log')
    changed=['src/ops.f90','src/diagnostics.f90','src/parallel_block_primitive_kernels.f90',
             'src/primitive_oracle_232.f90','test/climate/climate.f90','test/climate/test_case_module.f90','Makefile']
    save(out/'build-232.json',dict(status='passed',test_only=True,command=command,
         sources={n:digest(out/n) for n in changed},binary_sha256=digest(out/'bin-230-checked/climate')))
    print('PASS test-only primitive arithmetic build',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    build(**vars(p.parse_args()))
