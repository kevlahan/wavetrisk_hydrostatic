"""Serial numerical development gate; never a production ownership qualification.

Extract the original accepted step1 expressions for an independent primitive
oracle. Existing source/gradient operators are reused, not claimed independent.
No MPI, fixture rewriting, geometry-owner cache or production selector.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,time

ACCEPTED='ff740744608731a05d53771c32a3055309f3a586'
ROOT=Path(__file__).resolve().parents[3]

def sha(data): return hashlib.sha256(data).hexdigest()

def oracle(source):
    start=source.index('elseif (itype == 0) then ! standard')
    start=source.index('         u_dual_RT    =',start)
    stop=source.index('         ! Interpolate geopotential',start)
    code=source[start:stop]
    sw_start=source.index('         u_prim_RT_W  =',source.index('      elseif (itype == 0) then ! standard',stop))
    sw_stop=source.index('         ! Scalar fluxes',sw_start)
    sw_code=source[sw_start:sw_stop]
    return '''module legacy_primitive_mod
  use kind_mod, only : dp
  implicit none
  type :: Array
     real(dp) :: elts(21)=0.0_dp
  end type
  type :: Area
     real(dp) :: part(6)=0.0_dp,hex_inv=0.0_dp
  end type
  type :: Areas
     type(Area) :: elts(7)
  end type
  type :: Geometry
     type(Array) :: len,pedlen,triarea,coriolis
     type(Areas) :: areas
  end type
contains
  pure function interp(a,b) result(c)
    real(dp), intent(in) :: a,b
    real(dp) :: c
    c=0.5_dp*(a+b)
  end function
  subroutine legacy_primitive(x,g,y,sw_x,sw_g,sw_y)
    real(dp), intent(in) :: x(16),g(37)
    real(dp), intent(out) :: y(11)
    real(dp), optional, intent(in) :: sw_x(3),sw_g(12)
    real(dp), optional, intent(out) :: sw_y(8)
    integer, parameter :: EDGE=3,TRIAG=2,RT=0,DG=1,UP=2,LORT=0,UPLT=1,S_MASS=1
    integer, parameter :: NORTH=1,EAST=2,SOUTH=3,WEST=4,NORTHEAST=5,SOUTHWEST=7
    integer, parameter :: id=0,idN=1,idE=2,idS=3,idW=4,idNE=5,idSW=6
    integer, parameter :: id_i=1,idN_i=2,idE_i=3,idS_i=4,idW_i=5,idNE_i=6,idSW_i=7
    integer, parameter :: vnode(10)=[0,1,4,0,3,6,4,0,2,3],vedge(10)=[1,1,1,2,2,2,2,3,3,3]
    integer :: j,a
    type(Geometry) :: dom
    real(dp) :: velo(21),qe(21),vort(21),ke(7),rho_dz(0:7,1),pv_mass(4)
    real(dp) :: u_dual_RT,u_dual_DG,u_dual_UP,u_dual_RT_W,u_dual_DG_SW,u_dual_UP_S
    real(dp) :: u_prim_RT,u_prim_RT_N,u_prim_RT_W,u_prim_DG,u_prim_DG_S,u_prim_DG_W
    real(dp) :: u_prim_UP,u_prim_UP_E,u_prim_UP_S,u_prim_DG_SW
    real(dp) :: u_prim_RT_SW,u_prim_UP_SW,circ_SW_LORT,circ_SW_UPLT,pv_SW_LORT,pv_SW_UPLT
    real(dp) :: circ_LORT,circ_UPLT,circ_W_LORT,circ_S_UPLT,pv_LORT,pv_UPLT,pv_W_LORT,pv_S_UPLT
    velo=0;qe=0;vort=0;ke=0;rho_dz(0:5,1)=x(1:6)
    do j=1,10
       a=EDGE*vnode(j)+vedge(j)
       velo(a)=x(6+j);dom%len%elts(a)=g(j)
    end do
    dom%pedlen%elts([1,2,3,13,20,12])=g(11:16)
    dom%areas%elts(1)%part(1)=g(17);dom%areas%elts(3)%part(3)=g(18);dom%areas%elts(6)%part(5)=g(19)
    dom%areas%elts(1)%part(2)=g(20);dom%areas%elts(6)%part(4)=g(21);dom%areas%elts(2)%part(6)=g(22)
    dom%areas%elts(5)%part(1)=g(23);dom%areas%elts(1)%part(3)=g(24);dom%areas%elts(2)%part(5)=g(25)
    dom%areas%elts(4)%part(2)=g(26);dom%areas%elts(3)%part(4)=g(27);dom%areas%elts(1)%part(6)=g(28)
    dom%coriolis%elts([1,2,9,8])=g(29:32)
    dom%areas%elts(1)%hex_inv=g(33);dom%triarea%elts([1,2,9,8])=g(34:37)
'''+code+'''
    y(1)=u_dual_RT*interp(rho_dz(0,1),rho_dz(EAST,1))
    y(2)=u_dual_DG*interp(rho_dz(0,1),rho_dz(NORTHEAST,1))
    y(3)=u_dual_UP*interp(rho_dz(0,1),rho_dz(NORTH,1))
    y(4:6)=qe(1:3);y(7)=ke(1);y(8:11)=vort([1,2,9,8])
    if (present(sw_y)) then
       rho_dz(SOUTHWEST,1)=sw_x(1)
       velo([19,21])=sw_x(2:3);dom%len%elts([19,21])=sw_g(1:2)
       dom%areas%elts(7)%part(1)=sw_g(3);dom%areas%elts(4)%part(3)=sw_g(4)
       dom%areas%elts(1)%part(5)=sw_g(5);dom%areas%elts(7)%part(2)=sw_g(6)
       dom%areas%elts(1)%part(4)=sw_g(7);dom%areas%elts(5)%part(6)=sw_g(8)
       dom%coriolis%elts([13,14])=sw_g(9:10);dom%triarea%elts([13,14])=sw_g(11:12)
'''+sw_code+'''
       sw_y(1)=u_dual_RT_W*interp(rho_dz(0,1),rho_dz(WEST,1))
       sw_y(2)=u_dual_DG_SW*interp(rho_dz(0,1),rho_dz(SOUTHWEST,1))
       sw_y(3)=u_dual_UP_S*interp(rho_dz(0,1),rho_dz(SOUTH,1))
       sw_y(4:6)=qe([13,20,12]);sw_y(7:8)=vort([13,14])
    end if
  end subroutine
end module
'''

def pentagon_oracle(source):
    import re
    body=source.split('  subroutine post_step1 (',1)[1].split('  end subroutine post_step1',1)[0]
    body=body[body.index('    ! Parts 4, 5'):]
    body='\n'.join(line for line in body.splitlines() if not re.match(
        r'\s*(?:id(?:[A-Z]+)?|id_rhodz(?:\([^)]*\))?|rho_dz\([^)]*\))\s*=',line))
    return '''
  subroutine legacy_pentagon(c,rho,velocity,length,parts,area,coriolis,qe,vort)
    integer, intent(in) :: c
    real(dp), intent(in) :: rho(7),velocity(3,7),length(3,7),parts(6,7),area(2,7),coriolis(2,7)
    real(dp), intent(inout) :: qe(21),vort(14)
    integer, parameter :: EDGE=3,TRIAG=2,RT=0,DG=1,UP=2,LORT=0,UPLT=1
    integer, parameter :: NORTH=1,EAST=2,SOUTH=3,WEST=4,NORTHEAST=5,SOUTHWEST=7
    integer, parameter :: IJPLUS=5,IPLUSJMINUS=6,IJMINUS=7,IMINUSJPLUS=8
    integer, parameter :: id=0,idN=1,idE=2,idS=3,idW=4,idNE=5,idSW=6
    type(Geometry) :: dom
    real(dp) :: velo(21),rho_dz(0:7),pv_mass(3)
    real(dp) :: circ_LORT,circ_UPLT,circ_S_UPLT,circ_SW_LORT,circ_SW_UPLT,circ_W_LORT
    real(dp) :: pv_LORT,pv_UPLT,pv_S_UPLT,pv_SW_LORT,pv_SW_UPLT,pv_W_LORT
    real(dp) :: u_prim_RT,u_prim_RT_N,u_prim_RT_SW,u_prim_RT_W
    real(dp) :: u_prim_DG_SW,u_prim_UP,u_prim_UP_S,u_prim_UP_SW
    integer :: n
    velo=reshape(velocity,[21]);dom%len%elts=reshape(length,[21])
    dom%triarea%elts(1:14)=reshape(area,[14]);dom%coriolis%elts(1:14)=reshape(coriolis,[14])
    do n=1,7
       dom%areas%elts(n)%part=parts(:,n)
    end do
    rho_dz=rho([1,2,3,4,5,6,1,7])
'''+body+'''
  end subroutine legacy_pentagon
'''

def run(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    original=subprocess.check_output(['git','show',ACCEPTED+':src/ops.f90'],cwd=ROOT)
    if (ROOT/'src/ops.f90').read_bytes()!=original: raise ValueError('Primitive oracle source is not accepted Stage266')
    (out/'legacy_primitive.f90').write_text(oracle(original.decode()).replace('end module\n',pentagon_oracle(original.decode())+'end module\n'))
    sources=[ROOT/'src/kind.f90',ROOT/'src/parallel_block_mass.f90',ROOT/'src/parallel_block_velocity.f90',
             ROOT/'src/parallel_block_rhs.f90',out/'legacy_primitive.f90',Path(__file__).with_name('unit.f90')]
    result={'accepted_source':ACCEPTED,'production_active':False,'ownership_qualified':False,
            'source_sha256':{str(p):sha(p.read_bytes()) for p in sources},'modes':{}}
    for mode,flags in [('optimized',['-O3']),('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'])]:
        build=out/mode;build.mkdir()
        command=['gfortran','-ffp-contract=off','-ffree-line-length-132',*flags,*map(str,sources),'-o','unit']
        start=time.monotonic()
        with (build/'build.log').open('w') as log: subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT,check=True)
        tests={}
        for case,expected in [('normal',0),('bad_column',1),('bad_stencil',1),('bad_rk',1),('missing_sw',1),('bad_sw',1),('bad_pentagon',1)]:
            run_command=['./unit']+([] if case=='normal' else [case])
            proc=subprocess.run(run_command,cwd=build,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            (build/(case+'.log')).write_text(proc.stdout)
            if expected==0 and proc.returncode!=0: raise ValueError(mode+' '+case+' failed: '+proc.stdout)
            marker={'bad_column':'block RHS column extents differ','bad_stencil':'block RHS primitive stencil invalid','bad_rk':'block RHS RK extents differ','missing_sw':'block RHS SW physics missing','bad_sw':'block RHS SW stencil invalid','bad_pentagon':'block RHS pentagon stencil invalid'}.get(case)
            if expected and (not proc.returncode or marker not in proc.stdout): raise ValueError('Expected guard did not reject '+case)
            tests[case]={'exit_code':proc.returncode,'log_sha256':sha(proc.stdout.encode())}
        result['modes'][mode]={'command':command,'tests':tests,'elapsed_seconds':time.monotonic()-start,
                               'binary_sha256':sha((build/'unit').read_bytes())}
    result['status']='passed'
    result['scope']='Numerical draft only. Original accepted step1 direct/comp_SW and all four post_step1 pentagon oracles; pressure/injection/diffusion formula checks and reused source/gradient/RK connectivity. No MPI routes, production cutover, adaptive science, full-step timing or scaling qualification.'
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args();run(args.output)
