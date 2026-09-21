#!/usr/bin/env python3
"""Compile the real architecture guards and production configuration predicate.

Run with MPI socket access. Two ranks and a one-Domain stub isolate the actual
rank-limit policy; no full climate mesh or checkpoint is needed for rejection.
"""
import json
from pathlib import Path
import subprocess
import tempfile
import argparse


def run(repo,out):
    out.mkdir(exist_ok=False)
    arch=(repo/'src/arch.f90').read_text()
    main=(repo/'src/main.f90').read_text()
    start=main.index('  logical function domainless_block_configuration_supported')
    end=main.index('  end function domainless_block_configuration_supported',start)+len('  end function domainless_block_configuration_supported')
    predicate=main[start:end]
    (out/'stub.f90').write_text('''module kind_mod
use iso_fortran_env, only: dp=>real64
end module
module shared_mod
integer,parameter :: N_GLO_DOMAIN=1,N_CHDRN=4
integer,allocatable :: n_domain(:)
character(255) :: run_id='test'
end module
module policy_mod
character(255) :: test_case='climate',timeint_type='RK4',physics_type='Simple'
integer :: resume=2
logical :: compressible=.true.,mode_split=.false.,physics_model=.true.,NCAR_topo=.false.,sso=.false.
contains
'''+predicate+'''
end module
''')
    (out/'arch.f90').write_text(arch)
    (out/'driver.f90').write_text('''program guard
use arch_mod
use policy_mod
implicit none
character(32) :: mode
call get_command_argument(1,mode)
if(trim(mode)=='legacy') then
 call init_arch_mod
else
 call init_arch_mod(defer_domain_limit=.true.)
 select case(trim(mode))
 case('cold');resume=-1
 case('split');mode_split=.true.
 case('incompressible');compressible=.false.
 case('rk3');timeint_type='RK3'
 case('physics');physics_type='Held_Suarez'
 case('no-physics');physics_model=.false.
 case('topography');NCAR_topo=.true.
 case('sso');sso=.true.
 case('other-driver');test_case='ocean'
 end select
 call validate_domain_rank_count(domainless_block_configuration_supported())
 if(trim(mode)=='actual-capacity') call check_parallel_task_capacity(1)
 if(trim(mode)=='boundary') call check_parallel_task_capacity(2)
 call distribute_grid(-1)
 if(rank==1.and.size(glo_id,2)/=1) error stop 'Domain distribution changed'
 if(rank==0) print *, 'PASS supported production guard'
end if
call finalize
end program
''')
    cmd=['mpif90','-ffree-line-length-132','-Wall','-Wextra','-Werror','-fcheck=all','stub.f90','arch.f90','driver.f90','-o','guard']
    subprocess.run(cmd,cwd=out,check=True,capture_output=True,text=True)
    cases={}
    modes=['supported','boundary','legacy','cold','split','incompressible','rk3','physics','no-physics','topography','sso','other-driver','actual-capacity','upper-capacity']
    for mode in modes:
        n=5 if mode=='upper-capacity' else 2
        r=subprocess.run(['mpirun','--oversubscribe','--bind-to','none','-n',str(n),'./guard',mode],cwd=out,capture_output=True,text=True,timeout=45)
        log=r.stdout+r.stderr;(out/(mode+'.log')).write_text(log)
        passes=mode in ('supported','boundary')
        assert (r.returncode==0)==passes,(mode,log)
        marker='PASS supported production guard' if passes else ('exceeds available parallel blocks' if mode in ('actual-capacity','upper-capacity') else 'More MPI tasks than Domains requires restarted')
        assert marker in log,(mode,log)
        cases[mode]=dict(ranks=n,returncode=r.returncode,expected='accept' if passes else 'reject')
    result=dict(status='passed',cases=cases,build_command=cmd)
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 14 compiled production guard cases')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.repo.resolve(),a.out.resolve())
