#!/usr/bin/env python3
"""Checked full-solver regression: one-Domain restart source, empty receiver.

All ownership changes and metadata exports are confined to a disposable tree.
Source 0/empty rank 3 and source 1/empty rank 0 use the real block distributor.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace
import local_j4
from capability_223 import coverage, events, numerical, save
from candidate_capability_223 import once
from run_local_220 import SELECTORS

CHECK = '''
  subroutine stage224_export_prefix(fid,prefix_bytes,id)
    integer,intent(in) :: fid,id
    integer(int64),intent(in) :: prefix_bytes
    integer :: unit
    integer(int8),allocatable :: expected(:)
    character(128) :: filename

    if(rank==owner(1)) then
       if(size(grid)/=1) error stop 'Stage224 source must own exactly one Domain'
       allocate(expected(int(prefix_bytes)))
       ! Independent expected bytes at the start of the serialized record.
       read(fid,pos=1) expected
       write(filename,'(a,i0,a)') 'prefix224-expected-',id,'.bin'
       open(newunit=unit,file=trim(filename),status='new',form='unformatted',access='stream')
       write(unit) expected
       close(unit)
    end if
    write(filename,'(a,i0,a,i0,a)') 'prefix224-reload-',id,'-rank-',rank,'.bin'
    open(newunit=unit,file=trim(filename),status='new',form='unformatted',access='stream')
    write(unit) istep_cumul
    write(unit) time
    write(unit) itime
    write(unit) iwrite
    write(unit) threshold
    call dump(unit)
    close(unit)
  end subroutine stage224_export_prefix

'''


def instrument(tree, checkpoint):
    p=tree/'src/arch.f90'
    s=once(p.read_text(), "    case('native')", """    case('single-source0','single-source1')
       if(n_process/=4) error stop 'Stage224 local regression requires four ranks'
       if(rank==0) then
          if(trim(layout223)=='single-source0') then
             owner(1)=0
             do d=2,n_glo_block
                owner(d)=1+mod(d-2,2)
             end do
          else
             owner(1)=1
             do d=2,n_glo_block
                owner(d)=2+mod(d-2,2)
             end do
          end if
       end if
    case('native')""")
    p.write_text(s)
    s=checkpoint.read_text()
    # Support both the original Stage223 negative control and the corrected source.
    if 'use arch_mod, only : stage223_event' not in s:
        s=once(s,'  use mpi_f08','  use mpi_f08\n  use arch_mod, only : stage223_event')
        call='         call distribute_checkpoint_global_state(fid,global_prefix_position,global_prefix_bytes)'
        s=once(s,call,call+"\n    call stage223_event('restart',size(grid),-1)")
    s=once(s,"    call stage223_event('restart',size(grid),-1)",
           "    call stage224_export_prefix(fid,global_prefix_bytes,id)\n    call stage223_event('restart',size(grid),-1)")
    s=once(s,'  subroutine read_checkpoint_global_state (fid, position)',CHECK+'  subroutine read_checkpoint_global_state (fid, position)')
    (tree/'src/checkpoint.f90').write_text(s)


def run(a):
    out=a.out.resolve();out.mkdir(exist_ok=False)
    tree=out/'tree'
    shutil.copytree(a.prototype,tree,symlinks=True,ignore=shutil.ignore_patterns('build-*','bin-*','__pycache__'))
    instrument(tree,a.checkpoint)
    if getattr(a,'mpi_source',None):
        text=a.mpi_source.read_text()
        text=once(text,'module parallel_block_mpi_mod\n  use parallel_block_scalar_storage_mod',
                  'module parallel_block_mpi_mod\n  use arch_mod, only : stage223_event\n  use parallel_block_scalar_storage_mod')
        text=once(text,'    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)',
                  "    call install_local_blocks(size(block_catalog),local_seen,transfer_storage)\n    call stage223_event('install',size(grid),n_local_blocks())")
        (tree/'src/parallel_block_mpi.f90').write_text(text)
    check_optimization=getattr(a,'check_optimization','O0')
    if check_optimization!='O0':
        with (tree/'local-flags.mk').open('a') as flags:
            # GCC 15 warns about an allocatable descriptor allocated on pass 1
            # of a fixed two-pass loop. Keep the warning visible and all runtime
            # checks enabled; do not make this optimization-only warning fatal.
            flags.write('\nFFLAGS += -'+check_optimization+' -Wno-error=maybe-uninitialized -Wno-error=function-elimination\n')
    command=['make','-f','Makefile','-f','local-flags.mk','-j1','TEST_CASE=climate','PARAM=param_J4',
             'DEBUG=check','FP_CONTRACT=off','STAGE202=on','STAGE203=off','BUILD_DIR=build-224-checked','BIN_DIR=bin-224-checked']
    result=dict(status='building',checkpoint_source=str(a.checkpoint),
                checkpoint_sha256=hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),build_command=command,
                negative_control=a.negative_control,oracles=getattr(a,'oracles',False),check_optimization=check_optimization,cases={})
    if getattr(a,'mpi_source',None):
        result['mpi_source']=str(a.mpi_source)
        result['mpi_source_sha256']=hashlib.sha256(a.mpi_source.read_bytes()).hexdigest()
    save(out/'regression.json',result)
    if a.built_tree:
        previous=json.loads((a.built_tree.parent/'regression.json').read_text())
        assert previous.get('check_optimization','O0')==check_optimization, 'Checked compiler options differ'
        for path in (tree/'src').rglob('*'):
            if path.is_file() and path.suffix.lower()=='.f90':
                assert path.read_bytes()==(a.built_tree/path.relative_to(tree)).read_bytes(),path
        (tree/'bin-224-checked').mkdir()
        shutil.copy2(a.built_tree/'bin-224-checked/climate',tree/'bin-224-checked/climate')
        shutil.copy2(a.built_tree.parent/'build.log',out/'build.log')
        result['reused_checked_build']=str(a.built_tree)
    else:
        print('Building checked regression: '+str(out/'build.log'),flush=True)
        with (out/'build.log').open('w') as log:subprocess.run(command,cwd=tree,stdout=log,stderr=subprocess.STDOUT,check=True)
    result['binary_sha256']=hashlib.sha256((tree/'bin-224-checked/climate').read_bytes()).hexdigest()
    save(out/'regression.json',result)
    result['status']='running';save(out/'regression.json',result)
    original=local_j4.environment
    try:
        for layout,source,empty in [('single-source0',0,3),('single-source1',1,0)]:
            def env(**kw):
                v=original(**kw);v.update(SELECTORS,WAVETRISK_STAGE223_ALLOW='1',WAVETRISK_STAGE223_LAYOUT=layout)
                return v
            local_j4.environment=env
            case=out/layout;print('Running '+layout,flush=True)
            try:
                local_j4.run(SimpleNamespace(build=tree,fixture=a.fixture,out=case,ranks=4,
                                            binary='bin-224-checked/climate',oracles=getattr(a,'oracles',False),profile=False))
            except RuntimeError:
                log=(case/'run.log').read_text()
                failing_source=(tree/'src/checkpoint.f90').read_text().splitlines()
                expected_line=next(i for i,line in enumerate(failing_source,1) if 'read(fid,pos=domain_pos(1)) prefix' in line) if a.negative_control else -1
                if (not a.negative_control or 'Fortran runtime error: End of file' not in log or
                    f'At line {expected_line} of file src/checkpoint.f90' not in log):raise
                result['cases'][layout]=dict(expected_failure='EOF in distribute_checkpoint_global_state')
                save(out/'regression.json',result)
                continue
            if a.negative_control:raise AssertionError('Negative control unexpectedly passed')
            identity=json.loads((case/'results.json').read_text())
            if getattr(a,'oracles',False):
                log=(case/'run.log').read_text()
                assert 'Block dynamics validation oracle enabled' in log
                assert 'Block adaptation validation oracle enabled' in log
                assert 'Block vertical-remap validation oracle enabled' in log
                assert all(identity['environment'][f'WAVETRISK_VALIDATE_BLOCK_{name}']=='1' for name in ('DYNAMICS','ADAPTATION','REMAP'))
            analysis,cov=coverage(case,'J4',identity['parameters']['run_id'])
            life=events(case,4,40,len(analysis['step_records']),'empty-root' if empty==0 else 'native')
            assert life['restart_source_ranks']==[source]
            assert all(z==[empty] for z in life['zero_domain_rank_sets'])
            for line in (case/f'capability223-rank-{source}.txt').read_text().splitlines():
                assert int(line.split()[3])==1,'Source owns more than one Domain'
            prefix={}
            for index in cov['reloads']:
                expected=(case/f'prefix224-expected-{index}.bin').read_bytes()
                for rank in range(4):
                    assert (case/f'prefix224-reload-{index}-rank-{rank}.bin').read_bytes()==expected, (layout,index,rank)
                prefix[str(index)]=dict(bytes=len(expected),sha256=hashlib.sha256(expected).hexdigest(),all_ranks_exact=True)
            gate=numerical(a.reference,case,cov['new_checkpoint'],40,False)
            assert gate['passed'],gate['failures']
            reference=json.loads((a.reference/'results.json').read_text())['analysis']
            assert [x['state'] for x in analysis['step_records']]==[x['state'] for x in reference['step_records']]
            result['cases'][layout]=dict(coverage=cov,events=life,prefix=prefix,numerical=gate)
            save(out/'regression.json',result)
            print('PASS '+layout+': exact restart metadata on all ranks; full adaptive/restart trajectory',flush=True)
    except BaseException as error:
        result.update(status='failed',error=repr(error));save(out/'regression.json',result)
        raise
    finally:local_j4.environment=original
    result['status']='passed';save(out/'regression.json',result)
    print('PASS: '+str(out/'regression.json'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('prototype','checkpoint','fixture','reference','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--negative-control',action='store_true')
    p.add_argument('--oracles',action='store_true')
    p.add_argument('--check-optimization',choices=['O0','O1'],default='O0',help='Retain all DEBUG=check checks with selected optimization')
    p.add_argument('--mpi-source',type=Path,help='Overlay current production MPI source with test lifecycle evidence')
    p.add_argument('--built-tree',type=Path,help='Reuse a checked build only after matching every Fortran source')
    run(p.parse_args())
