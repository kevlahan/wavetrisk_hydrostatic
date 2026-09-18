#!/usr/bin/env python3
"""Validate final J4 candidate: optimized off/on, route-audit-check off, restart and legacy gates."""
import argparse,json
from pathlib import Path
from local_j4 import compare_windows
from cluster_199 import workload_gate
from compare_restart_checkpoints import compare
from checkpoint_gate import accept
from build_restriction_208 import verify,digest
from candidate_boundary_215 import candidate
import tempfile,shutil

def validate(a):
    m=json.loads((a.parent_kit/'stage208-sources.json').read_text())
    ref=json.loads((a.reference/'results.json').read_text());legacy=json.loads((a.legacy/'results.json').read_text())
    if ref['status']!='passed' or legacy['status']!='passed' or ref['ranks']!=4 or legacy['ranks']!=4:raise ValueError('Invalid references')
    verify(Path(ref['build']),m['reference_sources']);verify(Path(legacy['build']),m['legacy_io_sources'])
    if digest(a.reference/'climateJ4')!=ref['binary_sha256'] or digest(a.legacy/'climateJ4')!=legacy['binary_sha256']:raise ValueError('Reference binary changed')
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);(root/'src').mkdir()
        names=['src/comm_mpi.f90']
        for name in names:shutil.copy2(Path(ref['build'])/name,root/name)
        candidate(root);sources=dict(m['reference_sources'])
        for n in ('test/climate/test_case_module.f90','test/climate/climate.f90'):sources[n]=digest(Path(ref['build'])/n)
        for name in names:sources[name]=digest(root/name)
    result=dict(stage=215,status='passed',candidate_sources=sources,runs={},
                scope='Local J4/J6 four-rank correctness only; no performance or J5 claim.',
                reference=dict(path=str(a.reference),binary_sha256=ref['binary_sha256']))
    for name in ('optimized','checked-off','oracles-on'):
        root=getattr(a,name.replace('-','_'));r=json.loads((root/'results.json').read_text())
        if r['status']!='passed' or r['ranks']!=4 or r['reloads']!=[2,3] or len(r['analysis']['step_records'])!=12:raise ValueError('Run coverage differs')
        verify(Path(r['build']),sources)
        if r['source_sha256']['src/comm_mpi.f90']!=sources['src/comm_mpi.f90']:raise ValueError('Run recorded different candidate source')
        if digest(Path(r['binary']))!=r['binary_sha256']:raise ValueError('Built executable differs from tested snapshot')
        if digest(root/'climateJ4')!=r['binary_sha256']:raise ValueError('Binary changed')
        expected='bin-215-candidate/climate' if name!='checked-off' else 'bin-215-checked/climate'
        if not str(r['binary']).endswith(expected):raise ValueError('Wrong local build variant')
        if (r['environment'].get('WAVETRISK_VALIDATE_BLOCK_DYNAMICS')=='1')!=(name=='oracles-on'):raise ValueError('Oracle mode differs')
        equality=compare_windows(a.reference,root,'localJ4J6_checkpoint_0002.bin.zst',expected_new=[3])
        if not equality['passed'] or not workload_gate(ref['analysis'],r['analysis'])['passed']:raise ValueError('Numerical equality/workload failed')
        report=compare(root/'localJ4J6_checkpoint_0003.bin.zst',a.legacy/'localJ4J6_checkpoint_0003.bin.zst',expected_domains=40)
        if accept(report,exact=False):raise ValueError('Legacy tolerance failed')
        result['runs'][name]=dict(path=str(root),binary_sha256=r['binary_sha256'],exact_comparison=equality,legacy_comparison=report,
                                 reloads=r['reloads'],steps=12,oracles=(name=='oracles-on'),source_inventory=sources)
    build=Path(json.loads((a.optimized/'results.json').read_text())['build'])
    checked_flags=(build/'stage215-audit.mk').read_text()
    if checked_flags!='FFLAGS += -DSTAGE215_AUDIT\n':raise ValueError('Route audit build flag differs')
    result['build_flags']={n:(build/n).read_text() for n in ('stage215-flags.mk','stage215-audit.mk')}
    result['build_logs']={n:digest(build.parent/n) for n in ('build-candidate.log','build-checked.log','rebuild-candidate.log')}
    result['halo_contract']=json.loads((a.halo_contract/'results.json').read_text())
    if result['halo_contract']['status']!='passed':raise ValueError('Halo contract missing')
    result['physics_objects']={p.name:digest(p) for p in (build/'src/physics/simple_physics/phyparam/obj').glob('*.o')}
    expected={p.name:digest(p) for p in (Path(ref['build'])/'src/physics/simple_physics/phyparam/obj').glob('*.o')}
    if len(expected)!=23 or result['physics_objects']!=expected:raise ValueError('Local physics objects differ')
    a.out.write_text(json.dumps(result,indent=2)+'\n');print('PASS: optimized off/on and route-audited checked off; exact restarts and workload; fixed legacy tolerance.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('parent-kit','reference','legacy','optimized','checked-off','oracles-on','halo-contract','out'):p.add_argument('--'+name,type=Path,required=True)
    validate(p.parse_args())
