"""Validate all-stage owner tendency execution; local timing is not speed evidence."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from primitives_234.run_primitives import run as primitive_run
from capability_223 import digest,save

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    primitive_run(binary,fixture,out,ranks,oracles,full,reference,control)
    result=json.loads((out/'primitives-results.json').read_text())
    log=(out/'run.log').read_text()
    initial=[list(map(int,m)) for m in re.findall(r'Initial tendency alias epoch rank=(\d+) registered aliases=(\d+)',log)]
    programs=[list(map(int,m)) for m in re.findall(r'Block owner program rank=(\d+) local/owned/depth=\s*(\d+)\s+(\d+)\s+(\d+)',log)]
    mass=[list(map(int,m)) for m in re.findall(r'Block mass owner rank=(\d+) checked=(\d+)',log)]
    source=[list(map(int,m)) for m in re.findall(r'Block velocity owner rank=(\d+) checked=(\d+)',log)]
    if control:
        assert not initial and not programs and not mass and not source
    else:
        assert {r for r,n in initial}==set(range(ranks)),initial
        assert {r for r,*rest in programs}==set(range(ranks)),programs
        assert {r for r,n in mass}==set(range(ranks)),mass
        assert all(sum(r==rank for r,n in initial)==result['steps'] for rank in range(ranks))
        # Every RK stage, including stage one, must execute the owner columns.
        if full:assert result['steps']==12
    result.update(stage=235,scope='All RK stages: hydrostatic columns, primitives, scalar diffusion/flux, mass restriction/divergence, rotational diffusion, velocity source and final gradients; existing block temperature restriction/divergence.',
        complete_tendency_kernel_migration=not control,
        complete_domain_storage_removal=False,
        input_boundary='Stage-one upstream aliases and reference fields remain explicit inputs; interiors read resident blocks.',
        output_boundary='Native mass/velocity RK closure and existing temperature scalar assembly remain explicit returns.',
        initial_epochs=initial,owner_programs=programs,mass_checks=mass,source_checks=source,
        performance_claim=None)
    result.pop('complete_tendency_migration',None)
    save(out/'tendency-results.json',result)
    print('PASS tendency '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4);p.add_argument('--oracles',action='store_true')
    p.add_argument('--full',action='store_true');p.add_argument('--reference',type=Path)
    p.add_argument('--control',action='store_true')
    run(**vars(p.parse_args()))
