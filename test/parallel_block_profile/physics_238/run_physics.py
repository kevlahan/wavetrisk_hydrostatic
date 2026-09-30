"""Validate explicit velocity read versions, preserving the legacy column order."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from physics_237.run_physics import run as run_previous
from capability_223 import save
from analyze import analyze

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    run_previous(binary,fixture,out,ranks,oracles,full,reference,False)
    r=json.loads((out/'physics-results.json').read_text())
    rows=[list(map(int,m)) for m in re.findall(r'Physics velocity plan rank=(\d+) columns/entry/updated/checks=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',(out/'run.log').read_text())]
    if control:assert not rows
    else:
        for rank in range(ranks):assert sum(row[0]==rank for row in rows)==r['steps']
        assert sum(updated for rank,n,entry,updated,checks in rows)>0
        assert all(entry+updated==6*n for rank,n,entry,updated,checks in rows)
        assert all(checks==(30*n if oracles else 0) for rank,n,entry,updated,checks in rows)
        if ranks>40:assert sum(row==[40,0,0,0,0] for row in rows)==r['steps']
    if not full:
        assert reference is not None
        assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[analyze(reference/'run.log')['step_records'][0]['state']]
    r.update(stage=238,velocity_plan=rows,velocity_plan_control=control,scope='Explicit column velocity operands and versions; serial Domain-owned schedule preserved')
    save(out/'velocity-results.json',r)
    print('PASS physics velocity '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
