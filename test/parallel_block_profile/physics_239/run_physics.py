"""Validate ordered package execution on final block owners."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from physics_238.run_physics import run as run_previous
from capability_223 import save
from analyze import analyze

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    run_previous(binary,fixture,out,ranks,oracles,full,reference,False)
    r=json.loads((out/'velocity-results.json').read_text())
    rows=[list(map(int,m)) for m in re.findall(r'Physics owner rank=(\d+) executed/remote/scaffold/checks=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',(out/'run.log').read_text())]
    if control:assert not rows
    else:
        for rank in range(ranks):assert sum(row[0]==rank for row in rows)==r['steps']
        assert sum(remote for rank,n,remote,scaffold,checks in rows)>0
        assert sum(n+scaffold for rank,n,remote,scaffold,checks in rows)==sum(n for rank,n,checks in r['physics_columns'])
        assert sum(row[4] for row in rows)==(sum(n for rank,n,checks in r['physics_columns']) if oracles else 0)
        if ranks>40:assert all(n>0 and remote==n and scaffold==0 for rank,n,remote,scaffold,checks in rows if rank==40)
    if not full:
        assert reference is not None
        assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[analyze(reference/'run.log')['step_records'][0]['state']]
    r.update(stage=239,owner_rows=rows,owner_control=control,scope='Ordered block-owner package execution; Domain operands and return compatibility remain')
    save(out/'owner-results.json',r)
    print('PASS physics owner '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
