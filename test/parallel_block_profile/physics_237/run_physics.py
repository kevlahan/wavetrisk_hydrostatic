"""Validate column-local Simple-physics hydrostatics with retained RK enabled."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from rk_236.run_rk import run as run_rk
from capability_223 import save

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    run_rk(binary,fixture,out,ranks,oracles,full,reference,False)
    r=json.loads((out/'rk-results.json').read_text());log=(out/'run.log').read_text()
    rows=[list(map(int,m)) for m in re.findall(r'Local physics columns rank=(\d+) columns/layer checks=\s*(\d+)\s+(\d+)',log)]
    if control:assert not rows
    else:
        for rank in range(ranks):assert sum(row[0]==rank for row in rows)==r['steps']
        assert sum(row[1] for row in rows)>0
        if oracles:assert all(checks==30*columns for rank,columns,checks in rows)
        else:assert all(checks==0 for rank,columns,checks in rows)
        if ranks>40:assert any(rank==40 and columns==0 for rank,columns,checks in rows)
    r.update(stage=237,control=control,physics_columns=rows,scope='Simple physics hydrostatic scratch only; physics remains on geometry owners',complete_physics_owner_migration=False)
    save(out/'physics-results.json',r);print('PASS physics '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
