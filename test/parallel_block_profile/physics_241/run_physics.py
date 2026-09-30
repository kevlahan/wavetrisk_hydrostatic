"""Qualify owner-resident dependency waves against the accepted serial path."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from physics_240.run_physics import run as run_previous
from capability_223 import save

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    run_previous(binary,fixture,out,ranks,oracles,full,reference,False)
    r=json.loads((out/'batch-results.json').read_text())
    rows=[list(map(int,m)) for m in re.findall(r'Physics native input rank=(\d+) canonical requests=(\d+) alias requests=(\d+)',(out/'run.log').read_text())]
    if control:assert not rows
    else:
        for rank in range(ranks):assert sum(row[0]==rank for row in rows)==r['steps']
        assert sum(c for rank,c,a in rows)>0 and sum(a for rank,c,a in rows)>0
    r.update(stage=241,native_inputs=rows,native_input_control=control,scope='Canonical owner inputs and saved atmospheric aliases; alias mean/soil adapter and publication remain')
    save(out/'native-input-results.json',r);print('PASS native physics inputs '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
