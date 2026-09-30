"""Require exact soil entry/output oracles and the accepted atmospheric lifecycle."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from physics_242.run_physics import run as run_previous
from capability_223 import save,digest
from analyze import analyze

def review(out,ranks,soil_check,control,reference):
    r=json.loads((out/'native-state-results.json').read_text())
    log=(out/'run.log').read_text()
    rows=[tuple(map(int,m)) for m in re.findall(r'Native soil transform rank=(\d+) layers=(\d+)',log)]
    entry=[int(m) for m in re.findall(r'Native soil entry exact rank=(\d+)',log)]
    output=[int(m) for m in re.findall(r'Native soil output exact rank=(\d+)',log)]
    if control:assert not rows and not entry and not output
    else:
        for rank in range(ranks):assert sum(x==rank for x,n in rows)==r['steps']
        assert rows and all(n==11 for x,n in rows)
        for rank in range(ranks):
            assert entry.count(rank)==(r['steps'] if soil_check else 0)
            assert output.count(rank)==(r['steps'] if soil_check else 0)
    if reference:
        actual=analyze(out/'run.log')['step_records']
        expected=analyze(reference/'run.log')['step_records']
        assert [x['state'] for x in actual]==[x['state'] for x in expected[:len(actual)]]
    r.update(stage=243,soil=rows,soil_check=soil_check,soil_entry_checks=entry,log_sha256=digest(out/'run.log'))
    save(out/'soil-results.json',r)
    print('PASS native soil '+str(out),flush=True)

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False,soil_check=False):
    run_previous(binary,fixture,out,ranks,oracles,full,reference,control)
    review(out,ranks,soil_check,control,reference)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control','soil-check'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
