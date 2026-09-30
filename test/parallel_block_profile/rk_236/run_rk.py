"""Source-bound RK cycle validation; no performance or scaling claims."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from primitives_234.run_primitives import run as run_base
from capability_223 import digest,save

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    run_base(binary,fixture,out,ranks,oracles,full,reference,False)
    report=json.loads((out/'primitives-results.json').read_text())
    text=(out/'run.log').read_text()
    closure=[list(map(int,m)) for m in re.findall(r'Retained RK closure rank=(\d+) stage=(\d+) aliases/owned operands=\s*(\d+)\s+(\d+)',text)]
    exact=[list(map(int,m)) for m in re.findall(r'Retained RK exact aliases rank=(\d+) stage=(\d+) count=(\d+)',text)]
    q1=[int(n) for n in re.findall(r'Retained RK q1 absent rank=(\d+)',text)]
    returns=[int(n) for n in re.findall(r'Retained RK native return arrays absent rank=(\d+)',text)]
    if control:
        assert not closure and not exact
    else:
        for rank in range(ranks):
            for stage in range(1,5):
                assert sum(row[:2]==[rank,stage] for row in closure)==report['steps'],(rank,stage,closure)
                if oracles:assert sum(row[:2]==[rank,stage] for row in exact)==report['steps']
        if not oracles:
            for rank in range(ranks):
                assert q1.count(rank)==report['steps']
                assert returns.count(rank)==4*report['steps']
        if ranks>40:assert any(rank==40 and owned>0 for rank,stage,aliases,owned in closure)
    report.pop('complete_tendency_migration',None)
    report.update(stage=236,control=control,scope='Retained RK closure and joined forward/inverse boundaries',rk_closure=closure,exact_alias_checks=exact,
                  q1_absence=q1,return_storage_absence=returns,performance_claim=None)
    save(out/'rk-results.json',report)
    print('PASS RK cycle '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('oracles','full','control'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
