"""Qualify numerical adaptation against the accepted complete trajectory."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from adaptation_244.run_adaptation import run as previous
from capability_223 import digest,save

def run(binary,fixture,out,ranks=4,full=False,reference=None,control=False,checks=False,oracles=False):
    previous(binary,fixture,out,ranks,full,reference,control,checks,oracles)
    r=json.loads((out/'adaptation-results.json').read_text());log=(out/'run.log').read_text()
    if not control:
        bootstrap=list(map(int,re.findall(r'Native lifecycle bootstrap rank=(\d+)',log)))
        output=list(map(int,re.findall(r'Native output alias adapter rank=(\d+)',log)))
        for rank in range(ranks):
            assert bootstrap.count(rank)>=(2 if full else 1)
            assert output.count(rank)>=1
    save(out/'lifecycle-results.json',dict(stage=245,status='passed',ranks=ranks,steps=r['steps'],checks=checks,control=control,oracles=oracles,
        binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log'),adaptation_results_sha256=digest(out/'adaptation-results.json')))
    print('PASS native lifecycle '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('full','control','checks','oracles'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
