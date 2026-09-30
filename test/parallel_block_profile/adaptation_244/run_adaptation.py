"""Qualify numerical adaptation against the accepted complete trajectory."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from soil_243.run_soil import run as previous
from capability_223 import digest,save

def run(binary,fixture,out,ranks=4,full=False,reference=None,control=False,checks=False,oracles=False):
    previous(binary,fixture,out,ranks,oracles,full,reference,control,checks)
    r=json.loads((out/'soil-results.json').read_text());log=(out/'run.log').read_text()
    for marker,enabled in [('Native adaptation reconstruction rank=',not control),('Native adaptation exact rank=',checks)]:
        rows=list(map(int,re.findall(re.escape(marker)+r'(\d+)',log)))
        for rank in range(ranks):assert rows.count(rank)==(r['steps'] if enabled else 0)
    save(out/'adaptation-results.json',dict(stage=244,status='passed',ranks=ranks,steps=r['steps'],checks=checks,control=control,oracles=oracles,
        binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log'),soil_results_sha256=digest(out/'soil-results.json')))
    print('PASS native adaptation '+str(out),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('full','control','checks','oracles'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
