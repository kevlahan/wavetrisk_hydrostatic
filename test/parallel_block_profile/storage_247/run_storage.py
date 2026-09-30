"""Validate exact trajectories and the absence of persistent legacy numerical views."""
import argparse,json,re,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from lifecycle_245.run_lifecycle import run as previous
from capability_223 import digest,save

def review(out,ranks,views):
    r=json.loads((out/'lifecycle-results.json').read_text());log=(out/'run.log').read_text()
    absence=list(map(int,re.findall(r'Domain prognostic storage absent rank=(\d+)',log)))
    output=list(map(int,re.findall(r'Native output alias adapter rank=(\d+)',log)))
    scratch=list(map(int,re.findall(r'Legacy scalar staging and RK snapshots absent rank=(\d+)',log)))
    if views:assert not absence and not scratch
    else:
        for rank in range(ranks):
            assert absence.count(rank)==r['steps']+output.count(rank),(rank,absence.count(rank))
            assert scratch.count(rank)==output.count(rank)>0
    save(out/'storage-results.json',dict(stage=247,status='passed',ranks=ranks,steps=r['steps'],legacy_views=views,
         absence_checks=len(absence),bounded_output_checks=len(scratch),
         binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log'),
         lifecycle_results_sha256=digest(out/'lifecycle-results.json')))
    print('PASS persistent storage qualification '+str(out),flush=True)

def run(binary,fixture,out,ranks=4,full=False,reference=None,control=False,checks=False,oracles=False,views=False):
    previous(binary,fixture,out,ranks,full,reference,control,checks,oracles)
    review(out,ranks,views)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    for n in ('full','control','checks','oracles','views'):p.add_argument('--'+n,action='store_true')
    p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
