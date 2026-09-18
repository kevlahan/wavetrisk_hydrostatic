#!/usr/bin/env python3
"""Package Stage215 from the verified Stage210 kit's uninstrumented Stage202 sources."""
import argparse,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from build_restriction_208 import verify,digest,verify_helpers
from prepare_cluster_195 import patch
from candidate_boundary_215 import candidate
BASE=Path(__file__).parent
HELPERS=('candidate_boundary_215.py','protocol_215.py','cluster_215.py','build_restriction_215.py',
         'submit_restriction_215.py','prepare_cluster_215.py','build_215.sbatch','cluster_215.sbatch',
         'bootstrap_215.sh','launch_215.sh','collect_cluster_logs_215.py','test_protocol_215.py','test_bootstrap_215.py','test_collector_215.py','validate_local_215.py','test_boundary_215.f90','run_boundary_test_215.py','collect_215.sh','test_collect_215.py')

def prepare(previous,out,repo):
    previous=previous.resolve(strict=True);out=out.resolve();repo=repo.resolve(strict=True)
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit required')
    verify_helpers(previous);old=json.loads((previous/'stage210-sources.json').read_text())
    if old['stage']!=210 or old['baseline_commit']!='9bc3a19556dab79b3cba3944cf5df53517389da3':raise ValueError('Unexpected parent kit')
    shutil.copytree(previous,out,symlinks=True,ignore=shutil.ignore_patterns('__pycache__'))
    for name in HELPERS:shutil.copy2(BASE/name,out/name)
    shutil.copy2(previous/'stage210-baseline.patch',out/'stage215-baseline.patch')
    m={k:old[k] for k in ('baseline_commit','baseline_sources','reference_sources','makefile_sha256')}
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        with tarfile.open(previous/'stage196-source.tar.gz') as t:t.extractall(root,filter='data')
        verify(root,m['baseline_sources'])
        for n in ('test/climate/test_case_module.f90','test/climate/climate.f90'):m['baseline_sources'][n]=digest(root/n)
        subprocess.run(['git','apply',str(out/'stage215-baseline.patch')],cwd=root,check=True);verify(root,m['reference_sources'])
        for n in ('test/climate/test_case_module.f90','test/climate/climate.f90'):m['reference_sources'][n]=digest(root/n)
        names=['src/comm_mpi.f90']
        before={n:(root/n).read_bytes() for n in names};candidate(root)
        m['candidate_sources']=dict(m['reference_sources'])
        for n in names:m['candidate_sources'][n]=digest(root/n)
        (out/'stage215-candidate.patch').write_text(patch(before,root,names))
        subprocess.run(['git','apply','--reverse',str(out/'stage215-candidate.patch')],cwd=root,check=True);verify(root,m['reference_sources'])
        subprocess.run(['git','apply',str(out/'stage215-candidate.patch')],cwd=root,check=True);verify(root,m['candidate_sources'])
    m.update(stage=215,baseline_stage=202,diagnosis_only=False,
             reference_patch_sha256=digest(out/'stage215-baseline.patch'),candidate_patch_sha256=digest(out/'stage215-candidate.patch'),
             parent_manifest_sha256=digest(previous/'stage210-sources.json'))
    (out/'stage215-sources.json').write_text(json.dumps(m,indent=2)+'\n')
    for a,b in [('parallel-blocks-cluster-215.md','README.md'),('parallel-blocks-stage215.md','EXPERIMENT.md'),
                ('parallel-blocks-stage215-local-results.json','LOCAL_VALIDATION.json'),
                ('parallel-blocks-stage215-preparation-results.json','PREPARATION.json'),
                ('parallel-blocks-stage215-source-evidence.json','SOURCE_EVIDENCE.json')]:
        if (repo/'docs'/a).exists():shutil.copy2(repo/'docs'/a,out/b)
        else:(out/b).unlink(missing_ok=True)
    local=json.loads((out/'LOCAL_VALIDATION.json').read_text())
    if local['status']!='passed' or local['candidate_sources']!=m['candidate_sources']:raise ValueError('Local candidate gates missing/different')
    hashes=json.loads((previous/'helpers-sha256.json').read_text())
    for name in (*HELPERS,'stage215-sources.json','stage215-baseline.patch','stage215-candidate.patch','README.md','EXPERIMENT.md','LOCAL_VALIDATION.json','PREPARATION.json','SOURCE_EVIDENCE.json'):
        hashes[name]=digest(out/name)
    (out/'helpers-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n');verify_helpers(out)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(archive) as t:t.extractall(tmp,filter='data')
        verify_helpers(Path(tmp)/out.name)
    result=dict(stage=215,kit=str(out),archive=str(archive),sha256=digest(archive),scope='Three builds; adaptive correctness/16 AB/83 AB/83 ABBA funnel. No automatic promotion.')
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('previous','out','repo'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.previous,a.out,a.repo)
