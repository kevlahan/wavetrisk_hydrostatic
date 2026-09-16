#!/usr/bin/env python3
"""Freeze Stage200 numerical sources with whole-phase AMD PMU observers."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile,statistics as st
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
HELPERS=('build_restriction_201.py','build_201.sbatch','submit_restriction_201.py','cluster_201.py','cluster_201.sbatch',
 'analyze_hardware_201.py','rank_launch_201.py','preflight_201.py','collect_cluster_logs_201.py','workflow_clock_201.c')

def cohort(result):
    """Freeze common ranks covering each side's low/median/high local dynamics work on each socket."""
    mapping=result['cases']['n83-s2']['mapping'];selected=set();selections=[]
    for side in ('candidate','legacy'):
        cells=[c for c in result['runs']['n83-s2-'+side+'-probe']['scaling']['cells'] if c['phase'] in ('dynamics WT','dynamics outside WT')]
        scores={rank:sum(c['rank_rows'][rank]['outside_mpi_seconds'] for c in cells)/len({c['step'] for c in cells}) for rank in range(83)}
        for socket in (0,1):
            ranks=sorted([m['rank'] for m in mapping if m['socket']==socket],key=lambda rank:(scores[rank],rank))
            for label,index in (('low',0),('median',len(ranks)//2),('high',len(ranks)-1)):
                rank=ranks[index];selected.add(rank);selections.append(dict(side=side,socket=socket,stratum=label,rank=rank,local_dynamics_seconds=scores[rank]))
    return {'ranks':sorted(selected),'selections':selections,'source':'Stage200 job54689 n83-s2 probe; deterministic union shared by both implementations','scope':'Stratified diagnostic cohort, not an unbiased node-wide average.'}

def prepare(repo,out,previous,review):
    repo=repo.resolve();out=out.resolve();previous=previous.resolve();base=Path(__file__).parent
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit/archive required')
    old=json.loads((previous/'stage200-sources.json').read_text())
    for name,h in old['candidate_sources'].items():
        if file_digest(repo/name)!=h:raise ValueError('Production numerics differ: '+name)
    shutil.copytree(previous,out);manifest=dict(old)
    for name in HELPERS:shutil.copy2(base/name,out/name)
    shutil.copy2(previous/'stage200-candidate.patch',out/'stage201-candidate.patch')
    for side in ('block','legacy'):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            if side=='block':
                with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(out/'stage201-candidate.patch')],cwd=root,check=True)
                original=old['candidate_sources'];oldpatch='stage200-observer.patch';key='profiled_sources';target='stage201-observer.patch';hashkey='observer_patch_sha256'
            else:
                with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
                original=old['legacy_io_sources'];oldpatch='legacy-profile-200.patch';key='legacy_profile_sources';target='legacy-profile-201.patch';hashkey='legacy_profile_patch_sha256'
            before={n:(root/n).read_bytes() for n in original if not (root/n).is_symlink()}
            if file_digest(previous/oldpatch)!=old[hashkey]:raise ValueError('Frozen observer changed')
            subprocess.run(['git','apply','--whitespace=nowarn',str(previous/oldpatch)],cwd=root,check=True)
            for n,h in old[key].items():
                if file_digest(root/n)!=h:raise ValueError('Frozen observer source differs: '+n)
            (root/'src/workflow_profile_195.f90').write_bytes((base/'workflow_profile_201.f90').read_bytes())
            (root/'src/workflow_clock_199.c').write_bytes((base/'workflow_clock_201.c').read_bytes())
            manifest[key]={n:file_digest(root/n) for n in old[key]}
            (out/target).write_text(patch(before,root,[n for n,h in manifest[key].items() if h!=original.get(n)]))
            manifest[hashkey]=file_digest(out/target)
    evidence=json.loads(review.read_text())
    if evidence.get('stage')!=200 or evidence.get('status')!='passed' or not all(g['passed'] for g in evidence['gates'].values()):raise ValueError('Passed Stage200 evidence required')
    selected=cohort(evidence);selected['source_results_sha256']=file_digest(review)
    (out/'sample-cohort-201.json').write_text(json.dumps(selected,indent=2)+'\n')
    manifest.update(stage=201,observer_sha256=file_digest(base/'workflow_profile_201.f90'),diagnosis_only=True)
    (out/'stage201-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for source,target in [('parallel-blocks-cluster-201.md','README.md'),('parallel-blocks-stage201.md','EXPERIMENT.md')]:shutil.copy2(repo/'docs'/source,out/target)
    helpers=json.loads((previous/'helpers-sha256.json').read_text())
    for name in (*HELPERS,'sample-cohort-201.json'):helpers[name]=file_digest(out/name)
    (out/'helpers-sha256.json').write_text(json.dumps(helpers,indent=2)+'\n')
    validation=repo/'docs/parallel-blocks-stage201-local-results.json'
    if validation.exists():shutil.copy2(validation,out/'LOCAL_VALIDATION.json')
    else:(out/'LOCAL_VALIDATION.json').unlink(missing_ok=True)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    report={'stage':201,'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],'numerical_candidate_patch_sha256':manifest['candidate_patch_sha256'],'observer_patch_sha256':manifest['observer_patch_sha256'],'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'sample_ranks':selected['ranks'],'scope':'Five fresh builds; capability preflight and 20 solver runs at fixed83 two-socket placement. No numerical change.'}
    print(json.dumps(report,indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('repo','out','stage200-kit','stage200-results'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.stage200_kit,a.stage200_results)
