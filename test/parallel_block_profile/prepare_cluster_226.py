#!/usr/bin/env python3
"""Freeze Stage226 production integration against the reviewed Stage225 source."""
import argparse,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from capability_223 import digest,save,verify_kit,verify_sources
from capability_226 import schedule
from production_226 import candidate,audit
from prepare_cluster_223 import inventory
BASE=Path(__file__).parent
HELPERS=(
 'capability_223.py','capability_226.py','production_226.py','placement_226.py','placement_223.py','placement_200.py',
 'build_capability_226.py','cluster_226.py','build_226.sbatch','cluster_226.sbatch',
 'bootstrap_226.sh','launch_226.sh','collect_226.sh','collect_cluster_logs_226.py',
 'rank_launch_190.py','rank_launch_223.py','rank_metrics_190.py','profile_protocol.py',
 'experiment.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
 'run_local_220.py','local_j4.py','prepare_physics_dependencies.sh',
 'test_launch_226.py','test_collect_226.py','test_campaign_226.py','test_capability_223.py','test_restart_checkpoints.py',
 'collect_cluster_logs_223.py','verify_workflow_226.py','test_task_capacity_226.py','run_local_226.py',
 'prepare_cluster_226.py','prepare_cluster_223.py','candidate_capability_223.py','stage226-production.patch')

def prepare(repo,previous,fixture,local_tree,out):
    for p in (repo,previous,fixture,local_tree):assert p.is_dir(),p
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external output required')
    verify_kit(previous);old=json.loads((previous/'stage225-sources.json').read_text())
    local=json.loads((repo/'docs/parallel-blocks-stage226-local-results.json').read_text())
    assert local['status']=='passed' and local['guards']['status']=='passed' and local['production']['status']=='passed'
    for name,h in local['production']['production_sources'].items():assert digest(repo/name)==h,name
    main=subprocess.check_output(['git','rev-parse','main'],cwd=repo,text=True).strip();assert main==old['main_commit']
    out.mkdir()
    for name in HELPERS:shutil.copy2(BASE/name,out/name)
    shutil.copy2(previous/'stage225-source.tar.gz',out/'stage226-source.tar.gz')
    with tempfile.TemporaryDirectory(prefix='stage226-package-') as tmp:
        root=Path(tmp)/'source';root.mkdir()
        with tarfile.open(out/'stage226-source.tar.gz') as t:t.extractall(root,filter='data')
        reference=inventory(root);assert reference==old['reference_sources']
        candidate(root);production=inventory(root)
        changed=[n for n in production if production[n]!=reference[n]]
        assert set(changed)==set(local['production']['production_sources'])
        for n in changed:assert digest(repo/n)==production[n],n
        verify_sources(local_tree,production)
        for reverse in (True,False):
            subprocess.run(['git','apply',*(['--reverse'] if reverse else []),str(out/'stage226-production.patch')],cwd=root,check=True)
            verify_sources(root,reference if reverse else production)
        audit(root);checked=inventory(root)
        assert all('WAVETRISK_STAGE223_ALLOW' not in (root/n).read_text() for n in ('src/arch.f90','src/main.f90'))
    m=dict(stage=226,reference_stage=225,source_commit=old['source_commit'],main_commit=main,
           reference_sources=reference,candidate_sources=production,checked_sources=checked,
           source_archive_sha256=digest(out/'stage226-source.tar.gz'),production_patch_sha256=digest(out/'stage226-production.patch'),
           previous_source_manifest_sha256=digest(previous/'stage225-sources.json'),
           scope='Production has no campaign bypass or lifecycle output. Separate checked observer records lifecycle; frozen baseline retains Domain limit.')
    save(out/'stage226-sources.json',m)
    (out/'fixture-J4').mkdir();(out/'grids').mkdir()
    for n in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copy2(fixture/n,out/'fixture-J4'/n)
    for n in ('HR95JT_WT_003','HR95JT_WT_004'):shutil.copy2(fixture/'grids'/n,out/'grids'/n)
    for n,label in [('parallel-blocks-cluster-226.md','README.md'),('parallel-blocks-stage226.md','EXPERIMENT.md'),('parallel-blocks-stage226-local-results.json','LOCAL_VALIDATION.json')]:shutil.copy2(repo/'docs'/n,out/label)
    save(out/'helpers-sha256.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    verify_kit(out);archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='stage226-reextract-') as tmp:
        with tarfile.open(archive) as t:t.extractall(tmp,filter='data')
        kit=Path(tmp)/out.name;verify_kit(kit)
        root=Path(tmp)/'source';root.mkdir()
        with tarfile.open(kit/'stage226-source.tar.gz') as t:t.extractall(root,filter='data')
        verify_sources(root,reference);candidate(root);verify_sources(root,production);audit(root);verify_sources(root,checked)
    result=dict(stage=226,kit=str(out),archive=str(archive),sha256=digest(archive),bytes=archive.stat().st_size,
                archive_reextracted_and_all_sources_reproduced=True,main_commit=main,gates=len(list(schedule())),cluster_status='prepared; not submitted')
    save(repo/'docs/parallel-blocks-stage226-kit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ('repo','previous','fixture','local-tree','out'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();prepare(a.repo.resolve(),a.previous.resolve(),a.fixture.resolve(),a.local_tree.resolve(),a.out.resolve())
