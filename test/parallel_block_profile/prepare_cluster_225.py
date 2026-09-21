#!/usr/bin/env python3
"""Freeze the Stage224 reference plus the validated empty-Domain oracle fix."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from capability_223 import digest,save,verify_kit,verify_sources
from candidate_capability_223 import candidate,NAMES
from prepare_cluster_223 import inventory
from test_checkpoint_prefix_224 import instrument

BASE=Path(__file__).parent
HELPERS=(
 'capability_223.py','candidate_capability_223.py','placement_223.py','placement_200.py',
 'build_capability_225.py','cluster_225.py','build_225.sbatch','cluster_225.sbatch',
 'bootstrap_225.sh','launch_225.sh','collect_225.sh','collect_cluster_logs_225.py',
 'rank_launch_190.py','rank_launch_223.py','rank_metrics_190.py','profile_protocol.py',
 'experiment.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
 'run_local_220.py','local_j4.py','prepare_physics_dependencies.sh',
 'test_launch_225.py','test_collect_225.py','test_capability_223.py','test_restart_checkpoints.py',
 'collect_cluster_logs_223.py','test_checkpoint_prefix_224.py','verify_workflow_225.py',
 'prepare_cluster_225.py','prepare_cluster_223.py','stage225-oracle-extent.patch')


def prepare(repo,previous,fixture,out):
    repo=repo.resolve(strict=True);previous=previous.resolve(strict=True);fixture=fixture.resolve(strict=True);out=out.resolve()
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external output required')
    verify_kit(previous)
    old=json.loads((previous/'stage224-sources.json').read_text())
    local=json.loads((repo/'docs/parallel-blocks-stage225-local-results.json').read_text())
    assert local['status']=='passed' and local['regression']['status']=='passed'
    assert local['regression']['oracles'] and local['predicate_test_cases']==10
    for name,h in local['build_identity']['production_source_sha256'].items():
        assert digest(repo/name)==h,name
    main=subprocess.check_output(['git','rev-parse','main'],cwd=repo,text=True).strip()
    assert main==old['main_commit']==local['main_commit']
    assert set(local['regression']['cases'])=={'single-source0','single-source1'}
    for case in local['regression']['cases'].values():
        assert case['numerical']['passed'] and case['numerical']['exact_observed']
        assert case['events']['steps_per_rank']==12 and case['events']['reloads_per_rank']==2
        assert case['coverage']['remaps']==2
        assert len(case['oracle_switches'])==3 and set(case['oracle_switches'].values())=={'1'}
        assert all(p['all_ranks_exact'] for p in case['prefix'].values())
    out.mkdir()
    for name in HELPERS:shutil.copy2(BASE/name,out/name)
    with tempfile.TemporaryDirectory(prefix='stage225-package-') as tmp:
        scratch=Path(tmp);root=scratch/'source';root.mkdir()
        with tarfile.open(previous/'stage224-source.tar.gz') as t:t.extractall(root,filter='data')
        verify_sources(root,old['reference_sources'])
        shutil.copy2(repo/'src/parallel_block_mpi.f90',root/'src/parallel_block_mpi.f90')
        reference=inventory(root)
        assert set(reference)==set(old['reference_sources'])
        assert [n for n in reference if reference[n]!=old['reference_sources'][n]]==['src/parallel_block_mpi.f90']
        with tarfile.open(out/'stage225-source.tar.gz','w:gz') as t:
            for p in sorted(root.iterdir()):t.add(p,arcname=p.name)
        before={n:(root/n).read_text() for n in NAMES};candidate(root);trial=inventory(root)
        patch=''.join(''.join(difflib.unified_diff(before[n].splitlines(True),(root/n).read_text().splitlines(True),fromfile='a/'+n,tofile='b/'+n)) for n in NAMES)
        (out/'stage225-candidate.patch').write_text(patch)
        for reverse in (True,False):
            subprocess.run(['git','apply',*(['--reverse'] if reverse else []),str(out/'stage225-candidate.patch')],cwd=root,check=True)
            verify_sources(root,reference if reverse else trial)
        # Reapply only the test instrumentation and require the exact locally
        # validated Fortran sources, including the nonzero-source regression.
        instrument(root,repo/'src/checkpoint.f90')
        expected=local['build_identity']['sources']
        checked=0
        for name,h in expected.items():
            if name in trial:
                assert digest((root/name).resolve(strict=True))==h,name
                checked+=1
        assert checked==136
        # The local build also records generated physics files and the two
        # selected test-case links; these do not exist in a clean source tree.
        assert all(n.startswith('src/physics/simple_physics/phyparam/tmp/') or
                   n in {'src/test.f90','src/test_case_module.f90'}
                   for n in expected if n not in trial)
    m=dict(stage=225,reference_stage='224 plus allocation-consistent empty-Domain oracle extent correction',
           source_commit=old['source_commit'],main_commit=main,reference_sources=reference,candidate_sources=trial,
           source_archive_sha256=digest(out/'stage225-source.tar.gz'),candidate_patch_sha256=digest(out/'stage225-candidate.patch'),
           previous_source_manifest_sha256=digest(previous/'stage224-sources.json'),
           locally_validated_fortran_files_matched=checked,
           scope='Fresh corrected reference and opt-in capability candidate; no production guard change or branch mutation.')
    save(out/'stage225-sources.json',m)
    (out/'fixture-J4').mkdir();(out/'grids').mkdir()
    for name in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copy2(fixture/name,out/'fixture-J4'/name)
    for name in ('HR95JT_WT_003','HR95JT_WT_004'):shutil.copy2(fixture/'grids'/name,out/'grids'/name)
    for name,label in [('parallel-blocks-cluster-225.md','README.md'),('parallel-blocks-stage225-campaign.md','EXPERIMENT.md'),
                       ('parallel-blocks-stage225-local-results.json','LOCAL_VALIDATION.json'),
                       ('parallel-blocks-stage225-oracle-fix.md','CORRECTION.md')]:
        shutil.copy2(repo/'docs'/name,out/label)
    save(out/'helpers-sha256.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    verify_kit(out)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='stage225-reextract-') as tmp:
        with tarfile.open(archive) as t:t.extractall(tmp,filter='data')
        extracted=Path(tmp)/out.name;verify_kit(extracted)
        root=Path(tmp)/'source';root.mkdir()
        with tarfile.open(extracted/'stage225-source.tar.gz') as t:t.extractall(root,filter='data')
        verify_sources(root,reference);candidate(root);verify_sources(root,trial)
    result=dict(stage=225,kit=str(out),archive=str(archive),sha256=digest(archive),bytes=archive.stat().st_size,
                archive_reextracted_and_candidate_reproduced=True,main_commit=main,
                locally_validated_fortran_files_matched=checked,cluster_status='prepared; not submitted',
                scope='Eleven sequential correctness gates on nk, two 128-core sockets; no speedup claim.')
    save(repo/'docs/parallel-blocks-stage225-kit.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','previous','fixture','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.previous,a.fixture,a.out)
