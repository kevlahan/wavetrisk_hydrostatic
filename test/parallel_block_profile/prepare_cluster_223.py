#!/usr/bin/env python3
"""Package the validated Stage222 snapshot plus the bounded Stage223 prototype."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from capability_223 import digest, save, verify_kit, verify_sources
from candidate_capability_223 import candidate, NAMES

BASE=Path(__file__).parent
HELPERS=('capability_223.py','candidate_capability_223.py','placement_223.py','placement_200.py',
         'build_capability_223.py','cluster_223.py','prepare_cluster_223.py','build_223.sbatch','cluster_223.sbatch',
         'bootstrap_223.sh','launch_223.sh','collect_223.sh','collect_cluster_logs_223.py',
         'rank_launch_190.py','rank_launch_223.py','rank_launch_201.py','rank_metrics_190.py','profile_protocol.py',
         'experiment.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
         'run_local_220.py','local_j4.py','prepare_physics_dependencies.sh','run_local_223.py',
         'test_capability_223.py','test_collect_223.py','test_restart_checkpoints.py','verify_workflow_223.py')


def inventory(root):
    paths=[root/'Makefile',*(root/'src').rglob('*'),*(root/'test/climate').rglob('*')]
    return {str(p.relative_to(root)):digest(p) for p in sorted(paths) if p.is_file() or p.is_symlink()}


def prepare(repo, fixture, out):
    repo=repo.resolve(strict=True);fixture=fixture.resolve(strict=True);out=out.resolve()
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):
        raise ValueError('Fresh external kit required')
    gate=json.loads((repo/'docs/parallel-blocks-stage222-cleanup-results.json').read_text())
    local=json.loads((repo/'docs/parallel-blocks-stage223-local-results.json').read_text())
    if gate['status']!='passed' or local['status']!='passed':raise ValueError('Local correctness gates missing')
    verify_sources(repo,gate['workspace_sources'])
    out.mkdir()
    for name in HELPERS:shutil.copy2(BASE/name,out/name)
    with tempfile.TemporaryDirectory(prefix='stage223-package-') as tmp:
        scratch=Path(tmp);root=scratch/'source';root.mkdir()
        subprocess.run(['git','archive','--format=tar','--output='+str(scratch/'head.tar'),'HEAD'],cwd=repo,check=True)
        with tarfile.open(scratch/'head.tar') as t:t.extractall(root,filter='data')
        for name in gate['workspace_sources']:shutil.copy2(repo/name,root/name)
        reference=inventory(root)
        with tarfile.open(out/'stage223-source.tar.gz','w:gz') as t:
            for p in sorted(root.iterdir()):t.add(p,arcname=p.name)
        before={n:(root/n).read_text() for n in NAMES};candidate(root)
        trial=inventory(root)
        for name,h in local['sources'].items():
            # Local runtime records hash dereferenced source bytes; the kit's
            # tree inventory also preserves the identity of symlink targets.
            if name in trial and digest((root/name).resolve(strict=True))!=h:
                raise ValueError('Packaged candidate differs from local validation: '+name)
        if any(name not in trial for name in NAMES):raise ValueError('Incomplete prototype manifest')
        patch=''.join(''.join(difflib.unified_diff(before[n].splitlines(True),(root/n).read_text().splitlines(True),
                                                 fromfile='a/'+n,tofile='b/'+n)) for n in NAMES)
        (out/'stage223-candidate.patch').write_text(patch)
        for reverse in (True,False):
            subprocess.run(['git','apply',*(['--reverse'] if reverse else []),str(out/'stage223-candidate.patch')],cwd=root,check=True)
            verify_sources(root,reference if reverse else trial)
    m=dict(stage=223,reference_stage=222,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
           main_commit=subprocess.check_output(['git','rev-parse','main'],cwd=repo,text=True).strip(),
           reference_sources=reference,candidate_sources=trial,
           source_archive_sha256=digest(out/'stage223-source.tar.gz'),candidate_patch_sha256=digest(out/'stage223-candidate.patch'),
           scope='Uncommitted Stage222 cleanup frozen above HEAD; prototype changes only four files, with explicit opt-in. No branch mutation.')
    save(out/'stage223-sources.json',m)
    (out/'fixture-J4').mkdir();(out/'grids').mkdir()
    for name in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copy2(fixture/name,out/'fixture-J4'/name)
    for name in ('HR95JT_WT_003','HR95JT_WT_004'):shutil.copy2(fixture/'grids'/name,out/'grids'/name)
    for name,label in [('parallel-blocks-cluster-223.md','README.md'),('parallel-blocks-stage223.md','EXPERIMENT.md'),
                       ('parallel-blocks-stage223-local-results.json','LOCAL_VALIDATION.json'),
                       ('parallel-blocks-stage222-cleanup-results.json','STAGE222_VALIDATION.json')]:
        shutil.copy2(repo/'docs'/name,out/label)
    hashes={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
    save(out/'helpers-sha256.json',hashes);verify_kit(out)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='stage223-reextract-') as tmp:
        with tarfile.open(archive) as t:t.extractall(tmp,filter='data')
        extracted=Path(tmp)/out.name;verify_kit(extracted)
        source=Path(tmp)/'source';source.mkdir()
        with tarfile.open(extracted/'stage223-source.tar.gz') as t:t.extractall(source,filter='data')
        verify_sources(source,reference);candidate(source);verify_sources(source,trial)
    result=dict(stage=223,kit=str(out),archive=str(archive),sha256=digest(archive),bytes=archive.stat().st_size,
                archive_reextracted_and_candidate_reproduced=True,main_commit=m['main_commit'],
                scope='Eleven sequential correctness runs on nk, one exclusive 256-core node; no timing promotion.')
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','fixture','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.fixture,a.out)
