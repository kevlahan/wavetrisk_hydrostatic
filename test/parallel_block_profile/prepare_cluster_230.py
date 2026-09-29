#!/usr/bin/env python3
"""Seal validated Stage230 source and a self-contained baseline, without Git mutations."""
import argparse,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from capability_223 import digest,save,verify_kit,verify_sources
BASE=Path(__file__).resolve().parent
HELPERS=('build_capability_230.py','cluster_230.py','protocol_230.py','protocol_229.py','validation_230.py',
 'placement_230.py','placement_200.py','bootstrap_230.sh','launch_230.sh','collect_230.sh',
 'build_230.sbatch','cluster_230.sbatch','collect_cluster_logs_230.py',
 'test_campaign_230.py','test_launch_230.py','test_collect_230.py','verify_workflow_230.py','verify_build_230.py',
 'prepare_cluster_230.py','run_local_230.py',
 'capability_223.py','rank_launch_190.py','rank_launch_223.py','rank_metrics_190.py','profile_protocol.py',
 'experiment.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py','run_local_220.py','local_j4.py',
 'prepare_physics_dependencies.sh')


def inventory(root):
    paths=[root/'Makefile',*(root/'src').rglob('*'),*(root/'test/climate').rglob('*')]
    return {str(p.relative_to(root)):digest(p) for p in sorted(paths) if p.is_file() or p.is_symlink()}


def prepare(repo,local,out):
    if out.exists() or out.with_suffix('.tar.gz').exists():raise ValueError('Fresh kit destination required')
    gate=json.loads((repo/'docs/parallel-blocks-stage230-local-results.json').read_text())
    proof=json.loads((local/'build-integrity.json').read_text())
    if gate['status']!='passed' or proof['status']!='passed':raise ValueError('Missing local proof')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    main=subprocess.check_output(['git','rev-parse','main'],cwd=repo,text=True).strip()
    if commit!=gate['reference_commit'] or commit!='5c457a5dc869ce18938dd24db5727053b8e8dc2a':raise ValueError('Wrong baseline')
    verify_sources(repo,gate['source_sha256'])
    out.mkdir(parents=True)
    for name in HELPERS:shutil.copy2(BASE/name,out/name)
    with tempfile.TemporaryDirectory(prefix='stage230-package-') as tmp:
        scratch=Path(tmp);ref=scratch/'reference';ref.mkdir()
        subprocess.run(['git','archive','--format=tar','--output='+str(scratch/'head.tar'),commit],cwd=repo,check=True)
        with tarfile.open(scratch/'head.tar') as t:t.extractall(ref,filter='data')
        reference=inventory(ref)
        with tarfile.open(out/'stage230-source.tar.gz','w:gz') as t:
            for p in sorted(ref.iterdir()):t.add(p,arcname=p.name)
        patch=subprocess.check_output(['git','diff','--binary',commit,'--',*sorted(gate['source_sha256'])],cwd=repo)
        (out/'stage230-candidate.patch').write_bytes(patch)
        subprocess.run(['git','apply',str(out/'stage230-candidate.patch')],cwd=ref,check=True)
        candidate=inventory(ref)
        changed={n for n in set(reference)|set(candidate) if reference.get(n)!=candidate.get(n)}
        if changed!=set(gate['source_sha256']):raise ValueError('Unexpected patch scope')
        for name,h in proof['candidate_sources'].items():
            if digest((ref/name).resolve(strict=True))!=h:raise ValueError('Candidate differs from validated source: '+name)
        verify_sources(ref,gate['source_sha256'])
        subprocess.run(['git','apply','--reverse',str(out/'stage230-candidate.patch')],cwd=ref,check=True)
        verify_sources(ref,reference)
    m=dict(stage=230,reference_stage=229,source_commit=commit,main_commit=main,
           reference_sources=reference,candidate_sources=candidate,
           source_archive_sha256=digest(out/'stage230-source.tar.gz'),candidate_patch_sha256=digest(out/'stage230-candidate.patch'),
           local_proof_sha256=digest(repo/'docs/parallel-blocks-stage230-local-results.json'),
           scope='Seven validated source edits remove Domain tendency storage for the qualified production path. Both versions STAGE202/203 on. No main changes; no legacy/observer build.')
    save(out/'stage230-sources.json',m)
    (out/'grids').mkdir();shutil.copy2(local/'assets/grids/HR95JT_WT_004',out/'grids/HR95JT_WT_004')
    for name,label in [('parallel-blocks-cluster-230.md','README.md'),('parallel-blocks-stage230.md','EXPERIMENT.md'),
                       ('parallel-blocks-stage230-local-results.json','LOCAL_VALIDATION.json')]:
        shutil.copy2(repo/'docs'/name,out/label)
    save(out/'helpers-sha256.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    verify_kit(out);archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as ar:ar.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='stage230-unpack-') as tmp:
        with tarfile.open(archive) as ar:ar.extractall(tmp,filter='data')
        kit=Path(tmp)/out.name;verify_kit(kit)
        root=Path(tmp)/'root';ref=root/'reference';ref.mkdir(parents=True)
        with tarfile.open(kit/'stage230-source.tar.gz') as ar:ar.extractall(ref,filter='data')
        from build_capability_230 import materialize,source_key
        for kind in ('baseline','candidate','checked'):verify_sources(materialize(kit,root,kind),m[source_key(kind)])
    result=dict(stage=230,kit=str(out),archive=str(archive),sha256=digest(archive),bytes=archive.stat().st_size,
                reference_commit=commit,main_commit=main,source_reproduction_verified=True,cluster_status='prepared; not submitted')
    save(repo/'docs/parallel-blocks-stage230-kit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('repo','local','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.repo.resolve(),a.local.resolve(),a.out.resolve())
