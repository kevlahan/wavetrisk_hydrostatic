#!/usr/bin/env python3
"""Freeze Stage203 against the exact Stage202 production sources and legacy."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from instrument_workflow_195 import instrument
from build_restriction_203 import verify
HELPERS=('build_restriction_203.py','submit_restriction_203.py','build_203.sbatch','cluster_203.sbatch','cluster_203.py',
         'analyze_batching_203.py','preflight_203.py','collect_cluster_logs_203.py')
CHANGES=('Makefile','src/parallel_block_scalar_storage.f90','src/parallel_block_mpi.f90','src/parallel_block_inverse.f90','src/parallel_block.f90')

def observer(root,base):
    instrument(root,True)
    (root/'src/workflow_profile_195.f90').write_bytes((base/'workflow_profile_201.f90').read_bytes())
    (root/'src/workflow_clock_199.c').write_bytes((base/'workflow_clock_201.c').read_bytes())
    p=root/'Makefile';s=p.read_text();old='OBJ := $(addprefix $(BUILD_DIR)/,$(SRC:.f90=.o))'
    if s.count(old)!=1:raise ValueError('Missing object anchor')
    s=s.replace(old,old+' $(BUILD_DIR)/workflow_clock_199.o')
    # Append an explicit C rule; using mpicc avoids the system-cc include issue.
    s+='\nFFLAGS += -DSTAGE203_LEDGER\nMPICC ?= mpicc\n$(BUILD_DIR)/workflow_clock_199.o: $(SRC_DIR)/workflow_clock_199.c | dirs\n\t$(MPICC) -O2 -c $< -o $@\n'
    p.write_text(s)

def prepare(repo,out,previous):
    repo=repo.resolve();out=out.resolve();previous=previous.resolve();base=Path(__file__).parent
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit/archive required')
    old=json.loads((previous/'stage202-sources.json').read_text())
    # The only production differences allowed are the explicit prototype files.
    for n,h in old['candidate_sources'].items():
        if n not in CHANGES and file_digest(repo/n)!=h:raise ValueError('Unrelated source difference: '+n)
    for n,h in json.loads((previous/'helpers-sha256.json').read_text()).items():
        if file_digest(previous/n)!=h:raise ValueError('Changed previous helper: '+n)
    if file_digest(previous/'stage202-candidate.patch')!=old['candidate_patch_sha256']:raise ValueError('Previous numerical patch differs')
    shutil.copytree(previous,out);manifest=dict(old)
    for n in HELPERS:shutil.copy2(base/n,out/n)
    # Reconstruct the exact combined Stage202 patch against the pinned Stage196 source.
    (out/'stage203-baseline.patch').write_text('')
    manifest.update(stage=203,reference_sources=old['candidate_sources'],reference_patch_sha256=None,diagnosis_only=False)
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
        baseline_original={n:(root/n).read_bytes() for n in old['baseline_sources'] if (root/n).exists() and not (root/n).is_symlink()}
        for patchname in ('stage202-baseline.patch','stage202-candidate.patch'):
            subprocess.run(['git','apply',str(previous/patchname)],cwd=root,check=True)
        (out/'stage203-baseline.patch').write_text(patch(baseline_original,root,old['candidate_sources']))
        manifest['reference_patch_sha256']=file_digest(out/'stage203-baseline.patch')
        verify(root,manifest['reference_sources'])
        original={n:(root/n).read_bytes() for n in old['candidate_sources'] if not (root/n).is_symlink()}
        for n in CHANGES:shutil.copy2(repo/n,root/n)
        manifest['candidate_sources']={n:file_digest(root/n) for n in old['candidate_sources']}
        (out/'stage203-candidate.patch').write_text(patch(original,root,CHANGES))
        manifest['candidate_patch_sha256']=file_digest(out/'stage203-candidate.patch')
        for side,key,patchname,hashkey in [('candidate','profiled_sources','stage203-observer.patch','observer_patch_sha256'),
                                         ('baseline','baseline_profile_sources','stage203-baseline-observer.patch','baseline_observer_patch_sha256')]:
            sources=manifest['candidate_sources'] if side=='candidate' else manifest['reference_sources']
            before={n:(root/n).read_bytes() for n in sources if not (root/n).is_symlink()}
            observer(root,base)
            names=[*sources,'src/workflow_profile_195.f90','src/workflow_clock_199.c']
            manifest[key]={n:file_digest(root/n) for n in names}
            (out/patchname).write_text(patch(before,root,[n for n in names if manifest[key][n]!=sources.get(n)]))
            manifest[hashkey]=file_digest(out/patchname)
            subprocess.run(['git','apply','--reverse',str(out/patchname)],cwd=root,check=True)
            verify(root,sources)
            if side=='candidate':
                subprocess.run(['git','apply','--reverse',str(out/'stage203-candidate.patch')],cwd=root,check=True)
                verify(root,manifest['reference_sources'])
        # Exercise the exact builder's forward/reverse patch sequence on clean sources.
        for n,reverse,expected in [('stage203-baseline-observer.patch',False,'baseline_profile_sources'),
                                  ('stage203-baseline-observer.patch',True,'reference_sources'),
                                  ('stage203-candidate.patch',False,'candidate_sources'),('stage203-observer.patch',False,'profiled_sources')]:
            subprocess.run(['git','apply',*(['--reverse'] if reverse else []),str(out/n)],cwd=root,check=True)
            verify(root,manifest[expected])
    (out/'stage203-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    helpers=json.loads((previous/'helpers-sha256.json').read_text())
    for n in (*HELPERS,'stage203-sources.json','stage203-baseline.patch','stage203-baseline-observer.patch','stage203-candidate.patch','stage203-observer.patch'):helpers[n]=file_digest(out/n)
    (out/'helpers-sha256.json').write_text(json.dumps(helpers,indent=2)+'\n')
    for a,b in [('parallel-blocks-cluster-203.md','README.md'),('parallel-blocks-stage203.md','EXPERIMENT.md'),('parallel-blocks-stage203-local-results.json','LOCAL_VALIDATION.json')]:
        if (repo/'docs'/a).exists():shutil.copy2(repo/'docs'/a,out/b)
        else:(out/b).unlink(missing_ok=True)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    result={'stage':203,'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],
            'reference_patch_sha256':manifest['reference_patch_sha256'],'candidate_patch_sha256':manifest['candidate_patch_sha256'],
            'scope':'Seven fresh J5 builds, 18 solver runs, 83 ranks/two sockets; frozen Stage202 baseline and legacy; >=0.10 s/step prototype screen.'}
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('repo','out','stage202-kit'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.stage202_kit)
