#!/usr/bin/env python3
"""Freeze Stage202 against the exact Stage201 production sources and legacy."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from instrument_workflow_195 import instrument
from build_restriction_202 import verify
HELPERS=('build_restriction_202.py','submit_restriction_202.py','build_202.sbatch','cluster_202.sbatch','cluster_202.py',
         'analyze_batching_202.py','preflight_202.py','collect_cluster_logs_202.py')
CHANGES=('Makefile','src/parallel_block_scalar_storage.f90','src/parallel_block_mpi.f90','src/parallel_block_inverse.f90','src/parallel_block.f90')

def observer(root,base):
    instrument(root,True)
    (root/'src/workflow_profile_195.f90').write_bytes((base/'workflow_profile_201.f90').read_bytes())
    (root/'src/workflow_clock_199.c').write_bytes((base/'workflow_clock_201.c').read_bytes())
    p=root/'Makefile';s=p.read_text();old='OBJ := $(addprefix $(BUILD_DIR)/,$(SRC:.f90=.o))'
    if s.count(old)!=1:raise ValueError('Missing object anchor')
    s=s.replace(old,old+' $(BUILD_DIR)/workflow_clock_199.o')
    # Append an explicit C rule; using mpicc avoids the system-cc include issue.
    s+='\nMPICC ?= mpicc\n$(BUILD_DIR)/workflow_clock_199.o: $(SRC_DIR)/workflow_clock_199.c | dirs\n\t$(MPICC) -O2 -c $< -o $@\n'
    p.write_text(s)

def prepare(repo,out,previous):
    repo=repo.resolve();out=out.resolve();previous=previous.resolve();base=Path(__file__).parent
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit/archive required')
    old=json.loads((previous/'stage201-sources.json').read_text())
    # The only production differences allowed are the explicit prototype files.
    for n,h in old['candidate_sources'].items():
        if n not in CHANGES and file_digest(repo/n)!=h:raise ValueError('Unrelated source difference: '+n)
    for n,h in json.loads((previous/'helpers-sha256.json').read_text()).items():
        if file_digest(previous/n)!=h:raise ValueError('Changed previous helper: '+n)
    if file_digest(previous/'stage201-candidate.patch')!=old['candidate_patch_sha256']:raise ValueError('Previous numerical patch differs')
    shutil.copytree(previous,out);manifest=dict(old)
    for n in HELPERS:shutil.copy2(base/n,out/n)
    shutil.copy2(previous/'stage201-candidate.patch',out/'stage202-baseline.patch')
    manifest.update(stage=202,reference_sources=old['candidate_sources'],reference_patch_sha256=file_digest(out/'stage202-baseline.patch'),diagnosis_only=False)
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
        subprocess.run(['git','apply',str(out/'stage202-baseline.patch')],cwd=root,check=True)
        verify(root,manifest['reference_sources'])
        original={n:(root/n).read_bytes() for n in old['candidate_sources'] if not (root/n).is_symlink()}
        for n in CHANGES:shutil.copy2(repo/n,root/n)
        manifest['candidate_sources']={n:file_digest(root/n) for n in old['candidate_sources']}
        (out/'stage202-candidate.patch').write_text(patch(original,root,CHANGES))
        manifest['candidate_patch_sha256']=file_digest(out/'stage202-candidate.patch')
        for side,key,patchname,hashkey in [('candidate','profiled_sources','stage202-observer.patch','observer_patch_sha256'),
                                         ('baseline','baseline_profile_sources','stage202-baseline-observer.patch','baseline_observer_patch_sha256')]:
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
                subprocess.run(['git','apply','--reverse',str(out/'stage202-candidate.patch')],cwd=root,check=True)
                verify(root,manifest['reference_sources'])
        # Exercise the exact builder's forward/reverse patch sequence on clean sources.
        for n,reverse,expected in [('stage202-baseline-observer.patch',False,'baseline_profile_sources'),
                                  ('stage202-baseline-observer.patch',True,'reference_sources'),
                                  ('stage202-candidate.patch',False,'candidate_sources'),('stage202-observer.patch',False,'profiled_sources')]:
            subprocess.run(['git','apply',*(['--reverse'] if reverse else []),str(out/n)],cwd=root,check=True)
            verify(root,manifest[expected])
    (out/'stage202-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    helpers=json.loads((previous/'helpers-sha256.json').read_text())
    for n in (*HELPERS,'stage202-sources.json','stage202-baseline.patch','stage202-baseline-observer.patch','stage202-candidate.patch','stage202-observer.patch'):helpers[n]=file_digest(out/n)
    (out/'helpers-sha256.json').write_text(json.dumps(helpers,indent=2)+'\n')
    for a,b in [('parallel-blocks-cluster-202.md','README.md'),('parallel-blocks-stage202.md','EXPERIMENT.md'),('parallel-blocks-stage202-local-results.json','LOCAL_VALIDATION.json')]:
        if (repo/'docs'/a).exists():shutil.copy2(repo/'docs'/a,out/b)
        else:(out/b).unlink(missing_ok=True)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    result={'stage':202,'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],
            'reference_patch_sha256':manifest['reference_patch_sha256'],'candidate_patch_sha256':manifest['candidate_patch_sha256'],
            'scope':'Six fresh J5 builds, 17 solver runs, 83 ranks/two sockets; frozen Stage201 baseline and legacy; >=0.10 s/step prototype screen.'}
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('repo','out','stage201-kit'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.stage201_kit)
