#!/usr/bin/env python3
"""Freeze the tested Stage202/203 sources; replace observers only."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from build_restriction_203 import verify
from instrument_operations_204 import instrument,REGIONS
BASE=Path(__file__).parent
HELPERS=('build_restriction_204.py','submit_restriction_204.py','build_204.sbatch','cluster_204.sbatch','cluster_204.py',
         'analyze_operations_204.py','instrument_operations_204.py','collect_cluster_logs_204.py')

def upgrade(root):
    inventory=instrument(root)
    (root/'src/workflow_profile_195.f90').write_bytes((BASE/'workflow_profile_204.f90').read_bytes())
    p=root/'Makefile';s=p.read_text().replace('FFLAGS += -DSTAGE203_LEDGER\n','')
    p.write_text(s)
    return inventory

def prepare(repo,out,previous):
    repo=repo.resolve();out=out.resolve();previous=previous.resolve()
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit required')
    old=json.loads((previous/'stage203-sources.json').read_text())
    for n,h in json.loads((previous/'helpers-sha256.json').read_text()).items():
        if file_digest(previous/n)!=h:raise ValueError('Previous helper changed '+n)
    shutil.copytree(previous,out,ignore=shutil.ignore_patterns('__pycache__'));manifest=dict(old);inventory={}
    for n in HELPERS:shutil.copy2(BASE/n,out/n)
    for name in ('baseline','candidate'):shutil.copy2(previous/f'stage203-{name}.patch',out/f'stage204-{name}.patch')
    for side in ('baseline','candidate','legacy'):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            if side=='legacy':
                with tarfile.open(previous/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(previous/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
                source='legacy_io_sources';profile='legacy_profile_sources';oldpatch='legacy-profile-201.patch';target='legacy-profile-204.patch';hashkey='legacy_profile_patch_sha256'
            else:
                with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(previous/'stage203-baseline.patch')],cwd=root,check=True)
                if side=='candidate':subprocess.run(['git','apply',str(previous/'stage203-candidate.patch')],cwd=root,check=True)
                source='candidate_sources' if side=='candidate' else 'reference_sources'
                profile='profiled_sources' if side=='candidate' else 'baseline_profile_sources'
                oldpatch='stage203-observer.patch' if side=='candidate' else 'stage203-baseline-observer.patch'
                target=oldpatch.replace('203','204');hashkey='observer_patch_sha256' if side=='candidate' else 'baseline_observer_patch_sha256'
            verify(root,old[source]);before={n:(root/n).read_bytes() for n in old[source] if not (root/n).is_symlink()}
            subprocess.run(['git','apply',str(previous/oldpatch)],cwd=root,check=True);verify(root,old[profile])
            inventory[side]=upgrade(root)
            names=set(old[profile])|{n for n in inventory[side]}
            manifest[profile]={n:file_digest(root/n) for n in sorted(names)}
            (out/target).write_text(patch(before,root,sorted(n for n in names if manifest[profile][n]!=old[source].get(n))))
            manifest[hashkey]=file_digest(out/target)
            subprocess.run(['git','apply','--reverse',str(out/target)],cwd=root,check=True);verify(root,old[source])
            subprocess.run(['git','apply',str(out/target)],cwd=root,check=True);verify(root,manifest[profile])
    manifest.update(stage=204,diagnosis_only=True,observer_sha256=file_digest(BASE/'workflow_profile_204.f90'))
    (out/'stage204-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'operation-sites-204.json').write_text(json.dumps({'regions':REGIONS,'implementations':inventory},indent=2)+'\n')
    hashes=json.loads((previous/'helpers-sha256.json').read_text())
    for n in (*HELPERS,'stage204-sources.json','operation-sites-204.json','stage204-baseline.patch','stage204-candidate.patch','stage204-baseline-observer.patch','stage204-observer.patch','legacy-profile-204.patch'):hashes[n]=file_digest(out/n)
    (out/'helpers-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    for a,b in [('parallel-blocks-cluster-204.md','README.md'),('parallel-blocks-stage204.md','EXPERIMENT.md'),('parallel-blocks-stage204-local-results.json','LOCAL_VALIDATION.json')]:
        if (repo/'docs'/a).exists():shutil.copy2(repo/'docs'/a,out/b)
        else:(out/b).unlink(missing_ok=True)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    result={'stage':204,'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':old['baseline_commit'],'scope':'Frozen Stage202/203/legacy; seven builds, 20 runs, clocks only; no numerical changes.'}
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('repo','out','stage203-kit'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.stage203_kit)
