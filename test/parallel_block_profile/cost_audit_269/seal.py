"""Seal only after four real source-bound observer replays and workflow tests."""
import json,shutil,sys,tarfile
from pathlib import Path
from common import save,sha,verify
from prepare import inventory
from cost_analysis import read_costs
from evidence import publication_evidence


def seal(kit,validation,tests):
    kit,validation,tests=map(lambda p:Path(p).resolve(strict=True),(kit,validation,tests))
    verify(kit)
    result=json.loads(validation.read_text());checks=json.loads(tests.read_text())
    cases={'main_observer','observer','checked','observer-oracle'}
    if result['status']!='passed' or result['ranks']!=4 or set(result['cases'])!=cases:raise ValueError('Require four real observer replays')
    if checks['status']!='passed' or checks['workflow_tests']!=9 or checks['observer_tests']!=1:raise ValueError('Require complete test evidence')
    if sha(kit/'source-files.json')!=result['source_manifest_sha256']:raise ValueError('Validated sources differ')
    builds_path=validation.parent/'builds.json';builds=json.loads(builds_path.read_text())
    if sha(builds_path)!=result['build_manifest_sha256'] or builds['status']!='passed' or len(builds['physics_objects'])!=23:raise ValueError('Build evidence differs')
    for kind,row in result['cases'].items():
        directory=validation.parent/kind;oracle=kind.endswith('-oracle');block=kind!='main_observer'
        if row['status']!='passed' or not row['numerical']['exact_observed'] or not row['step_states_exact'] or row['vtk_exact_members']!=62:raise ValueError('Local exact replay differs')
        if sha(directory/'run.log')!=row['log_sha256'] or read_costs(directory,4,block,oracle)!=row['costs']:raise ValueError('Local records changed')
        if {p.name:sha(p) for p in directory.glob('cost-rank-*.txt')}!=row['observer_file_sha256']:raise ValueError('Local raw files changed')
        row['publication']=publication_evidence((directory/'run.log').read_text(),'candidate' if block else 'main',oracle)
    here=Path(__file__).resolve().parent
    for p in here.iterdir():
        if p.is_file() and p.suffix in ('.py','.sh','.sbatch','.md','.f90'):shutil.copy2(p,kit/p.name)
    if (kit/'repeat.sh').exists():raise ValueError('No repeat launcher allowed')
    for name,expected in checks['tested_source_sha256'].items():
        if sha(kit/name)!=expected:raise ValueError('Tested helper changed: '+name)
    save(kit/'LOCAL_VALIDATION.json',dict(result,build_notes=builds,validation_sha256=sha(validation)))
    save(kit/'TEST_VALIDATION.json',checks)
    save(kit/'files-sha256.json',inventory(kit));verify(kit)
    archive=kit.with_suffix('.tar.gz')
    if archive.exists():raise ValueError('Fresh archive required')
    temporary=archive.with_suffix('.tmp')
    def keep(info):return None if '__pycache__' in Path(info.name).parts else info
    with tarfile.open(temporary,'w:gz',format=tarfile.PAX_FORMAT) as stream:stream.add(kit,arcname=kit.name,filter=keep)
    temporary.rename(archive)
    return dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,files=len(inventory(kit)),source_manifest_sha256=sha(kit/'source-files.json'),local_validation_sha256=sha(kit/'LOCAL_VALIDATION.json'),test_validation_sha256=sha(kit/'TEST_VALIDATION.json'))

if __name__=='__main__':print(json.dumps(seal(*sys.argv[1:]),indent=2))
