"""Seal only after real source-bound optimized/checked/oracle diagnostic replays."""
import json,shutil,sys,tarfile
from pathlib import Path
from common import save,sha,verify,source_scope,verify_packaged_grid
from prepare import inventory
from evidence import publication_evidence
from contract import evidence


def seal(kit,validation,tests):
    kit,validation,tests=map(lambda p:Path(p).resolve(strict=True),(kit,validation,tests))
    verify(kit);source_scope(kit);verify_packaged_grid(kit)
    result=json.loads(validation.read_text());checks=json.loads(tests.read_text())
    if result['status']!='passed' or result['ranks']!=4 or set(result['cases'])!={'candidate','checked','candidate-oracle','candidate-production'}:
        raise ValueError('Require four real source-bound production and diagnostic replays')
    if checks['status']!='passed' or checks['workflow_tests']!=11 or checks['contract_tests']!=4:
        raise ValueError('Require workflow and contract test evidence')
    if result['source_manifest_sha256']!=sha(kit/'source-files.json'):raise ValueError('Tested source differs')
    bp=validation.parent/'builds.json';builds=json.loads(bp.read_text())
    if sha(bp)!=result['build_manifest_sha256'] or builds['status']!='passed' or len(builds['physics_objects'])!=23:
        raise ValueError('Build evidence differs')
    for kind,row in result['cases'].items():
        directory=validation.parent/kind;oracle=kind.endswith('-oracle')
        check_assignments=kind!='candidate-production'
        if check_assignments and (not row['contract']['assignment_gate_passed'] or not row['contract']['reference_gate_passed']):
            raise ValueError('Local primitive assignment gate failed')
        if row['contract']!=evidence(directory,4,check_assignments,domains=40):raise ValueError('Raw contract evidence differs')
        marker='Stage277 cached fused primitive = on; assignment checks = '+('T' if check_assignments else 'F')
        if (directory/'run.log').read_text().count(marker)!=1:raise ValueError('Local primitive mode marker differs')
        if row['status']!='passed' or not row['numerical']['exact_observed'] or not row['step_states_exact'] or row['vtk_exact_members']!=62:
            raise ValueError('Local exact replay differs')
        if sha(directory/'run.log')!=row['log_sha256']:raise ValueError('Local raw log changed')
        if publication_evidence((directory/'run.log').read_text(),'candidate',oracle)!=row['publication']:
            raise ValueError('Installation/oracle evidence differs')
        binary=validation.parent/builds['builds']['candidate' if kind.startswith('candidate') else kind]['binary']
        if sha(binary)!=row['binary_sha256']:raise ValueError('Replayed binary changed')
        row['contract_file_sha256']={p.name:sha(p) for p in sorted(list(directory.glob('primitive277-rank-*.txt')))}
    here=Path(__file__).resolve().parent
    for p in here.iterdir():
        if p.is_file() and p.suffix in ('.py','.sh','.sbatch','.md','.f90','.json'):shutil.copy2(p,kit/p.name)
    if (kit/'repeat.sh').exists():raise ValueError('No repeat launcher allowed')
    for name,expected in checks['tested_source_sha256'].items():
        if sha(kit/name)!=expected:raise ValueError('Tested helper changed: '+name)
    save(kit/'LOCAL_VALIDATION.json',dict(result,build_notes=builds,validation_sha256=sha(validation)))
    save(kit/'TEST_VALIDATION.json',checks)
    save(kit/'files-sha256.json',inventory(kit));verify(kit)
    archive=kit.with_suffix('.tar.gz')
    if archive.exists():raise ValueError('Fresh archive required')
    temporary=archive.with_suffix('.tmp')
    def keep(info):
        if '__pycache__' in Path(info.name).parts:return None
        info.pax_headers={k:v for k,v in info.pax_headers.items() if 'xattr' not in k.lower()}
        return info
    with tarfile.open(temporary,'w:gz',format=tarfile.PAX_FORMAT) as stream:stream.add(kit,arcname=kit.name,filter=keep)
    temporary.rename(archive)
    return dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,files=len(inventory(kit)),
                source_manifest_sha256=sha(kit/'source-files.json'),local_validation_sha256=sha(kit/'LOCAL_VALIDATION.json'),test_validation_sha256=sha(kit/'TEST_VALIDATION.json'))

if __name__=='__main__':print(json.dumps(seal(*sys.argv[1:]),indent=2))
