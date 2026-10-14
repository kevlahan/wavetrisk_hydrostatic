"""Seal only the source bytes qualified by complete active production replays."""
import argparse,json,shutil,tarfile,tempfile
from pathlib import Path
from common import save,sha,verify,source_scope,verify_packaged_grid
from prepare import inventory
from campaign_gate import check
from contract import evidence

def seal(root,kit,evaluation,tests):
    root,kit,evaluation,tests=map(lambda p:Path(p).resolve(strict=True),(root,kit,evaluation,tests))
    result=json.loads(evaluation.read_text());test=json.loads(tests.read_text())
    gate=check(root,result)
    if not gate['sealing_ready']:raise ValueError(gate['errors'])
    if test.get('status')!='passed' or test.get('tests_run',0)<10 or not test.get('packaged_import_graph_passed'):raise ValueError('Campaign helper tests have not passed')
    verify(kit);source_scope(kit);verify_packaged_grid(kit)
    manifests=json.loads((kit/'source-files.json').read_text())
    if set(result['qualified_source_sha256'])!=set(json.loads((kit/'overlay-files.json').read_text())):
        raise ValueError('Qualification must cover the entire candidate overlay')
    if result.get('qualified_build_source_sha256')!=manifests['candidate']:
        raise ValueError('Qualification build differs from the complete packaged source')
    for name,wanted in result['qualified_source_sha256'].items():
        if manifests['candidate'].get(name)!=wanted:raise ValueError('Kit differs from qualified source: '+name)
    for mode,record in result['full_step_local_runs'].items():
        path=Path(record['validation_json']);directory=path.parent
        if sha(path)!=record['validation_sha256']:raise ValueError('Local validation changed')
        raw=json.loads(path.read_text())
        if raw['status']!='passed' or raw['steps']!=12 or not raw['step_states_exact'] or not raw['numerical']['exact_observed'] or raw['vtk_exact_members']!=62:
            raise ValueError('Incomplete exact local replay: '+mode)
        if raw['oracles']!=(mode=='oracle'):
            raise ValueError('Local oracle selection differs: '+mode)
        if ('DEBUG=check' in raw['build_record']['command'])!=(mode=='checked'):
            raise ValueError('Local compiler mode differs: '+mode)
        if sha(directory/'run.log')!=raw['log_sha256'] or sha(directory/'climate')!=raw['binary_sha256']:
            raise ValueError('Local log/binary changed: '+mode)
        if raw.get('connected_contract')!=evidence(directory,4,True,domains=40):
            raise ValueError('Local raw RHS/RK ownership evidence differs: '+mode)
        if raw.get('qualified_source_sha256')!=result['qualified_source_sha256']:
            raise ValueError('Replay was not bound to the delivered solver source: '+mode)
    here=Path(__file__).resolve().parent
    for name,wanted in test['tested_source_sha256'].items():
        if sha(here/name)!=wanted:raise ValueError('Tested campaign helper changed: '+name)
        shutil.copyfile(here/name,kit/name)
    save(kit/'LOCAL_VALIDATION.json',result);save(kit/'TEST_VALIDATION.json',test)
    save(kit/'ENTRY_GATE.json',gate)
    spec=json.loads((kit/'campaign-spec.json').read_text())
    spec.update(status='sealed_for_user_submission',launch_ready=True,sealed=True,submitted=False)
    save(kit/'campaign-spec.json',spec)
    save(kit/'files-sha256.json',inventory(kit));verify(kit)
    archive=kit.with_suffix('.tar.gz')
    if archive.exists():raise ValueError('Fresh archive required')
    with tarfile.open(archive,'w:gz',format=tarfile.PAX_FORMAT) as stream:stream.add(kit,arcname=kit.name)
    # Delivery is checked from freshly extracted bytes, including the packaged
    # grid and source manifest. A successful tar creation alone is insufficient.
    with tempfile.TemporaryDirectory(prefix='rhs278-delivery-') as tmp:
        with tarfile.open(archive) as stream:stream.extractall(tmp,filter='data')
        fresh=Path(tmp)/kit.name;verify(fresh);source_scope(fresh);verify_packaged_grid(fresh)
        with tarfile.open(fresh/'sources.tar.gz') as stream:stream.extractall(Path(tmp)/'source',filter='data')
        for kind,files in manifests.items():
            for name,wanted in files.items():
                if sha(Path(tmp)/'source'/kind/name)!=wanted:raise ValueError('Fresh source extraction differs')
    return dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,files=len(inventory(kit)),
                launch_authorized=False,scope='Sealed campaign for user submission; cluster performance and promotion pending')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','kit','evaluation','tests'):p.add_argument('--'+name,type=Path,required=True)
    print(json.dumps(seal(**vars(p.parse_args())),indent=2))
