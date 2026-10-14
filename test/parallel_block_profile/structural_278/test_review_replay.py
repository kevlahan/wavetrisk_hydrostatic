"""Workflow test using copied 55503 raw reports, never repair science/performance evidence.
Run with: python3 test_review_replay.py PREVIOUS_EXTRACTED_REVIEW
"""
from pathlib import Path
import copy,json,shutil,sys,tarfile,tempfile
from common import save,sha,VALIDATIONS
from diagnostic import schedule,summarize_diagnostic
from review import review
K=Path(__file__).resolve().parent

def exercise(fixture):
    fixture=Path(fixture);old=json.loads((fixture/'results/results.json').read_text())
    with tempfile.TemporaryDirectory(prefix='rk278-review-replay-') as tmp:
        root=Path(tmp)/'review';out=root/'results';out.mkdir(parents=True);(root/'kit').mkdir()
        for n in ('IDENTITY.json','files-sha256.json','source-files.json','LOCAL_VALIDATION.json','TEST_VALIDATION.json','ENTRY_GATE.json','overlay-files.json','parent-source-files.json'):
            shutil.copyfile(K/n,root/'kit'/n)
        shutil.copyfile(fixture/'accounting.json',root/'accounting.json');shutil.copytree(fixture/'results/input',out/'input')
        for p in (fixture/'results').glob('placement-*.json'):shutil.copyfile(p,out/p.name)
        result=copy.deepcopy(old);result.update(revision='structural-v1',identity=json.loads((K/'IDENTITY.json').read_text()),kit_manifest_sha256=sha(K/'files-sha256.json'),schedule=[],runs={},gates={})
        result.pop('fine_decisions',None)
        notes=copy.deepcopy(old['builds']);notes['builds']['previous']=notes['builds'].pop('v5')
        for k,r in notes['builds'].items():
            src='v5' if k=='previous' else k;shutil.copyfile(fixture/'results'/f'build-{src}.log',out/f'build-{k}.log')
        notes.update(identity=result['identity'],source_manifest_sha256=sha(K/'source-files.json'));save(out/'builds.json',notes);result['builds']=notes
        specs=[('validation',n,k,'off',o) for n,k,o in VALIDATIONS]+[('measurement',n,k,p,False) for n,k,p in schedule()]
        for index,(phase,n,kind,mode,oracle) in enumerate(specs):
            sourcekind='candidate' if kind=='previous' else kind
            name=next(name for name,r in old['runs'].items() if (r['phase'],r['ranks'],r['kind'],r['profile'],r['oracle'])==(phase,n,sourcekind,mode,oracle))
            row=copy.deepcopy(old['runs'][name]);row['kind']=kind
            # Bind the unchanged raw fixture to its actual source binary in this mock.
            if kind=='previous':notes['builds']['previous']['binary_sha256']=row['binary_sha256']
            label=f'replay-{index}-{n}-{kind}-{mode}';shutil.copytree(fixture/'results'/name,out/label)
            result['schedule'].append(label);result['runs'][label]=row;result['gates'][label]=copy.deepcopy(old['gates'][name])
        save(out/'builds.json',notes);result['comparison']=summarize_diagnostic(list(result['runs'].values()))
        def archive():
            save(out/'results.json',result)
            save(root/'files-sha256.json',{str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*')) if p.is_file() and p!=root/'files-sha256.json'})
            a=Path(tmp)/'synthetic-review.tar.gz'
            with tarfile.open(a,'w:gz') as t:t.add(root,arcname='review')
            return a
        report=review(archive());assert report['status']=='verified'
        result['comparison']['by_ranks']['160']['gap_to_stage266_seconds']+=.1
        try:review(archive())
        except ValueError as e:assert 'Timing comparison differs' in str(e)
        else:raise AssertionError('Tampered derived performance accepted')
        return dict(status='passed',cases=2,raw_fixture_job=old['job_id'],runs=48,scope='Synthetic workflow: reused old raw reports to exercise new schedule/reviewer and summary-tamper rejection. Does not test repair science or performance.')
if __name__=='__main__':print(json.dumps(exercise(sys.argv[1]),indent=2))
