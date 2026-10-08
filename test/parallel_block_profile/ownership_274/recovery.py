"""Recover only the known v1 postprocessing failure; preserve all original evidence."""
import copy,json,shutil
from pathlib import Path
from common import sha,verify,save,VALIDATIONS,J5_HASHES

EXPECTED_ERROR="ValueError('Duplicate Domain snapshot')"


def load_recovery(kit,previous):
    kit=Path(kit).resolve(strict=True);previous=Path(previous).resolve(strict=True)
    verify(kit)
    seal=json.loads((kit/'RECOVERY_SOURCE.json').read_text())
    if seal['source_manifest_sha256']!=sha(kit/'source-files.json'):
        raise ValueError('Recovery kit changed solver sources')
    notes=json.loads((previous/'builds.json').read_text())
    result=json.loads((previous/'results.json').read_text())
    if result.get('status')!='failed' or result.get('error')!=EXPECTED_ERROR:
        raise ValueError('Recovery requires the known v1 Duplicate Domain snapshot failure')
    if result.get('kit_manifest_sha256')!=seal['v1_manifest_sha256']:
        raise ValueError('Prior job is not the sealed v1 campaign')
    if (notes.get('status')!='passed' or notes.get('level')!='J5' or
            set(notes['builds'])!={'main','stage266','candidate','checked'} or
            len(notes['physics_objects'])!=23 or result.get('builds')!=notes or
            notes['source_manifest_sha256']!=sha(kit/'source-files.json')):
        raise ValueError('Cannot reuse incomplete or different builds')
    expected=[f'J5-n{n}-validation-{i}-{kind}' for i,(n,kind,oracle) in enumerate(VALIDATIONS[:3],1)]
    if result['schedule']!=expected or set(result['runs'])!=set(expected):
        raise ValueError('Recovery expects exactly the three completed solver runs from v1')
    source=json.loads((kit/'source-files.json').read_text())
    for family,files in source.items():
        for name,wanted in files.items():
            if sha(previous/'sources'/family/name)!=wanted:raise ValueError('Prior source changed: '+family+'/'+name)
    for name,wanted in J5_HASHES.items():
        if sha(previous/'input'/name)!=wanted:raise ValueError('Prior fixture changed: '+name)
    for kind,row in notes['builds'].items():
        binary=(previous/row['binary']).resolve(strict=True)
        if not binary.is_relative_to(previous) or sha(binary)!=row['binary_sha256']:
            raise ValueError('Prior binary changed: '+kind)
        if sha(previous/('build-'+kind+'.log'))!=row['log_sha256'] or row['physics_objects']!=notes['physics_objects']:
            raise ValueError('Prior build evidence changed: '+kind)
    for label,(n,kind,oracle) in zip(expected,VALIDATIONS[:3]):
        row=result['runs'][label];directory=previous/label
        if row['returncode']!=0 or (row['ranks'],row['kind'],row['oracle'])!=(n,kind,oracle):
            raise ValueError('Prior solver did not finish successfully: '+label)
        if sha(directory/'climate')!=notes['builds'][kind]['binary_sha256'] or sha(directory/'simple.in')!=J5_HASHES['simple.in']:
            raise ValueError('Prior run identity changed: '+label)
        if not (directory/'run.log').is_file():raise ValueError('Missing prior solver log')
        if 'log_sha256' in row and sha(directory/'run.log')!=row['log_sha256']:
            raise ValueError('Prior completed log changed: '+label)
    return dict(previous=str(previous),notes=notes,result=result,labels=expected,
                results_sha256=sha(previous/'results.json'),builds_sha256=sha(previous/'builds.json'))


def reuse_builds(kit,out,record):
    out=Path(out);previous=Path(record['previous']);notes=copy.deepcopy(record['notes'])
    for kind,row in notes['builds'].items():
        destination=out/'reused-binaries'/kind/'climate';destination.parent.mkdir(parents=True)
        shutil.copy2(previous/row['binary'],destination)
        if sha(destination)!=row['binary_sha256']:raise ValueError('Copied binary differs')
        row['binary']=str(destination.relative_to(out));row['reused_from_job']=record['result']['job_id']
        shutil.copy2(previous/('build-'+kind+'.log'),out/('build-'+kind+'.log'))
    shutil.copy2(previous/'results.json',out/'recovery-source-results.json')
    shutil.copy2(previous/'builds.json',out/'recovery-source-builds.json')
    notes['reused_builds']=True
    save(out/'builds.json',notes)
    print('Reusing all four verified binaries; no compilation required',flush=True)
    return notes


def submit(kit):
    import os,subprocess,tempfile
    import control
    kit=Path(kit).resolve(strict=True);home=Path.home()
    prior=json.loads(control.state_path().read_text())
    if prior.get('phase')=='recovery':raise ValueError('Recovery already submitted; use status.sh and collect.sh')
    if control.queued(prior['job']):raise ValueError('Prior job is still active')
    control.preflight(kit,Path(prior['fixture']))
    record=load_recovery(kit,prior['out'])
    if str(prior['job'])!=str(record['result']['job_id']):raise ValueError('Latest job/result identity differs')
    lock=home/'.rk274-launch.lock'
    try:lock.mkdir()
    except FileExistsError:raise ValueError('Another Stage274 submission is in progress') from None
    try:
        if json.loads(control.state_path().read_text())!=prior:raise ValueError('Latest job changed during recovery preflight')
        work=Path(tempfile.mkdtemp(prefix='rk274-recovery-',dir=home));out=work/'results'
        response=subprocess.check_output(['sbatch','--parsable','--chdir='+str(home),
            '--output='+str(home/'rk274-%j.log'),str(kit/'recovery.sbatch'),
            str(kit),prior['fixture'],str(out),prior['out']],text=True).strip().split(';')[0]
        if not response.isdigit():raise ValueError('Unexpected sbatch response: '+response)
        state=dict(phase='recovery',job=int(response),kit=str(kit),out=str(out),fixture=prior['fixture'],
                   log=str(home/f'rk274-{response}.log'),previous_job=prior['job'],previous_out=prior['out'])
        save(home/f'rk274-job-{response}.json',state);save(control.state_path(),state)
        print('Submitted Stage274 recovery job: '+response)
        print('Reuse four builds and three completed solver runs; run five remaining cases.')
        print('After completion: bash '+str(kit/'collect.sh'))
    finally:lock.rmdir()

if __name__=='__main__':
    import sys
    if len(sys.argv)!=2:raise SystemExit('Usage: recovery.py KIT')
    submit(sys.argv[1])
