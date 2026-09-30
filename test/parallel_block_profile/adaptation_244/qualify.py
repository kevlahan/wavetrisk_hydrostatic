"""Run one local solver at a time; 41 ranks is oversubscribed correctness coverage."""
import subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
HELPERS=Path(__file__).resolve().parent
ROOT=BASE/'.stage244-validation';FIXTURE=BASE/'.stage241-validation/production-4'
def invoke(command,name):
    with (ROOT/(name+'-driver.log')).open('x') as log:
        subprocess.run([sys.executable,*map(str,command)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('PASS '+name,flush=True)
for variant,name,ranks,flags in [
        ('control','control3-4',4,['--control']),
        ('native','oracle3-4',4,['--checks']),
        ('production','production3-4',4,['--full']),
        ('production','production3-41',41,['--full']),
        ('oracle','rk-oracle3-4',4,['--oracles','--checks'])]:
    if name=='production3-4' and (ROOT/name/'adaptation-results.json').exists():
        import json,hashlib
        r=json.loads((ROOT/name/'adaptation-results.json').read_text())
        assert r['status']=='passed' and r['steps']==12
        assert r['binary_sha256']==hashlib.sha256((ROOT/'candidate3'/'bin-production'/'climate').read_bytes()).hexdigest()
        print('PASS existing '+name,flush=True)
        continue
    invoke([HELPERS/'run_adaptation.py','--binary',ROOT/'candidate3'/('bin-'+variant)/'climate','--fixture',FIXTURE,
            '--out',ROOT/name,'--ranks',ranks,'--reference',BASE/('.stage243b-validation/production-'+str(ranks)),*flags],name)
for mode in ['epoch','generation','poison','owner-readback','owner-readback-41']:
    name='reject-'+mode
    invoke([HELPERS/'reject_invalid.py','--tree',ROOT/'candidate3','--fixture',FIXTURE,'--out',ROOT/name,'--mode',mode],name)
