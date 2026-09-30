"""Run one local solver at a time; 41 ranks is oversubscribed correctness coverage."""
import subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
HELPERS=Path(__file__).resolve().parent
ROOT=BASE/'.stage243b-validation';FIXTURE=BASE/'.stage241-validation/production-4'
def invoke(command,name):
    with (ROOT/(name+'-driver.log')).open('x') as log:
        subprocess.run([sys.executable,*map(str,command)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('PASS '+name,flush=True)
for variant,name,ranks,flags in [
        ('control','control-4',4,['--control']),
        ('native','oracle-4',4,['--soil-check']),
        ('production','production-4',4,['--full']),
        ('production','production-41',41,['--full'])]:
    invoke([HELPERS/'run_soil.py','--binary',ROOT/'candidate'/('bin-'+variant)/'climate','--fixture',FIXTURE,
            '--out',ROOT/name,'--ranks',ranks,'--reference',BASE/('.stage242-validation/production-'+str(ranks)),*flags],name)
for mode in ['epoch','range','alias','poison','owner-readback','owner-readback-41']:
    name='reject-poison-exact-capacity' if mode=='poison' else 'reject-'+mode
    invoke([HELPERS/'reject_invalid.py','--tree',ROOT/'candidate','--fixture',FIXTURE,'--out',ROOT/name,'--mode',mode],name)
