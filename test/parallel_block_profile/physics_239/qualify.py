"""Run local solver cases strictly sequentially on the four-core host."""
import subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
ROOT=BASE/'.stage239-validation'
FIXTURE=BASE/'.stage238-validation/production-4'
HELPERS=Path(__file__).parent

def invoke(args,label):
    with (ROOT/(label+'-driver.log')).open('x') as log:
        subprocess.run([sys.executable,*map(str,args)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('PASS '+label,flush=True)
for name,variant,ranks,flags in [('oracle-4','oracle',4,['--oracles']),('production-4','production',4,['--full']),('production-41','production',41,['--full']),('control-4','control',4,['--control'])]:
    invoke([HELPERS/'run_physics.py','--binary',ROOT/'candidate'/('bin-'+variant)/'climate','--fixture',FIXTURE,'--out',ROOT/name,'--ranks',ranks,'--reference',BASE/('.stage238-validation/production-'+str(ranks)),*flags],name)
for mode in ('route','order','version','output','soil'):
    invoke([HELPERS/'reject_invalid.py','--tree',ROOT/'candidate','--fixture',FIXTURE,'--out',ROOT/('reject-'+mode),'--mode',mode],'reject-'+mode)
