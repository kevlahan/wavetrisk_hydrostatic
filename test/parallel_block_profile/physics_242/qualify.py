"""Run local solver cases strictly sequentially on the four-core host."""
import subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
ROOT=BASE/'.stage242-validation'
FIXTURE=BASE/'.stage241-validation/production-4'
HELPERS=Path(__file__).parent

def invoke(args,label):
    with (ROOT/(label+'-driver.log')).open('x') as log:
        subprocess.run([sys.executable,*map(str,args)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('PASS '+label,flush=True)
for name,variant,ranks,flags in [('oracle-4','oracle',4,['--oracles']),('production-4','production',4,['--full']),('production-41','production',41,['--full']),('control-4','control',4,['--control'])]:
    invoke([HELPERS/'run_physics.py','--binary',ROOT/'candidate'/('bin-'+variant)/'climate','--fixture',FIXTURE,'--out',ROOT/name,'--ranks',ranks,'--reference',BASE/('.stage241-validation/production-'+str(ranks)),*flags],name)
for mode in ('route','dependency','version','output','poison'):
    invoke([BASE/'test/parallel_block_profile/physics_240/reject_invalid.py','--tree',ROOT/'candidate','--fixture',FIXTURE,'--out',ROOT/('reject-'+mode),'--mode',mode],'reject-'+mode)

invoke([BASE/'test/parallel_block_profile/physics_239/reject_invalid.py','--tree',ROOT/'candidate','--fixture',FIXTURE,'--out',ROOT/'soil-2','--mode','soil'],'soil-2')

for mode in ('native-route','native-epoch','native-value','entry-poison','retained-owner','retained-alias'):
    invoke([HELPERS/'reject_invalid.py','--tree',ROOT/'candidate','--fixture',FIXTURE,'--out',ROOT/('reject-'+mode),'--mode',mode],'reject-'+mode)
