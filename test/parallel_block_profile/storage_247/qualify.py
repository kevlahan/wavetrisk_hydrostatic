"""Sequential local qualification; oversubscribed ranks are correctness coverage only."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3];H=Path(__file__).resolve().parent

def run(tree,root,checked_case):
    tree,root,checked_case=(p.resolve() for p in (tree,root,checked_case))
    sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
    from storage_247.run_storage import review
    from capability_223 import digest
    r=json.loads((checked_case/'lifecycle-results.json').read_text())
    assert r['status']=='passed' and r['steps']==12 and r['ranks']==4
    assert r['binary_sha256']==digest(tree/'bin-storage/climate')
    review(checked_case,4,False)
    def invoke(name,args):
        with (root/(name+'-driver.log')).open('x') as f:
            subprocess.run([sys.executable,*map(str,args)],cwd=BASE,stdout=f,stderr=subprocess.STDOUT,check=True)
        print('PASS '+name,flush=True)
    invoke('remaining-builds',[H/'build_storage.py','--reference',BASE/'.stage246-validation/candidate','--out',tree,'--variants','production','oracle','control'])
    for name,variant,ranks,flags,reference in [
        ('oracle4-4','oracle',4,['--oracles','--checks','--views'],BASE/'.stage246-validation/production-4'),
        ('control4-4','control',4,['--control','--views'],BASE/'.stage246-validation/production-4'),
        ('production4-4','production',4,['--full'],BASE/'.stage246-validation/production-4'),
        ('production4-41','production',41,['--full'],BASE/'.stage244-validation/production3-41')]:
        invoke(name,[H/'run_storage.py','--binary',tree/('bin-'+variant)/'climate','--fixture',BASE/'.stage241-validation/production-4',
                     '--out',root/name,'--ranks',ranks,'--reference',reference,*flags])
    for mode in ['lifecycle','coefficients','mirror']:
        invoke('reject-'+mode,[H/'reject_invalid.py','--tree',tree,'--fixture',BASE/'.stage241-validation/production-4',
                              '--out',root/('reject-'+mode),'--mode',mode])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','root','checked-case'):p.add_argument('--'+n,type=Path,required=True)
    run(**vars(p.parse_args()))
