"""Verify width-independent route metadata and exact variable-width exchange."""
import argparse,hashlib,json,subprocess
from pathlib import Path

def run(output):
    root=Path(__file__).resolve().parents[3]
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    unit=Path(__file__).with_name('width_route_unit.f90')
    sources=[root/'src/kind.f90',root/'src/parallel_block_profile.f90',root/'src/parallel_block_rhs_exchange.f90',unit]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report=dict(status='running',source_sha256={str(p.relative_to(root)):sha(p) for p in sources},modes={})
    def save():(output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    try:
        for mode in ('checked','optimized'):
            directory=output/mode;directory.mkdir()
            build=['mpifort','-ffree-line-length-132','-ffp-contract=off','-Wall','-Wextra','-Werror','-J'+str(directory),'-I'+str(directory)]
            build+=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if mode=='checked' else ['-O2']
            build += [str(p) for p in sources]+['-o',str(directory/'unit')]
            with (directory/'build.log').open('x') as f:subprocess.run(build,cwd=directory,stdout=f,stderr=subprocess.STDOUT,check=True)
            runs=[]
            for ranks,case,marker in [(1,'normal','PASS variable widths'),(4,'normal','PASS variable widths'),(4,'variable-mismatch','inconsistent distributed plan epoch'),(4,'bad-width','runtime workspace extent')]:
                cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(ranks),str(directory/'unit'),case]
                result=subprocess.run(cmd,cwd=directory,text=True,capture_output=True,timeout=60)
                log=directory/f'{ranks}-{case}.log';log.write_text(result.stdout+result.stderr)
                ok=(result.returncode==0 if case=='normal' else result.returncode!=0) and marker in log.read_text()
                runs.append(dict(ranks=ranks,case=case,passed=ok,returncode=result.returncode,log_sha256=sha(log)))
                if not ok:raise ValueError((case,log.read_text()))
            report['modes'][mode]=dict(build=build,binary_sha256=sha(directory/'unit'),runs=runs);save()
        report['status']='passed';save();print('PASS width-independent route tests')
    except BaseException as error:
        report['status']='failed';report['error']=repr(error);save();raise

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
