"""Package a qualified candidate for a later bbserv run; never submit or contact a host."""
import argparse,json,shutil,subprocess,sys,tarfile,tempfile
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import checkpoint_name
from experiment import SWITCHES
from run_local_220 import SELECTORS

def ignored(directory,names):
    return [n for n in names if n in {'.git','__pycache__','obj','lib','xml','outputs'} or
            (n=='driver' and (Path(directory)/n).is_file()) or n.startswith(('build-','bin-','._')) or n=='.DS_Store' or n.endswith(('.o','.mod','.so','.dylib','.a','.pyc','.log'))]

def prepare(tree,out):
    tree,out=tree.resolve(),out.resolve()
    gate=BASE/'docs/parallel-blocks-stage247-results.json';r=json.loads(gate.read_text())
    assert r['status']=='passed' and (BASE/r['frozen_candidate']).resolve()==tree
    if out.exists() or out.with_suffix('.tar.gz').exists():raise ValueError('Use a fresh package destination')
    out.mkdir(parents=True);source=out/'source';source.mkdir()
    for name in ['Makefile','README.md']:shutil.copy2(tree/name,source/name)
    for name in ['src','test']:shutil.copytree(tree/name,source/name,symlinks=False,ignore=ignored)
    for key in ['changed_sources','package_sources','resolved_sources']:
        for name,sha in r['build'][key].items():
            original=tree/name
            if key=='resolved_sources':
                assert digest(original.resolve())==sha,name
                content_sha=sha
            else:
                assert digest(original)==sha,name
                content_sha=digest(original.resolve()) if original.is_symlink() else sha
            assert digest(source/name)==content_sha,name
    fixture=BASE/'.stage241-validation/production-4';target=out/'fixture-J4';target.mkdir()
    shutil.copy2(fixture/'simple.in',target/'simple.in')
    seed=checkpoint_name((fixture/'simple.in').read_text());shutil.copy2(fixture/seed,target/seed)
    shutil.copytree(fixture/'grids',target/'grids',symlinks=False,ignore=ignored)
    shutil.copy2(BASE/'test/parallel_block_profile/prepare_physics_dependencies.sh',out/'prepare_physics_dependencies.sh')
    shutil.copy2(gate,out/'LOCAL_VALIDATION.json')
    command=[v for v in r['build']['commands']['production'] if not v.startswith(('BUILD_DIR=','BIN_DIR='))]
    command+=['BUILD_DIR=build-bbserv','BIN_DIR=bin-bbserv']
    save(out/'build-command.json',command)
    (out/'build.py').write_text('''"""Build the sealed candidate locally on the invoking host; never submit a job."""
import hashlib,json,platform,subprocess
from pathlib import Path
kit=Path(__file__).resolve().parent
for name,sha in json.loads((kit/'files-sha256.json').read_text()).items():
    assert hashlib.sha256((kit/name).read_bytes()).hexdigest()==sha,name
subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(kit/'source')],check=True)
cmd=json.loads((kit/'build-command.json').read_text())
versions={}
for name,args in [('compiler',['gfortran','--version']),('mpi',['mpif90','--version']),('netcdf',['nf-config','--all'])]:
    versions[name]=subprocess.check_output(args,text=True)
with (kit/'build.log').open('x') as log:subprocess.run(cmd,cwd=kit/'source',stdout=log,stderr=subprocess.STDOUT,check=True)
binary=kit/'source/bin-bbserv/climate'
(kit/'REMOTE_BUILD.json').write_text(json.dumps(dict(command=cmd,host=platform.node(),versions=versions,
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()),indent=2)+'\\n')
print(binary)
''')
    switches={**dict.fromkeys(SWITCHES,'0'),**SELECTORS,'WAVETRISK_VALIDATE_BLOCK_DYNAMICS':'0','WAVETRISK_VALIDATE_BLOCK_ADAPTATION':'0','WAVETRISK_VALIDATE_BLOCK_REMAP':'0',
              **{n:'1' for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']}}
    (out/'runtime.env').write_text('\n'.join('export '+k+'='+v for k,v in sorted(switches.items()))+'\n')
    (out/'README.md').write_text('''# Qualified RK candidate for bbserv\n\nThis package is prepared only; no remote build, job or performance test has run.\n\nLoad GNU Fortran, MPI and NetCDF-Fortran/C on bbserv, then run `python3 build.py`. The build is serial and records compiler/library details and the resulting binary hash. Source hashes are checked before building. Do not invoke the physics package clean target: it reformats its sources.\n\n`fixture-J4` contains the local restarted climate RK4/Simple-physics qualification input. Source `runtime.env` in a clean job environment before invoking the candidate. Select an allocation and matching reference run before timing; avoid node bb04. Rank placement, problem sizes, repetitions and baseline comparison belong to the next performance-test step. No local timing is a scaling result.\n\nThe candidate disables persistent legacy SOL/WAV views during stepping and retains native vector arithmetic/stencil/alias workspaces and fixed references. Output and restart use bounded views. `LOCAL_VALIDATION.json` binds exact step states and VTK outputs, the existing numerical checkpoint gate, checked builds, independent oracle and storage/epoch guards.\n''')
    files={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
    save(out/'files-sha256.json',files)
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as tar:tar.add(out,arcname=out.name)
    with tempfile.TemporaryDirectory(prefix='rk-bbserv-verify-') as temp:
        with tarfile.open(archive) as tar:tar.extractall(temp,filter='data')
        unpacked=Path(temp)/out.name
        for name,sha in files.items():assert digest(unpacked/name)==sha,name
    save(BASE/'docs/parallel-blocks-stage247-bbserv-kit.json',dict(status='prepared-not-run',path=str(out),archive=str(archive),
         archive_sha256=digest(archive),bytes=archive.stat().st_size,local_validation_sha256=digest(gate),
         commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),unpack_verified=True))
    print('PASS sealed bbserv candidate '+str(archive),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','out'):p.add_argument('--'+n,type=Path,required=True)
    prepare(**vars(p.parse_args()))
