"""New serial candidate builds; source-bound controls reused, four MPI ranks only."""
import json,os,platform,shutil,subprocess,sys,tarfile
from pathlib import Path
from common import save,sha,verify,runtime
from build import objects
HERE=Path(__file__).resolve().parent
sys.path.insert(1,str(HERE.parent))
from run_local_230 import validate as validate_block,vtk_members
from capability_223 import coverage,numerical
from analyze import analyze
from cost_analysis import read_costs

def validate(kit,out,baseline,main):
    kit,baseline,main=[Path(p).resolve(strict=True) for p in (kit,baseline,main)]
    out=Path(out).resolve();out.mkdir(exist_ok=False)
    verify(kit);source=json.loads((kit/'source-files.json').read_text())
    refs={'stage266':baseline,'main':main}
    ref_notes={kind:json.loads((p/'builds.json').read_text()) for kind,p in refs.items()}
    for family,tree in (('stage266',baseline/'source'),('main',main/'sources/main')):
        for name,expected in source[family].items():
            if sha(tree/name)!=expected:raise ValueError('Reference source differs: '+family+'/'+name)
    physics=objects(baseline/'source')
    if physics!=objects(main/'sources/main'):raise ValueError('Reference physics differs')
    root=out/'sources';root.mkdir()
    with tarfile.open(kit/'sources.tar.gz') as stream:
        for kind in ('candidate',):
            members=[m for m in stream.getmembers() if m.name==kind or m.name.startswith(kind+'/')]
            stream.extractall(root,members=members,filter='data')
    notes=dict(status='building',host=platform.node(),level='J4',builds={},physics_objects=physics,source_manifest_sha256=sha(kit/'source-files.json'))
    for family,key in (('stage266','candidate'),('main','main')):
        row=ref_notes[family]['builds'][key].copy();binary=refs[family]/row['binary']
        if sha(binary)!=row['binary_sha256']:raise ValueError('Reference binary changed')
        row.update(binary=os.path.relpath(binary,out),reused_reference=True)
        notes['builds'][family]=row
    save(out/'builds.json',notes)
    for kind in ('candidate','checked'):
        family='candidate' if kind=='checked' else kind;tree=root/family
        if kind!='checked':
            dest=tree/'src/physics/simple_physics';shutil.rmtree(dest)
            shutil.copytree(baseline/'source/src/physics/simple_physics',dest,symlinks=True)
        (tree/'audit-flags.mk').write_text('FFLAGS += -ffp-contract=off\n')
        command=['make','-f','Makefile','-f','audit-flags.mk','-j1','TEST_CASE=climate','PARAM=param_J4','DEBUG='+('check' if kind=='checked' else 'false'),'BUILD_DIR=build-rk272-'+kind,'BIN_DIR=bin-rk272-'+kind]
        if family=='candidate':command+=['FP_CONTRACT=off','STAGE202=on','STAGE203=on']
        print('Local serial build: '+kind,flush=True)
        with (out/('build-'+kind+'.log')).open('x') as log:subprocess.run(command,cwd=tree,stdout=log,stderr=subprocess.STDOUT,check=True)
        if objects(tree)!=physics:raise ValueError('Shared physics changed')
        for name,expected in source[family].items():
            if sha(tree/name)!=expected:raise ValueError('Build changed sealed source: '+name)
        binary=tree/('bin-rk272-'+kind)/'climate'
        notes['builds'][kind]=dict(command=command,binary=str(binary.relative_to(out)),binary_sha256=sha(binary),physics_objects=physics,log_sha256=sha(out/('build-'+kind+'.log')))
        save(out/'builds.json',notes)
    notes['status']='passed';save(out/'builds.json',notes)
    return replay(kit,out,baseline,main,notes)

def replay(kit,out,baseline,main,notes):
    kit,out,baseline,main=map(lambda p:Path(p).resolve(strict=True),(kit,out,baseline,main))
    verify(kit)
    if notes['status']!='passed' or notes['level']!='J4' or notes['source_manifest_sha256']!=sha(kit/'source-files.json'):
        raise ValueError('Replay requires exact completed builds')
    for kind,row in notes['builds'].items():
        if sha(out/row['binary'])!=row['binary_sha256']:raise ValueError('Replay binary changed: '+kind)
    result=dict(status='running',ranks=4,cases={},source_manifest_sha256=sha(kit/'source-files.json'),build_manifest_sha256=sha(out/'builds.json'))
    fixture=kit/'fixture-J4'
    for kind,oracle in (('candidate',False),('checked',False),('candidate',True)):
        label=kind+('-oracle' if oracle else '');directory=out/label;directory.mkdir()
        for name in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copy2(fixture/name,directory/name)
        (directory/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
        binary=out/notes['builds'][kind]['binary'];shutil.copy2(binary,directory/'climate')
        env=runtime(oracle,observer=True);env['OMPI_MCA_mpi_yield_when_idle']='1'
        command=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
        print('Local four-rank replay: '+label,flush=True)
        with (directory/'run.log').open('x') as log:
            code=subprocess.run(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800).returncode
        if code:raise ValueError('Candidate replay failed: '+label)
        row=validate_block(directory,baseline/'opt',4,oracle,command,env)
        from evidence import publication_evidence
        row['publication']=publication_evidence((directory/'run.log').read_text(),'candidate',oracle)
        row['costs']=read_costs(directory,4,block=True,oracle=oracle,domains=40)
        result['cases'][label]=row;save(out/'local-validation.json',result)
    result.update(status='passed',limitation='Exact four-rank correctness and science exactness only; cluster scaling unmeasured.')
    save(out/'local-validation.json',result);print('PASS: optimized/checked/oracle candidate replays; exact source-bound Stage266 results',flush=True)
    return result

if __name__=='__main__':
    if len(sys.argv)==6 and sys.argv[5]=='--replay-only':
        kit,out,baseline,main=map(Path,sys.argv[1:5])
        replay(kit,out,baseline,main,json.loads((out/'builds.json').read_text()))
    else:validate(*sys.argv[1:])
