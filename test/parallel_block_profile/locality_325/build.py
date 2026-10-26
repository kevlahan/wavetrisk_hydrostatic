"""Two fresh serial candidate builds; native compiler/runtime gates precede solvers."""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess,tarfile,time
from common import sha,save,verify

def utc():return datetime.now(timezone.utc).isoformat()
def stamp(prefix):return {prefix+'_utc':utc(),prefix+'_monotonic_ns':time.monotonic_ns()}
def run_logged_build(command,tree,log_path):
 try:
  with Path(log_path).open('x') as log:subprocess.run(command,cwd=tree,stdout=log,stderr=subprocess.STDOUT,check=True)
 except subprocess.CalledProcessError:
  print('BUILD FAILED: '+str(log_path),flush=True)
  print('\n'.join(Path(log_path).read_text(errors='replace').splitlines()[-60:]),flush=True);raise

def extract_sources(kit,root):
 manifest=json.loads((kit/'source-files.json').read_text())
 if set(manifest)!={'candidate'}:raise ValueError('Only isolated candidate sources are built')
 root.mkdir()
 with tarfile.open(kit/'sources.tar.gz') as t:
  names=[]
  for m in t.getmembers():
   rel=Path(m.name)
   if rel.is_absolute() or '..' in rel.parts or not rel.parts or rel.parts[0]!='candidate' or not(m.isfile() or m.isdir() or m.issym()):raise ValueError('Unsafe source archive')
   if m.name in names:raise ValueError('Duplicate source archive member')
   names.append(m.name)
   if m.issym() and (Path(m.linkname).is_absolute() or not (root/rel.parent/m.linkname).resolve().is_relative_to((root/'candidate').resolve())):raise ValueError('Unsafe source link')
  t.extractall(root,filter='data')
 tree=root/'candidate'
 actual={str(p.relative_to(tree)) for p in tree.rglob('*') if p.is_file() or p.is_symlink()}
 if actual!=set(manifest['candidate']):raise ValueError('Source inventory differs')
 for name,h in manifest['candidate'].items():
  if sha(tree/name)!=h:raise ValueError('Source differs: '+name)
 return tree,manifest

def build(kit,out):
 from build_environment313 import controlled_environment
 with controlled_environment(kit,out):return build_in_environment(Path(kit),Path(out))

def build_in_environment(kit,out):
 verify(kit)
 record=dict(status='building',**stamp('started'),source_manifest_sha256=sha(kit/'source-files.json'),builds={},build_order=[])
 save(out/'builds.json',record)
 try:
  from toolchain313 import check,linked
  toolchain=check(kit,out);record['toolchain']=toolchain;save(out/'builds.json',record)
  tree,manifest=extract_sources(kit,out/'sources')
  from callback_gate import check as callback_check
  record['callback_gate']=callback_check(kit,tree,out,compiler=toolchain['wrapper_resolved'],native=True)
  if record['callback_gate'].get('status')!='passed' or len(record['callback_gate'].get('cases',[]))!=17:raise ValueError('Require17 native callback cases')
  subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(tree)],check=True)
  for kind,debug in [('checked','check'),('candidate','false')]:
   flags=tree/'stage325-flags.mk';flags.write_text('FFLAGS += -ffp-contract=off\n')
   command=['make','-f','Makefile','-f',flags.name,'-j1','TEST_CASE=climate','PARAM=param_J5','DEBUG='+debug,'BUILD_DIR=build-rk325-'+kind,'BIN_DIR=bin-rk325-'+kind,'MPIF90='+toolchain['wrapper_resolved'],'CC='+toolchain['c_compiler_resolved'],'FP_CONTRACT=off','STAGE202=on','STAGE203=on']
   row=dict(status='building',source_family='candidate',command=command,**stamp('started'))
   record['builds'][kind]=row;record['build_order'].append(kind);save(out/'builds.json',record)
   print('Building Stage325 '+kind+' with exact GCC13.2.0/MVAPICH4.0',flush=True)
   try:run_logged_build(command,tree,out/('build-'+kind+'.log'))
   except BaseException as exc:row.update(status='failed',error=repr(exc),returncode=getattr(exc,'returncode',None));raise
   finally:
    row.update(stamp('ended'));row['elapsed_seconds']=(row['ended_monotonic_ns']-row['started_monotonic_ns'])/1e9;save(out/'builds.json',record)
   binary=tree/('bin-rk325-'+kind)/'climate'
   for name,h in manifest['candidate'].items():
    if sha(tree/name)!=h:raise ValueError('Source changed during build: '+name)
   physics={p.name:sha(p) for p in (tree/'src/physics/simple_physics/phyparam/obj').glob('*.o')}
   if physics!=json.loads((kit/'reference-physics-objects.json').read_text()):raise ValueError('Physics objects differ from native job55551')
   libdir=out/('link-'+kind);libdir.mkdir()
   row.update(status='passed',binary=str(binary.relative_to(out)),binary_sha256=sha(binary),physics_objects=physics,libraries=linked(binary,libdir,toolchain['expected']),build_log_sha256=sha(out/('build-'+kind+'.log')))
   save(out/'builds.json',record)
  record['status']='passed'
 except BaseException as exc:record.update(status='failed',error=repr(exc));raise
 finally:
  record.update(stamp('ended'));record['elapsed_seconds']=(record['ended_monotonic_ns']-record['started_monotonic_ns'])/1e9;save(out/'builds.json',record)
 return record
