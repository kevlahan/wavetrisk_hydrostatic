"""Exact installed compiler/MPI gate, including runtime library identity."""
from pathlib import Path
import json,os,re,shlex,shutil,subprocess,platform
from common import save,sha

def python_dependencies():
 """Exercise the exact Python process used by the native runner/verifier."""
 import gzip,sys,zlib
 import numpy as np
 payload=b'Stage322 compressed qualification dependency gate\x00'*17
 compressed=gzip.compress(payload,mtime=0)
 if gzip.decompress(compressed)!=payload:raise ValueError('Python gzip/zlib roundtrip failed')
 values=np.frombuffer(bytes.fromhex('000000000000f03f0000000000000040'),dtype='<f8')
 if values.tolist()!=[1.,2.] or not bool(np.isfinite(values).all()):raise ValueError('Python numpy binary/finite check failed')
 return dict(passed=True,executable=str(Path(sys.executable).resolve()),python_version=platform.python_version(),
             numpy_version=np.__version__,zlib_compile_version=zlib.ZLIB_VERSION,
             zlib_runtime_version=zlib.ZLIB_RUNTIME_VERSION,gzip_roundtrip=True,numpy_binary_finite=True)

def require_version(text):
 if text.strip()!='13.2.0':raise ValueError('Require exact GCC 13.2.0; got '+repr(text.strip()))

def require_mpi_version(text):
 matches=re.findall(r'(?im)^\s*MVAPICH(?:2)?(?:\s+Version)?\s*:\s*([0-9]+(?:\.[0-9]+)+)\b',text)
 if matches!=['4.0']:raise ValueError('Require reported MVAPICH 4.0 library identity; full output retained')

def require_under(path,root):
 p=Path(path).resolve(strict=True);r=Path(root).resolve(strict=True)
 if not p.is_relative_to(r):raise ValueError('Unexpected toolchain path: '+str(p))
 return str(p)

def inspect_c_compiler(compiler):
 """Bind the helper to the exact Fortran compiler installation, allowing local aliases."""
 compiler=Path(compiler).resolve(strict=True)
 requested=compiler.with_name('gcc')
 resolved=Path(require_under(requested,compiler.parent.parent))
 version=subprocess.check_output([str(resolved),'-dumpfullversion'],text=True)
 require_version(version)
 target=subprocess.check_output([str(resolved),'-dumpmachine'],text=True).strip()
 if target!='x86_64-pc-linux-gnu':raise ValueError('Unexpected native C compiler target')
 return dict(c_compiler_requested=str(requested),c_compiler_resolved=str(resolved),
             c_compiler_sha256=sha(resolved),c_compiler_version=version,c_compiler_target=target)

def linked(binary,out,expected):
 p=subprocess.run(['ldd',str(binary)],text=True,capture_output=True)
 (out/(binary.name+'-ldd.log')).write_text(p.stdout+p.stderr)
 if p.returncode or 'not found' in p.stdout:raise ValueError('Unresolved dynamic libraries')
 libs={m[1]:m[2] for m in re.finditer(r'^\s*(\S+)\s+=>\s+(\S+)',p.stdout,re.M)}
 for prefix,root in [('libmpi.so',expected['mpi_root']),('libmpifort.so',expected['mpi_root']),('libgfortran.so',str(Path(expected['gcc_compiler']).parents[1]))]:
  values=[v for k,v in libs.items() if k.startswith(prefix)]
  if len(values)!=1:raise ValueError('Missing/unexpected linked '+prefix)
  require_under(values[0],root)
 return {k:dict(path=str(Path(v).resolve()),sha256=sha(Path(v).resolve())) for k,v in libs.items() if v.startswith('/')}

def inspect(kit,out):
 """Read-only compiler/wrapper identity check; no compile, MPI run or submission."""
 kit,out=Path(kit),Path(out);expected=json.loads((kit/'expected-toolchain.json').read_text())
 record=dict(status='checking_identity',expected=expected,host=platform.node())
 save(out/'toolchain313.json',record)
 try:
  record['phase']='platform_and_dependencies'
  if platform.system()!='Linux':raise ValueError('Native Linux toolchain gate required')
  record['python_dependencies']=python_dependencies()
  record['modules']=os.environ.get('LOADEDMODULES','')
  record['effective_dependencies']={k:os.environ.get(k) for k in expected['build_dependencies']}
  if record['effective_dependencies']!=expected['build_dependencies']:raise ValueError('Build environment was not normalized')
  for key in ('NETLIB_LAPACK_ROOT','NF_CONFIG','NC_CONFIG'):
   Path(expected['build_dependencies'][key]).resolve(strict=True)
  wrapper=shutil.which('mpif90')
  if not wrapper:raise ValueError('Missing mpif90')
  record['phase']='wrapper'
  record['wrapper_resolved']=require_under(wrapper,expected['mpi_root'])
  record['wrapper_sha256']=sha(Path(wrapper).resolve())
  command=subprocess.check_output([wrapper,'-show'],text=True,stderr=subprocess.STDOUT)
  record['wrapper_command']=command;tokens=shlex.split(command)
  compiler=Path(tokens[0]).resolve(strict=True)
  record['compiler_resolved']=str(compiler)
  if compiler!=Path(expected['gcc_compiler']).resolve(strict=True):raise ValueError('Wrapper selected unexpected compiler')
  record['compiler_sha256']=sha(compiler)
  record.update(inspect_c_compiler(compiler))
  record['phase']='compiler'
  record['compiler_version']=subprocess.check_output([str(compiler),'-dumpfullversion'],text=True)
  record['compiler_banner']=subprocess.check_output([wrapper,'--version'],text=True)
  record['target']=subprocess.check_output([str(compiler),'-dumpmachine'],text=True).strip()
  require_version(record['compiler_version'])
  if not re.fullmatch(r'GNU Fortran \(Spack GCC\) 13\.2\.0',record['compiler_banner'].splitlines()[0]):raise ValueError('Unexpected compiler banner')
  if record['target']!='x86_64-pc-linux-gnu' or record['c_compiler_target']!=record['target']:raise ValueError('Unexpected compiler target')
  if not set(expected['required_modules'])<=set(record['modules'].split(':')):raise ValueError('Established GCC/MPI/Fortran dependency modules missing')
  if '-lmpifort' not in tokens or '-lmpi' not in tokens:raise ValueError('Missing MPI Fortran linkage')
  for token in tokens:
   if token.startswith(('-I','-L')) and 'mvapich' in token.lower():require_under(token[2:],expected['mpi_root'])
  record['status']='identity_passed_runtime_pending'
 except Exception as exc:
  record.update(status='failed',error=repr(exc));raise
 finally:save(out/'toolchain313.json',record)
 return record

def check(kit,out):
 kit,out=Path(kit),Path(out)
 record=inspect(kit,out);expected=record['expected'];wrapper=record['wrapper_resolved']
 record.update(status='checking_runtime',phase='mpi_probe');save(out/'toolchain313.json',record)
 from build import run_logged_build
 binary=out/'mpi-identity'
 cmd=[wrapper,'-O0','-g','-Wall','-Wextra','-Werror','-fcheck=all',str(kit/'mpi_identity.f90'),'-o',str(binary)]
 record['probe_compile_command']=cmd;save(out/'toolchain313.json',record)
 run_logged_build(cmd,out,out/'mpi-identity-build.log')
 record['probe_libraries']=linked(binary,out,expected)
 cmd=['srun','--nodes=1','--ntasks=1','--cpu-bind=none','--kill-on-bad-exit=1',str(binary)]
 p=subprocess.run(cmd,text=True,capture_output=True,timeout=60)
 (out/'mpi-identity-run.log').write_text(p.stdout+p.stderr)
 record.update(probe_command=cmd,probe_returncode=p.returncode,library_version=p.stdout);save(out/'toolchain313.json',record)
 if p.returncode:raise ValueError('Native MPI identity probe failed')
 require_mpi_version(p.stdout)
 record['status']='passed';save(out/'toolchain313.json',record)
 return record
