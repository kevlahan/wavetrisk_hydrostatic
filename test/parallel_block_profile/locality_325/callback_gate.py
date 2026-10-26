"""Source-bound callback runtime gate, before expensive solver builds."""
from pathlib import Path
import subprocess
from common import sha,save

def check(kit,source,out,compiler='mpif90',native=True):
 kit,source,out=map(Path,(kit,source,out));text=(source/'src/wavelet.f90').read_text();test=(kit/'callback_check.f90').read_text()
 a=text.index('  subroutine subdivision_inverse_oracle(');b=text.index('  subroutine inverse_wavelet_transform',a)
 chunk=text[a:b]
 a=text.index('  logical :: subdivision_oracle_active');b=text.index('\n  \n  interface',a);state=text[a:b]
 if chunk not in test or state not in test or '\n  contains\n' in chunk:raise ValueError('Callback harness/source drift')
 work=out/'callback-gate';work.mkdir();binary=work/'callback-check'
 command=[compiler,'-O0','-g','-fcheck=all','-Wall','-Wextra','-Wno-unused-dummy-argument','-Wno-compare-reals','-Werror','-Werror=trampolines','-ffree-line-length-none']
 if native:command+=['-Wl,-z,noexecstack']
 else:command+=['-ftrampoline-impl=stack']
 command +=[str((kit/'callback_check.f90').resolve()),'-o',str(binary.resolve())]
 result=dict(status='building',source_sha256=sha(source/'src/wavelet.f90'),test_sha256=sha(kit/'callback_check.f90'),command=command,native=native,cases=[])
 save(out/'callback-gate.json',result)
 from build import run_logged_build
 run_logged_build(command,work,work/'build.log')
 if native:
  headers=subprocess.check_output(['readelf','-W','-l',str(binary.resolve())],text=True)
  stack=[l.split() for l in headers.splitlines() if 'GNU_STACK' in l]
  if len(stack)!=1 or any('E' in v or 'X' in v for v in stack[0][6:-1]):raise ValueError('Callback executable has executable stack')
  result['elf_stack_header']=' '.join(stack[0]);save(out/'callback-gate.json',result)
 for nd in (1,2,4):
  for nk in (1,3):
   for frozen in (0,1):
    p=subprocess.run([str(binary.resolve()),str(nd),str(nk),str(frozen),'0'],capture_output=True,text=True)
    result['cases'].append(dict(domains=nd,layers=nk,frozen=frozen,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr));save(out/'callback-gate.json',result)
    if p.returncode or 'PASS' not in p.stdout:raise ValueError('Native oracle callback failed; see callback-gate.json')
 for probe in range(1,6):
  p=subprocess.run([str(binary.resolve()),'1','1','0',str(probe)],capture_output=True,text=True)
  result['cases'].append(dict(invalid_context=probe,returncode=p.returncode,stderr=p.stderr));save(out/'callback-gate.json',result)
  if p.returncode==0 or 'ERROR STOP' not in p.stderr:raise ValueError('Callback context guard failed')
 result.update(status='passed',binary_sha256=sha(binary),scope='Exact oracle routine with deterministic small walker/numerical stubs. Does not replace full solver science/oracle validation.')
 save(out/'callback-gate.json',result);print('PASS source-bound oracle callback runtime gate',flush=True)
 return result
