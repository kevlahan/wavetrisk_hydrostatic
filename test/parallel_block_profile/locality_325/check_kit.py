"""Static fail-closed seal and scope gate; never compiles or launches a solver."""
from pathlib import Path
import hashlib,json,sys,tarfile,posixpath
from common import verify,sha
from config325 import SCHEDULE
REQUIRED=('sources.tar.gz','source-files.json','expected-toolchain.json','simple.in',
 'verify325.py','account325.py','qualification.sbatch','run325.py','build.py',
 'continuation_binding.py','replay-binding.json','rk322-semantic-review-55567.tar.gz',
 'current256.json','owner_label_permutation256.json')
MAP_HASHES={'current256.json':'d5519da50e82252d5f7dede09bc99633056e29e982770af559f0c836597ed75a','owner_label_permutation256.json':'0ec557fc4220e6c1b7c09fe53c8af25602403bbd78fd4c7280e3fbde540ced36'}

def check(kit):
 kit=Path(kit).resolve(strict=True);manifest=verify(kit)
 for name in REQUIRED:
  if name not in manifest:raise ValueError('Missing gate artifact: '+name)
 for name,want in MAP_HASHES.items():
  if sha(kit/name)!=want:raise ValueError('Audited exact map differs: '+name)
 if tuple((x['label'],x['kind'],x['owner'],x['ranks'],x['steps']) for x in SCHEDULE)!=(('checked-L256','checked','L256',256,4),('qual-C0','candidate','off',256,4),('qual-L256','candidate','L256',256,4)):raise ValueError('Sealed three-run schedule differs')
 want=json.loads((kit/'source-files.json').read_text())
 if set(want)!={'candidate'} or not isinstance(want['candidate'],dict) or len(want['candidate'])<100:raise ValueError('Candidate source manifest invalid')
 found={}
 with tarfile.open(kit/'sources.tar.gz') as t:
  for m in t:
   p=Path(m.name)
   if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='candidate':raise ValueError('Unsafe source path')
   if m.isdir():continue
   name=str(Path(*p.parts[1:]))
   if name in found:raise ValueError('Duplicate source file')
   if m.isfile():found[name]=hashlib.file_digest(t.extractfile(m),'sha256').hexdigest()
   elif m.issym():
    target=posixpath.normpath(posixpath.join(str(p.parent),m.linkname))
    if Path(m.linkname).is_absolute() or not target.startswith('candidate/'):raise ValueError('Unsafe source link')
    found[name]=hashlib.sha256(m.linkname.encode()).hexdigest()
   else:raise ValueError('Unexpected source member type')
 if found!=want['candidate']:raise ValueError('Packaged source identity differs')
 fields={}
 for line in (kit/'simple.in').read_text().splitlines():
  f=line.split('!',1)[0].split()
  if f:
   if len(f)!=2 or f[0] in fields:raise ValueError('Ambiguous solver input')
   fields[f[0]]=f[1]
 for key,value in {'max_level':'5','tol':'0','zlevels':'30','resume':'3','run_id':'test','physics_type':'Simple','NCAR_topo':'F','sso':'F'}.items():
  if fields.get(key)!=value:raise ValueError('Selected fixed configuration differs: '+key)
 expected=json.loads((kit/'expected-toolchain.json').read_text())
 if expected.get('gcc_version')!='13.2.0' or expected.get('mpi_version')!='4.0':raise ValueError('Exact native compiler/MPI version gate absent')
 from continuation_binding import load_qualification
 report,qualified,native,binding=load_qualification(kit)
 if report['passed'] is not True or qualified['passed'] is not True:raise ValueError('Prior qualification evidence absent')
 from verify325 import capture_environment
 for spec in SCHEDULE:
  env=capture_environment(spec['kind'],spec['owner'])
  if env.get('WAVETRISK_EXECUTION_AUDIT316')!='1' or env.get('WAVETRISK_COMPACT_AUDIT325')!='1':raise ValueError('Compact mode must be explicit')
 return dict(stage=325,revision='locality-qualification-v1',passed=True,files=len(manifest),source_files=len(found),package_manifest_sha256=sha(kit/'files-sha256.json'),source_manifest_sha256=sha(kit/'source-files.json'),builds=2,solver_launches=3,performance_launches=0,prior_qualification_job=55567)
if __name__=='__main__':
 result=check(Path(__file__).resolve().parent);print('PASS sealed Stage325 compact locality qualification:',result['files'],'files;2builds3solverlaunches0timing')
