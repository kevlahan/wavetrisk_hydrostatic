"""Exercise the real frozen patch/build/submit sequence with a compiler stand-in.

Enable with STAGE204_TEST_KIT=/path/to/frozen/kit. Git patches and source hashes
are real; only compilation, scheduler topology and solver launch are replaced.
"""
import json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import build_restriction_204 as build
import submit_restriction_204 as submit
ROOT=Path(__file__).resolve().parents[2]
@unittest.skipUnless(os.environ.get('STAGE204_TEST_KIT'),'requires frozen Stage204 kit')
class Build204(unittest.TestCase):
 def test_actual_patch_order_flags_observers_and_submission(self):
  kit=Path(os.environ['STAGE204_TEST_KIT']).resolve()
  manifest=json.loads((kit/'stage204-sources.json').read_text())
  real_run=subprocess.run
  with tempfile.TemporaryDirectory() as tmp:
   tree=Path(tmp)/'tree';tree.mkdir()
   import io,tarfile
   with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(ROOT),'archive',manifest['baseline_commit']]))) as t:t.extractall(tree,filter='data')
   commands=[]
   def run(command,*args,**kwargs):
    if command[0]=='make':
     commands.append(command);cwd=Path(kwargs['cwd'])
     opts=dict(v.split('=',1) for v in command if '=' in v)
     dest=cwd/opts['BIN_DIR'];dest.mkdir();(dest/'climate').write_text(json.dumps(command));(dest/'climate').chmod(0o755)
     obj=cwd/'src/physics/simple_physics/phyparam/obj';obj.mkdir(exist_ok=True)
     for n in range(23):(obj/f'common{n}.o').write_bytes(bytes([n]))
     return subprocess.CompletedProcess(command,0)
    if command[0]=='lscpu':return subprocess.CompletedProcess(command,0,stdout='test topology')
    if len(command)>1 and command[1]==str(kit/'cluster_204.py'):
     commands.append(command);return subprocess.CompletedProcess(command,0)
    return real_run(command,*args,**kwargs)
   with patch.object(subprocess,'run',run):
    build.build(kit,tree)
    self.assertEqual(len(commands),7)
    for cmd in commands[:5]:
     self.assertIn('STAGE202=on',cmd);self.assertIn('FP_CONTRACT=off',cmd)
    self.assertIn('STAGE203=off',commands[0]);self.assertIn('STAGE203=on',commands[2])
    fixture=Path(tmp)/'fixture';fixture.mkdir()
    submit.submit(kit,tree,fixture,Path(tmp)/'out')
    self.assertIn('--legacy_profile',commands[-1])
    record=json.loads((tree/'restriction-204-builds.json').read_text())
    binary=tree/record['builds']['candidate']['binary'];binary.write_text('changed')
    with self.assertRaisesRegex(ValueError,'Executable changed'):submit.submit(kit,tree,fixture,Path(tmp)/'out2')
if __name__=='__main__':unittest.main()
