#!/usr/bin/env python3
"""Regression gates for the previous missing-source/login-shell failure."""
import hashlib,io,json,os,subprocess,tarfile,tempfile,unittest
from pathlib import Path
SCRIPT=Path(__file__).with_name('bootstrap_215.sh')

class Bootstrap(unittest.TestCase):
    def setup_home(self,root):
        home=root/'home with spaces';home.mkdir();kit=home/'wavetrisk-stage215-cluster-kit';kit.mkdir()
        (kit/'bootstrap_215.sh').write_bytes(SCRIPT.read_bytes())
        data=b'complete frozen source\n'
        with tarfile.open(kit/'stage196-source.tar.gz','w:gz') as t:
            i=tarfile.TarInfo('src/example.f90');i.size=len(data);t.addfile(i,io.BytesIO(data))
        (kit/'stage215-sources.json').write_text(json.dumps(dict(baseline_commit='9bc3a19556dab79b3cba3944cf5df53517389da3',baseline_sources={'src/example.f90':hashlib.sha256(data).hexdigest()})))
        bin=root/'bin';bin.mkdir();stub=bin/'sbatch'
        stub.write_text('#!/bin/bash\nprintf "submitted\\n" >> "$HOME/submissions"\nprintf "12345\\n"\n');stub.chmod(0o755)
        env=dict(os.environ,HOME=str(home),PATH=str(bin)+os.pathsep+os.environ['PATH'])
        return home,kit,env

    def test_fresh_source_then_single_submission_and_safe_retry(self):
        with tempfile.TemporaryDirectory() as d:
            home,kit,env=self.setup_home(Path(d));state=home/'stage215-state.sh';state.write_text('old state\n')
            r=subprocess.run(['bash',str(kit/'bootstrap_215.sh')],env=env,capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertTrue((home/'wt-stage215/src/example.f90').is_file())
            self.assertEqual(list(home.glob('stage215-state.previous.*'))[0].read_text(),'old state\n')
            loaded=subprocess.check_output(['bash','-c','source "$HOME/stage215-state.sh"; printf "%s\\n" "$build_job" "$tree"'],env=env,text=True)
            self.assertEqual(loaded.splitlines(),['12345',str(home/'wt-stage215')])
            before=state.read_bytes()
            retry=subprocess.run(['bash',str(kit/'bootstrap_215.sh')],env=env,capture_output=True)
            self.assertNotEqual(retry.returncode,0);self.assertEqual((home/'submissions').read_text(),'submitted\n')
            self.assertEqual(state.read_bytes(),before)

    def test_bad_source_never_submits_or_creates_tree(self):
        with tempfile.TemporaryDirectory() as d:
            home,kit,env=self.setup_home(Path(d));p=kit/'stage215-sources.json';m=json.loads(p.read_text());m['baseline_sources']['src/example.f90']='bad';p.write_text(json.dumps(m))
            r=subprocess.run(['bash',str(kit/'bootstrap_215.sh')],env=env,capture_output=True)
            self.assertNotEqual(r.returncode,0);self.assertFalse((home/'submissions').exists())
            self.assertFalse((home/'wt-stage215').exists());self.assertFalse((home/'stage215-state.sh').exists())

if __name__=='__main__':unittest.main()
