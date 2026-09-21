"""Recovery regression: purged Slurm jobs, interrupted scp and real queue errors."""
import os,subprocess,tempfile,unittest
from pathlib import Path

class Retrieval(unittest.TestCase):
    def fixture(self,root,queue=''):
        home=root/'home';home.mkdir();kit=home/'kit';kit.mkdir();bin=root/'bin';bin.mkdir()
        (home/'stage223-state.sh').write_text(f'kit="{kit}"\njob=54799\nout="{home}/results"\ntree="{home}/tree"\n')
        (kit/'collect_cluster_logs_223.py').write_text('''import os,sys
from pathlib import Path
h=Path(os.environ['HOME'])
with (h/'collector-calls').open('a') as f:f.write(' '.join(sys.argv[1:])+'\\n')
if '--verify' not in sys.argv:(h/'stage223-review-54799.tar.gz').write_bytes(b'archive')
''')
        for name,text in {'squeue':'#!/bin/bash\n'+queue+'\n','sha256sum':'#!/bin/bash\nexit 0\n','scp':'#!/bin/bash\nprintf "%s\\n" "$*" >> "$HOME/scp-calls"\n'}.items():
            p=bin/name;p.write_text(text);p.chmod(0o755)
        return home,dict(os.environ,HOME=str(home),USER='test',PATH=str(bin)+os.pathsep+os.environ['PATH'])
    def run_script(self,env):
        return subprocess.run(['bash',str(Path(__file__).with_name('collect_223.sh'))],env=env,capture_output=True,text=True)
    def test_purged_job_and_repeated_transfer(self):
        with tempfile.TemporaryDirectory() as d:
            home,env=self.fixture(Path(d),'exit 0')
            self.assertEqual(self.run_script(env).returncode,0)
            self.assertEqual(self.run_script(env).returncode,0)
            calls=(home/'collector-calls').read_text();self.assertEqual(calls.count('--stage'),1);self.assertEqual(calls.count('--verify'),2)
            self.assertEqual((home/'scp-calls').read_text().count('kevlahan@murier.math.mcmaster.ca:~/Desktop/'),2)
    def test_running_job_does_not_collect(self):
        with tempfile.TemporaryDirectory() as d:
            home,env=self.fixture(Path(d),'echo 54799')
            self.assertNotEqual(self.run_script(env).returncode,0);self.assertFalse((home/'scp-calls').exists())
    def test_queue_error_does_not_mean_completed(self):
        with tempfile.TemporaryDirectory() as d:
            home,env=self.fixture(Path(d),'exit 1')
            self.assertNotEqual(self.run_script(env).returncode,0);self.assertFalse((home/'collector-calls').exists())
if __name__=='__main__':unittest.main()
