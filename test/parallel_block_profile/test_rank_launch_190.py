"""Run real child processes; verify exit propagation and per-child accounting."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from rank_metrics_190 import rank_resources


class RankLaunch190(unittest.TestCase):
    def launch(self, root, command, **kwargs):
        return subprocess.run([sys.executable,str(Path(__file__).with_name('rank_launch_190.py')),str(root),*command],
            env={**os.environ,'SLURM_PROCID':'0','PATH':'/nonexistent'},capture_output=True,text=True,**kwargs)

    def test_inherited_launcher_descriptor_reaches_child(self):
        read_fd,write_fd=os.pipe()
        try:
            os.write(write_fd,b'launcher descriptor')
            os.close(write_fd)
            with tempfile.TemporaryDirectory() as tmp:
                r=self.launch(Path(tmp),[sys.executable,'-c',
                    f'import os; print(os.read({read_fd},64).decode())'],pass_fds=(read_fd,))
                self.assertEqual(r.returncode,0,r.stderr)
                self.assertIn('launcher descriptor',r.stdout)
        finally:
            os.close(read_fd)

    def test_real_child_exit_and_resources_without_external_time(self):
        for code in (0,23):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                script='import sys; x=bytearray(8*1024*1024); print("solver ran"); sys.exit('+str(code)+')'
                r=self.launch(root,[sys.executable,'-c',script])
                self.assertEqual(r.returncode,code,r.stderr)
                self.assertIn('solver ran',r.stdout)
                record=rank_resources(root)[0]
                self.assertEqual(record['returncode'],code)
                self.assertEqual(record['resource_collector'],'os.wait4')
                self.assertGreater(record['max_rss_kib'],0)
                self.assertNotEqual(record['solver_pid'],record['launcher_pid'])

    def test_signal_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            r=self.launch(root,[sys.executable,'-c','import os,signal; os.kill(os.getpid(), signal.SIGTERM)'])
            self.assertEqual(r.returncode,128+signal.SIGTERM)
            self.assertEqual(rank_resources(root)[0]['returncode'],-signal.SIGTERM)

    def test_missing_executable_reports_original_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            r=self.launch(root,['/missing-stage190-executable'])
            self.assertEqual(r.returncode,127)
            self.assertIn('/missing-stage190-executable',r.stderr)
            record=json.loads((root/'rank-0.json').read_text())
            self.assertEqual(record['returncode'],127)
            self.assertIn('launch_error',record)

if __name__=='__main__': unittest.main()
