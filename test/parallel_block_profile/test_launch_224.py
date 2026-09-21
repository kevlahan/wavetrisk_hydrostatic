import contextlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import cluster_224 as c
import rank_launch_223 as r


class LaunchTests(unittest.TestCase):
    def test_whole_node_step(self):
        command=c.launch(Path('/kit'),Path('/placement'),40,'--probe')
        self.assertIn('--cpu-bind=none',command)
        self.assertIn('/kit/rank_launch_223.py',command)
        self.assertFalse(any(x.startswith(('--cpus-per-task','--exact','--exclusive','--hint')) for x in command))
        batch=Path(c.__file__).with_name('cluster_224.sbatch').read_text()
        self.assertNotIn('#SBATCH --cpus-per-task',batch)
        self.assertNotIn('#SBATCH --hint',batch)
        self.assertIn('#SBATCH --exclusive',batch)

    def test_inherited_restrictions(self):
        env={key:'1' for key in c.STEP_RESTRICTIONS}
        env.update(SLURM_JOB_ID='54724',SLURM_JOB_CPUS_PER_NODE='256')
        removed=c.clear_step_restrictions(env)
        self.assertEqual(set(removed),set(c.STEP_RESTRICTIONS))
        self.assertEqual(env,dict(SLURM_JOB_ID='54724',SLURM_JOB_CPUS_PER_NODE='256'))

    def test_failed_probe_preserves_diagnostic(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);error=io.StringIO()
            result=subprocess.CompletedProcess([],1,'partial rank output\n','Requested CPU unavailable\n')
            with patch.object(c.subprocess,'run',return_value=result),contextlib.redirect_stderr(error):
                with self.assertRaisesRegex(RuntimeError,'Placement probe for 40 ranks exited 1'):
                    c.probe_placement(Path('/kit'),out,40,{'mapping':[]})
            self.assertEqual((out/'placement-40-probe.stderr.log').read_text(),result.stderr)
            self.assertEqual((out/'placement-40-probe.stdout.log').read_text(),result.stdout)
            self.assertEqual(json.loads((out/'placement-40-probe.json').read_text())['returncode'],1)
            self.assertIn(result.stderr,error.getvalue())

    def test_affinity_error_is_actionable(self):
        with patch.object(r.os,'sched_getaffinity',return_value={0,1},create=True):
            with self.assertRaisesRegex(ValueError,r'Rank 3 .* requested CPU 128 unavailable; allowed CPUs=\[0, 1\]'):
                r.pin(3,128)

    def test_pin_before_solver(self):
        with patch.object(r.os,'sched_getaffinity',return_value={0,128},create=True), \
             patch.object(r.os,'sched_setaffinity',create=True) as pin,patch.dict(r.os.environ,{},clear=True):
            r.pin(3,128)
            pin.assert_called_once_with(0,{128})
            self.assertEqual(r.os.environ['MV2_ENABLE_AFFINITY'],'0')


if __name__=='__main__':unittest.main()
