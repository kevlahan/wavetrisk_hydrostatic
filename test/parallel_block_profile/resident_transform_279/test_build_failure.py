import contextlib,io,subprocess,sys,tempfile,unittest
from pathlib import Path
from build import run_logged_build

class BuildFailure(unittest.TestCase):
    def test_compiler_diagnostic_and_original_failure_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);log=root/'build.log';output=io.StringIO()
            cmd=[sys.executable,'-c',"import sys; print('source.f90: compiler diagnostic',file=sys.stderr); sys.exit(2)"]
            with contextlib.redirect_stdout(output),self.assertRaises(subprocess.CalledProcessError) as caught:
                run_logged_build(cmd,root,log)
            self.assertEqual(caught.exception.returncode,2)
            self.assertIn('source.f90: compiler diagnostic',output.getvalue())
            self.assertIn(str(log),output.getvalue())
            self.assertIn('source.f90: compiler diagnostic',log.read_text())
    def test_success_stays_in_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);log=root/'build.log';output=io.StringIO()
            with contextlib.redirect_stdout(output):
                run_logged_build([sys.executable,'-c',"print('successful build')"],root,log)
            self.assertEqual(output.getvalue(),'')
            self.assertIn('successful build',log.read_text())
