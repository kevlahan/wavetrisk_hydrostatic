import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("analyze", Path(__file__).with_name("analyze.py"))
analyzer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analyzer)


class LogTests(unittest.TestCase):
    def test_windows_states_and_checkpoint(self):
        row = lambda i, name, self_time, inclusive: (
            f"{i:3d} {name:32s} {self_time} 0.1 0.3 0.1 0.2 {inclusive} 2 83")
        text = "\n".join([
            "Minimum relative mass = 9.55E-01",
            "00000322 0.3021 d dt = 103.6 s Jmax = 7 dof = 62782 cpu = 1.71E+00",
            "Saving checkpoint 4 at time [day] = 0.31",
            "00000322 0.3104 d dt = 103.6 s Jmax = 6 dof = 64222 cpu = 6.41E+00",
            "Block detail profile: self times exclude nested instrumented regions; inclusive times overlap.",
            row(1, "timestep residual", 0.2, 0.5),
            row(5, "source extraction", 0.3, 0.3),
            "  source proof payload bytes global = 1000",
            "  expanded scalar records bytes sum/max = 4000 1200",
            "Total cpu time = 8.1200E+00"])
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "test.log"
            log.write_text(text)
            result = analyzer.analyze(log)
        self.assertTrue(result["completed"])
        self.assertEqual(result["non_checkpoint_step_times"], [1.71])
        self.assertEqual(result["checkpoint_step_times"], [6.41])
        window = result["detail_windows"][0]
        self.assertAlmostEqual(window["outside_or_rounding_seconds"], 0)
        self.assertEqual(window["ranked_self_wall"], ["source extraction", "timestep residual"])
        self.assertEqual(window["counters"]["source proof payload bytes"], 1000)
        self.assertEqual(window["regions"][0]["max_wall_rank"], 2)
        self.assertNotIn("cpu =", result["printed_states"][1])

    def test_failure_not_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "bad.log"
            log.write_text("ERROR STOP test\nTotal cpu time = 1.0\n")
            self.assertFalse(analyzer.analyze(log)["completed"])


if __name__ == "__main__":
    unittest.main()
