import json
from pathlib import Path
import tempfile
import unittest
import contextlib
import io
import os
import sys

from analyze import analyze
from compare import summarize, resources
from experiment import schedule, validate_launch_mode, launch_command, execute_with_tee, SWITCHES
from samples import read_stat


class ExperimentTests(unittest.TestCase):
    def test_every_numerical_oracle_is_explicit(self):
        for name in ('DYNAMICS', 'ADAPTATION', 'REMAP'):
            self.assertIn('WAVETRISK_VALIDATE_BLOCK_'+name, SWITCHES)

    def test_launch_modes(self):
        validate_launch_mode(True, {})
        validate_launch_mode(False, {"SLURM_JOB_ID": "123"})
        with self.assertRaises(ValueError):
            validate_launch_mode(False, {})
        with self.assertRaises(ValueError):
            validate_launch_mode(True, {"SLURM_JOB_ID": "123"})

    def test_standalone_command(self):
        direct = launch_command(83, "/usr/bin/python3", Path("/helpers"), Path("/placement"),
                                ["./climateJ5", "simple.in"], "bb", direct=True)
        self.assertEqual(direct, ["srun", "-n", "83", "--partition=bb", "./climateJ5", "simple.in"])
        command = launch_command(83, "/usr/bin/python3", Path("/helpers"), Path("/placement"),
                                 ["/legacy/climate", "simple.in"])
        self.assertEqual(command[0:3], ["srun", "--ntasks", "83"])
        self.assertFalse(any(x.startswith("--partition") for x in command))
        self.assertFalse(any(x.startswith("--nodelist") for x in command))
        self.assertNotIn("salloc", command)
        self.assertNotIn("sbatch", command)
        self.assertEqual(command[-2:], ["/legacy/climate", "simple.in"])
        explicit = launch_command(83, "/usr/bin/python3", Path("/helpers"), Path("/placement"),
                                  ["/block/climate", "simple.in"], "bb", "bb02")
        self.assertIn("--partition=bb", explicit)
        self.assertIn("--nodelist=bb02", explicit)

    def test_tee_keeps_failure_status(self):
        console, log = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(console):
            code = execute_with_tee([sys.executable, "-c", "print('test output'); raise SystemExit(7)"],
                                    Path.cwd(), os.environ.copy(), log)
        self.assertEqual(code, 7)
        self.assertEqual(console.getvalue(), "test output\n")
        self.assertEqual(log.getvalue(), console.getvalue())

    def test_order_and_separate_sampling(self):
        self.assertEqual(schedule(2, 0, None), [("timing", 0, "legacy"), ("timing", 0, "block"),
                                               ("timing", 1, "block"), ("timing", 1, "legacy")])
        self.assertTrue(all(r[0] == "sample" for r in schedule(6, 2, "stat")))

    def test_perf_counters(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rank-0.stat"
            path.write_text("1000;;instructions;100;100.00\n500;;cycles;100;75.00\n"
                            "<not supported>;;cache-misses;0;0.00\n")
            result = read_stat(path)
            self.assertEqual(result["instructions_per_cycle"], 2)
            self.assertEqual(len(result["warnings"]), 2)

    def test_resource_scope_and_missing_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "rank-placement").mkdir()
            (root / "rank-placement/rank-0.resources").write_text(
                "\tMaximum resident set size (kbytes): 1024\n\tUser time (seconds): 1.2\n"
                "\tSystem time (seconds): 0.1\n\tMajor (requiring I/O) page faults: 0\n"
                "\tMinor (reclaiming a frame) page faults: 100\n\tInvoluntary context switches: 2\n")
            self.assertEqual(resources(root, 1)["max_rank_peak_rss_kib"], 1024)
            self.assertFalse(resources(root, 2)["available"])

    def fixture(self, root):
        runs = []
        for group, pair, label in schedule(3, 0, None):
            name = f"{pair}-{label}"
            directory = root / name
            (directory / "rank-placement").mkdir(parents=True)
            (directory / "rank-placement/rank-0.json").write_text(json.dumps({
                "rank": 0, "hostname": "node", "cpu_model": "test", "affinity": [0], "memory_affinity": "0"}))
            seconds = 2 if label == "block" else 1
            (directory / "run.log").write_text(
                f"00000322 0.3 d dt = 100 s Jmax = 7 dof = 10 cpu = {seconds}\n"
                "Saving checkpoint 4\n00000322 0.31 d dt = 100 s Jmax = 6 dof = 9 cpu = 20\n"
                f"Total cpu time = {seconds + 20}\n")
            runs.append(dict(group=group, pair=pair, label=label, directory=name, returncode=0))
        (root / "experiment.json").write_text(json.dumps(dict(ranks=1, runs=runs, fixture={})))

    def test_comparison_and_placement_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            report = summarize(root)
            self.assertTrue(report["valid_for_comparison"])
            timing = report["comparisons"]["timing"]
            self.assertEqual(timing["mean_ratio"], 2)
            self.assertEqual(timing["mean_excess_seconds_per_ordinary_step"], 1)
            self.assertAlmostEqual(timing["fraction_of_block_time_to_remove_for_20pct_target"], .4)
            path = root / "1-block/rank-placement/rank-0.json"
            data = json.loads(path.read_text())
            data["affinity"] = [1]
            path.write_text(json.dumps(data))
            report = summarize(root)
            self.assertFalse(report["valid_for_comparison"])
            self.assertFalse(report["comparisons"])

    def test_direct_measurements_are_provisional(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            path = root / "experiment.json"
            data = json.loads(path.read_text())
            data["launch_mode"] = "independent-srun"
            path.write_text(json.dumps(data))
            for p in root.glob("*/rank-placement/*.json"):
                p.unlink()
            report = summarize(root)
            self.assertTrue(report["measurements_usable"])
            self.assertFalse(report["valid_for_comparison"])
            self.assertFalse(report["placement_verified"])
            self.assertEqual(report["comparisons"]["timing"]["mean_ratio"], 2)

    def test_step_snapshots_and_boundary_callers(self):
        text = "\n".join([
            "Block detail profile: test", "  detail schema = 2; regions = 2",
            "  detail rank 0 self-wall by region id = .2 .8",
            "  detail boundary caller 2 dynamics residual wall-avg CPU-avg wall-max calls = .5 .4 .5 8",
            "  detail dropped step records = 0",
            "  detail step 7 rank 0 restarts 0 remaps 1 wall 1 CPU .9",
            "  detail step 7 rank 0 self-wall = .2 .8",
            "  detail step 7 rank 0 self-CPU = .2 .7",
            "Total cpu time = 1"])
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "detail.log"
            log.write_text(text)
            window = analyze(log)["detail_windows"][0]
            step = window["step_summaries"][0]
            self.assertTrue(step["complete"])
            self.assertTrue(step["remap"])
            self.assertEqual(step["max_self_conservation_error"], 0)
            self.assertEqual(window["boundary_callers"][0]["wall_avg"], .5)
            log.write_text(text.replace("records = 0", "records = 1"))
            self.assertFalse(analyze(log)["detail_windows"][0]["step_summaries"][0]["complete"])
            log.write_text(text.replace("self-CPU = .2 .7", "self-CPU = .2"))
            self.assertFalse(analyze(log)["detail_windows"][0]["step_summaries"][0]["complete"])


if __name__ == "__main__":
    unittest.main()
