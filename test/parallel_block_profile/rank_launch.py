#!/usr/bin/env python3
"""Record application-rank placement before exec; never samples the launcher."""
import json
import os
from pathlib import Path
import socket
import subprocess
import sys


def main():
    directory, *command = sys.argv[1:]
    if not command or "SLURM_PROCID" not in os.environ:
        raise SystemExit("Use under srun: rank_launch.py NEW_DIRECTORY executable [arguments]")
    rank = int(os.environ["SLURM_PROCID"])
    result = {"rank": rank, "hostname": socket.gethostname(), "pid_before_exec": os.getpid(),
              "affinity": sorted(os.sched_getaffinity(0)), "command": command,
              "slurm": {k: os.environ[k] for k in ("SLURM_JOB_ID", "SLURM_LOCALID", "SLURM_NODEID")
                        if k in os.environ}}
    result["cpu_model"] = next((line.split(":", 1)[1].strip() for line in
        Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")), "unavailable")
    result["memory_affinity"] = next((line for line in Path("/proc/self/status").read_text().splitlines()
                                     if line.startswith("Mems_allowed_list:")), "unavailable")
    with (Path(directory) / f"rank-{rank}.json").open("x") as output:
        json.dump(result, output, indent=2)
    # GNU time observes the actual application (or perf wrapper), not srun.
    # This adds no sampling/timers inside the numerical program.
    timer = "/usr/bin/time"
    version = subprocess.check_output([timer, "--version"], stderr=subprocess.STDOUT, text=True)
    if "GNU" not in version:
        raise SystemExit("GNU /usr/bin/time required for comparable per-rank RSS/fault reports")
    os.execv(timer, [timer, "-v", "-o", str(Path(directory) / f"rank-{rank}.resources"), *command])


if __name__ == "__main__":
    main()
