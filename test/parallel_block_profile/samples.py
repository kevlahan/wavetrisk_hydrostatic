#!/usr/bin/env python3
"""External perf evidence from both unchanged executables, never a job-speedup estimate."""
import argparse
import json
from pathlib import Path
import subprocess


def read_stat(path):
    events, warnings = {}, []
    for line in path.read_text().splitlines():
        fields = line.split(";")
        if len(fields) < 3 or line.startswith("#"):
            continue
        value, unit, event = (s.strip() for s in fields[:3])
        if event not in ("task-clock", "cycles", "instructions", "cache-references", "cache-misses"):
            continue
        try:
            events[event] = {"value": float(value), "unit": unit}
            if len(fields) > 4 and fields[4].strip():
                running = float(fields[4].strip().rstrip("%"))
                events[event]["running_percent"] = running
                if running < 99:
                    warnings.append(f"{event}: multiplexed ({running}% running)")
        except ValueError:
            warnings.append(f"{event}: unavailable or unparsable: {line}")
    if not events:
        warnings.append("No supported counter records found")
    ipc = None
    if all(e in events for e in ("instructions", "cycles")) and events["cycles"]["value"] > 0:
        ipc = events["instructions"]["value"] / events["cycles"]["value"]
    return {"events": events, "instructions_per_cycle": ipc, "warnings": warnings}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", type=Path)
    args = parser.parse_args()
    root = args.experiment.resolve()
    manifest = json.loads((root / "experiment.json").read_text())
    result = {}
    for run in manifest["runs"]:
        directory = root / run["directory"] / "samples"
        counters = {}
        for path in sorted(directory.glob("rank-*.stat")):
            counters[path.stem] = read_stat(path)
        for path in sorted(directory.glob("rank-*.record")):
            output = path.with_suffix(".report.txt")
            with output.open("x") as stream:
                subprocess.run(["perf", "report", "--stdio", "--no-children", "--percent-limit", "0.5",
                                "--sort", "comm,dso,symbol", "-i", str(path)], stdout=stream, check=True)
            # Call-chain view identifies owners of shared allocation/copy/MPI helpers.
            with path.with_suffix(".callgraph.txt").open("x") as stream:
                subprocess.run(["perf", "report", "--stdio", "--children", "--percent-limit", "1",
                                "-i", str(path)], stdout=stream, check=True)
        result[run["label"]] = counters
    print(json.dumps({"samples": result, "scope": "Whole application, including startup/output and MPI spin. "
          "Use symbols/call paths to explain costs, not whole-run percentages as ordinary-step shares. "
          "Do not interpret generic cache misses as DRAM bytes or all extra instructions as block overhead."}, indent=2))


if __name__ == "__main__":
    main()
