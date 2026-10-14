#!/usr/bin/env python3
"""Summarize Stage 179 detail windows and compare printed states, not CPU fields.

No third-party dependencies. Input logs are read-only; JSON goes to stdout.
"""
import argparse
import json
import re
import statistics
from pathlib import Path

NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"
ROW = re.compile(r"^\s*(\d+) (.{32})\s+(.+)$")
STEP = re.compile(r"^\s*\d+\s+" + NUMBER + r"\s+d dt =.*\bcpu =\s*(" + NUMBER + r")")


def number(text):
    return float(text.replace("D", "E").replace("d", "e"))


def analyze(path):
    lines = path.read_text().splitlines()
    windows, states, ordinary, checkpoint, records = [], [], [], [], []
    pending_checkpoint = False
    pending_remap = False
    current = None
    for line in lines:
        if "Saving checkpoint" in line:
            pending_checkpoint = True
        if 'Remapping vertical coordinates' in line:
            pending_remap = True
        match = STEP.match(line)
        if match:
            records.append({"state": re.sub(r"\s+cpu =.*", "", line.strip()),
                            "seconds": number(match[1]), "checkpoint": pending_checkpoint,
                            "remap": pending_remap})
            (checkpoint if pending_checkpoint else ordinary).append(number(match[1]))
            pending_checkpoint = False
            pending_remap = False
            states.append(re.sub(r"\s+cpu =.*", "", line.strip()))
        elif line.strip().startswith(("Minimum relative mass =", "time [d] =")):
            states.append(line.strip())
        if line.startswith("Block detail profile:"):
            current = {"regions": [], "counters": {}, "storage": {}, "steps": {},
                       "boundary_callers": [], "rank_self_wall": {}, "dropped_steps": 0}
            windows.append(current)
            continue
        if current is None:
            continue
        stripped = line.strip()
        if stripped.startswith("detail schema = "):
            current["region_count"] = int(stripped.split("regions = ")[1])
            continue
        if stripped.startswith("detail step "):
            fields = stripped.split()
            key = fields[2] + ":" + fields[4]
            step = current["steps"].setdefault(key, {"sequence": int(fields[2]), "rank": int(fields[4])})
            if fields[5] == "restarts":
                step.update(restarts=int(fields[6]), remaps=int(fields[8]),
                            wall=number(fields[10]), cpu=number(fields[12]))
            else:
                step[fields[5]] = list(map(number, fields[7:]))
            continue
        if stripped.startswith("detail dropped step records = "):
            current["dropped_steps"] = int(stripped.split(" = ")[1])
            continue
        if stripped.startswith("detail boundary caller "):
            head, values = stripped.split(" = ")
            values = values.split()
            current["boundary_callers"].append({"caller_id": int(head.split()[3]),
                "wall_avg": number(values[0]), "cpu_avg": number(values[1]),
                "wall_max": number(values[2]), "calls": int(values[3])})
            continue
        if stripped.startswith("detail rank ") and " self-wall by region id = " in stripped:
            head, values = stripped.split(" = ")
            current["rank_self_wall"][int(head.split()[2])] = list(map(number, values.split()))
            continue
        match = ROW.match(line)
        if match:
            values = match[3].split()
            if len(values) != 8 or not match[2].strip()[0].isalpha():
                continue
            current["regions"].append(dict(zip(
                ("id", "region", "self_wall_avg", "self_wall_min", "self_wall_max", "self_cpu_avg",
                 "self_cpu_on_max_wall_rank", "inclusive_wall_avg", "max_wall_rank", "global_calls"),
                (int(match[1]), match[2].strip(), *map(number, values[:6]), int(values[6]), int(values[7])))))
        elif " global = " in line:
            key, value = line.strip().split(" global = ")
            current["counters"][key] = int(value)
        elif " sum/max = " in line:
            key, value = line.strip().split(" sum/max = ")
            current["storage"][key] = dict(zip(("sum_rank_peaks", "max_rank_peak"), map(int, value.split())))
    for window in windows:
        regions = window["regions"]
        window["sum_self_wall_avg"] = sum(r["self_wall_avg"] for r in regions)
        window["ranked_self_wall"] = [r["region"] for r in sorted(regions, key=lambda r: r["self_wall_avg"], reverse=True)]
        root = next((r for r in regions if r["id"] == 1), None)
        if root:
            window["timestep_inclusive_avg"] = root["inclusive_wall_avg"]
            # Nonzero residual can include instrumented work between timestep
            # scopes; printed timer rounding also contributes. Never rescale it.
            window["outside_or_rounding_seconds"] = window["sum_self_wall_avg"] - root["inclusive_wall_avg"]
        window["step_summaries"] = summarize_steps(window)
    total = [number(m[1]) for line in lines if (m := re.search(r"Total cpu time =\s*(" + NUMBER + ")", line))]
    failures = [line for line in lines if re.search(r"ERROR STOP|runtime error:|MPI_Abort|signal SIG", line)]
    return {"log": str(path.resolve()), "completed": bool(total) and not failures,
            "total_cpu_reported": total, "failures": failures,
            "non_checkpoint_step_times": ordinary, "checkpoint_step_times": checkpoint,
            "non_checkpoint_median": statistics.median(ordinary) if ordinary else None,
            "detail_windows": windows, "printed_states": states, "step_records": records}


def summarize_steps(window):
    """Rank-local attribution. Never sum independent phase maxima into a path."""
    groups = {}
    for step in window["steps"].values():
        groups.setdefault(step["sequence"], []).append(step)
    result = []
    expected_ranks = set(window["rank_self_wall"])
    for sequence, steps in sorted(groups.items()):
        complete = all(all(k in s for k in ("wall", "cpu", "self-wall", "self-CPU")) for s in steps)
        complete &= {s["rank"] for s in steps} == expected_ranks and bool(expected_ranks)
        count = window.get("region_count")
        if count is not None:
            complete &= all(len(s.get("self-wall", [])) == count and len(s.get("self-CPU", [])) == count for s in steps)
        summary = {"sequence": sequence, "complete": complete and not window["dropped_steps"]}
        if complete:
            slow = max(steps, key=lambda s: s["wall"])
            summary.update(wall_avg=statistics.mean(s["wall"] for s in steps),
                           wall_max=slow["wall"], max_wall_rank=slow["rank"],
                           restart=any(s["restarts"] for s in steps), remap=any(s["remaps"] for s in steps),
                           max_self_conservation_error=max(abs(sum(s["self-wall"])-s["wall"]) for s in steps),
                           self_on_max_wall_rank=slow["self-wall"],
                           self_avg=[statistics.mean(v) for v in zip(*(s["self-wall"] for s in steps))])
        result.append(summary)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path)
    args = parser.parse_args()
    reports = [analyze(path) for path in args.logs]
    reference = reports[0]["printed_states"]
    for report in reports:
        report["printed_states_match_first_log"] = bool(reference) and report["printed_states"] == reference
    print(json.dumps(reports, indent=2))


if __name__ == "__main__":
    main()
