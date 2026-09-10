#!/usr/bin/env python3
"""Compare paired runs, not overlapping timer sums or independently slow ranks.

Outputs measured excess over legacy and an explicitly hypothetical removal
scenario. No observed excess is labelled irreducible block overhead.
"""
import argparse
import json
import math
from pathlib import Path
import random
import statistics

from analyze import analyze


def interval(values):
    """Descriptive paired-run bootstrap; not independent timestep resampling."""
    rng = random.Random(179)
    means = sorted(statistics.mean(rng.choices(values, k=len(values))) for _ in range(4000))
    return [means[100], means[3899]]


def placement(directory, ranks):
    result = []
    for rank in range(ranks):
        item = json.loads((directory / "rank-placement" / f"rank-{rank}.json").read_text())
        if item["rank"] != rank or not item["affinity"]:
            raise ValueError("Invalid rank placement")
        result.append((item["hostname"], item["cpu_model"], item["affinity"], item["memory_affinity"]))
    return result


def resources(directory, ranks):
    keys = {"Maximum resident set size (kbytes)": "peak_rss_kib",
            "User time (seconds)": "user_seconds", "System time (seconds)": "system_seconds",
            "Major (requiring I/O) page faults": "major_faults",
            "Minor (reclaiming a frame) page faults": "minor_faults",
            "Involuntary context switches": "involuntary_switches"}
    values = []
    for rank in range(ranks):
        path = directory / "rank-placement" / f"rank-{rank}.resources"
        if not path.exists():
            return {"available": False}
        row = {}
        for line in path.read_text().splitlines():
            if ": " in line:
                key, value = line.strip().rsplit(": ", 1)
                if key in keys:
                    row[keys[key]] = float(value)
        if set(row) != set(keys.values()):
            return {"available": False, "reason": "Incomplete GNU time report"}
        values.append(row)
    return {"available": True, "sum_rank_peaks_rss_kib": sum(r["peak_rss_kib"] for r in values),
            "max_rank_peak_rss_kib": max(r["peak_rss_kib"] for r in values),
            "whole_run_rank_sums": {k: sum(r[k] for r in values) for k in values[0] if k != "peak_rss_kib"},
            "scope": "Whole run, including initialization/output and MPI spin; RSS peaks need not be simultaneous"}


def summarize(root, removals=()):
    manifest = json.loads((root / "experiment.json").read_text())
    groups, warnings, failures, detailed, memory_reports = {}, [], [], [], []
    direct = manifest.get("launch_mode") == "independent-srun"
    if direct:
        warnings.append("Direct srun runs have no rank wrapper: node/core placement and per-rank RSS are unverified. "
                        "Timing ratios are provisional, not a validated same-hardware comparison. "
                        "Separate allocations may encounter different nodes, cores and background load.")
    reference_placement = None
    reference_states = None
    for run in manifest["runs"]:
        directory = root / run["directory"]
        report = analyze(directory / "run.log")
        usage = resources(directory, manifest["ranks"])
        memory_reports.append({"run": directory.name, **usage})
        if not usage["available"]:
            warnings.append(f"No complete RSS/fault report: {directory.name}")
        if not direct:
            identity = placement(directory, manifest["ranks"])
            if reference_placement is None:
                reference_placement = identity
            if identity != reference_placement:
                failures.append(f"Rank/node/CPU affinity or allowed memory nodes differ: {directory.name}")
        states = [" ".join(s.split()) for s in report["printed_states"]]
        if reference_states is None:
            reference_states = states
        if states != reference_states or not states:
            failures.append(f"Printed states differ or missing: {directory.name}")
        if run["returncode"] or not report["completed"]:
            failures.append(f"Run incomplete: {directory.name}")
        if not report["non_checkpoint_step_times"]:
            failures.append(f"No ordinary timestep records: {directory.name}")
            continue
        groups.setdefault((run["group"], run["pair"]), {})[run["label"]] = report
        if run["label"] == "detail":
            detailed.append(report)
    comparisons = {}
    for group, numerator, denominator in (("timing", "block", "legacy"), ("overhead", "detail", "block")):
        ratios, seconds, legacy_means = [], [], []
        for (name, pair), reports in groups.items():
            if name != group:
                continue
            if numerator not in reports or denominator not in reports:
                failures.append(f"Incomplete {group} pair {pair}")
                continue
            a, b = reports[numerator], reports[denominator]
            keys = lambda r: [(" ".join(s["state"].split()), s["checkpoint"]) for s in r["step_records"]]
            if keys(a) != keys(b):
                failures.append(f"Different timestep/checkpoint workloads in {group} pair {pair}")
                continue
            x, y = [statistics.mean(r["non_checkpoint_step_times"]) for r in (a, b)]
            if not all(math.isfinite(v) and v > 0 for v in (x, y)):
                failures.append(f"Nonpositive/nonfinite times in {group} pair {pair}")
                continue
            ratios.append(x / y)
            seconds.append(x - y)
            legacy_means.append(y)
        if not ratios:
            continue
        result = {"paired_run_ratios": ratios, "mean_ratio": statistics.mean(ratios),
                  "ratio_range": [min(ratios), max(ratios)],
                  "mean_excess_seconds_per_ordinary_step": statistics.mean(seconds)}
        if len(ratios) >= 3:
            result["paired_run_bootstrap_95pct_descriptive"] = interval(ratios)
        else:
            warnings.append(f"Only {len(ratios)} {group} pairs; insufficient to characterize variability")
        if group == "timing":
            result["fraction_of_block_time_to_remove_for_20pct_target"] = max(0, 1 - 1.2 / result["mean_ratio"])
            result["mean_legacy_seconds_per_ordinary_step"] = statistics.mean(legacy_means)
        comparisons[group] = result
    scenarios = []
    for report in detailed:
        checkpoint_sequences = {i for i, row in enumerate(report["step_records"], 1) if row["checkpoint"]}
        steps = [s for w in report["detail_windows"] for s in w["step_summaries"]
                 if s["complete"] and not s["restart"] and s["sequence"] not in checkpoint_sequences]
        if not steps:
            failures.append(f"Missing complete ordinary per-step rank records: {report['log']}")
            continue
        if any(w["dropped_steps"] for w in report["detail_windows"]):
            failures.append(f"Dropped timestep detail records: {report['log']}")
        regions = {r["id"]: r["region"] for w in report["detail_windows"] for r in w["regions"]}
        walls = [s["wall_avg"] for s in steps]
        if max(s["max_self_conservation_error"] for s in steps) > 1e-3:
            failures.append(f"Per-step self-time conservation failed: {report['log']}")
        ranking = sorted(((i + 1, statistics.mean(s["self_avg"][i] for s in steps))
                          for i in range(len(steps[0]["self_avg"]))), key=lambda x: x[1], reverse=True)
        entry = {"log": report["log"], "ordinary_detail_steps": len(steps),
                 "rank_average_step_seconds": statistics.mean(walls),
                 "ranked_self_seconds": [{"id": i, "region": regions.get(i, "unused"), "seconds": seconds}
                                         for i, seconds in ranking if seconds > 0]}
        if removals:
            if any(i not in regions for i, _ in removals):
                raise ValueError("Removal region not measured in detail log")
            saved = statistics.mean(sum(s["self_avg"][i - 1] * f for i, f in removals) for s in steps)
            original = statistics.mean(walls)
            entry["hypothetical_rank_average_scenario"] = {
                "removed_self_seconds": saved, "retained_self_seconds": original - saved,
                "ratio_before_over_after": original / (original - saved) if saved < original else None,
                "assumption": "Selected self costs vanish by specified fractions with no replacement cost; "
                              "MPI waits unchanged. Not a predicted job speedup or irreducible overhead bound."}
        scenarios.append(entry)
    if not comparisons:
        failures.append("No paired timing/overhead observations; sample runs are not timing comparisons")
    if any(not v.get("contents_hashed", True) for v in manifest["fixture"].values()):
        warnings.append("Directory-linked assets were identified by path, not content-hashed; keep grids immutable")
    warnings.append("Verify archived source identity, compiler/ISA flags and library versions in both build logs. "
                    "Hashes and matching placement alone do not prove equivalent builds.")
    warnings.append("Printed-state agreement is not bitwise checkpoint equality. MPI CPU includes spinning. "
                    "Observed block-minus-legacy time is excess cost, not necessarily inherent block overhead.")
    return {"valid_for_comparison": not failures and not direct, "measurements_usable": not failures,
            "placement_verified": not failures and not direct, "failures": failures, "warnings": warnings,
            "comparisons": comparisons if not failures else {}, "detail_attribution": scenarios,
            "whole_run_resources": memory_reports}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", type=Path)
    parser.add_argument("--remove", action="append", default=[], metavar="REGION_ID:FRACTION",
                        help="Explicit hypothetical self-cost removal, e.g. 6:1 --remove 7:1 --remove 9:1")
    args = parser.parse_args()
    removals = [(int(i), float(f)) for i, f in (x.split(":") for x in args.remove)]
    if len({i for i, _ in removals}) != len(removals) or any(i < 1 or not 0 <= f <= 1 for i, f in removals):
        parser.error("Unique positive region IDs and fractions in [0,1] required")
    result = summarize(args.experiment.resolve(), removals)
    print(json.dumps(result, indent=2))
    if result["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
