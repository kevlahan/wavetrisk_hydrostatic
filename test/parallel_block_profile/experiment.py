#!/usr/bin/env python3
"""Matched Slurm runs in fresh fixtures. No source changes, commits or overwrites.

Default: six alternating-order legacy/block pairs, then two off/detail pairs.
Use --sample for a separate paired external perf run, not a timing observation.
Use --independent-srun from a normal shell to allocate each run separately.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze import analyze

LEGACY_REF = "b82647b36eb825aed0c0aae4e92adcf4ab46d608"
SWITCHES = ("WAVETRISK_VALIDATE_BLOCK_DYNAMICS", "WAVETRISK_VALIDATE_BLOCK_ADAPTATION",
            "WAVETRISK_VALIDATE_BLOCK_REMAP",
            "WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS", "WAVETRISK_PROFILE_PARALLEL_BLOCKS",
            "WAVETRISK_PROFILE_BLOCK_DETAIL")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def manifest(root):
    """Hash regular input files, including symlinked files; identify directory links."""
    result = {}
    for path in sorted(root.rglob("*")):
        key = str(path.relative_to(root))
        if path.is_file():
            result[key] = {"sha256": digest(path), "bytes": path.stat().st_size}
        elif path.is_symlink():
            if not path.exists():
                raise ValueError(f"Broken fixture link: {path}")
            result[key] = {"directory_link": str(path.resolve()), "contents_hashed": False}
    return result


def schedule(pairs, detail_pairs, sample):
    if sample:
        return [("sample", 0, "legacy"), ("sample", 0, "block")]
    result = []
    for group, count, labels in (("timing", pairs, ("legacy", "block")),
                                 ("overhead", detail_pairs, ("block", "detail"))):
        for pair in range(count):
            order = labels if pair % 2 == 0 else labels[::-1]
            result.extend((group, pair, label) for label in order)
    return result


def validate_launch_mode(independent, env):
    if independent and env.get("SLURM_JOB_ID"):
        raise ValueError("--independent-srun must run from a normal shell outside an allocation; "
                         "omit this flag inside an existing job")
    if not independent and not env.get("SLURM_JOB_ID"):
        raise ValueError("Use --independent-srun for separate srun jobs from a normal shell, "
                         "or run inside an existing allocation")


def launch_command(ranks, python, helper, placement, application, partition=None, nodelist=None, direct=False):
    command = (["srun", "-n", str(ranks)] if direct else
               ["srun", "--ntasks", str(ranks), "--cpus-per-task=1", "--cpu-bind=cores",
                "--export=ALL", "--kill-on-bad-exit=1"])
    if partition:
        command.append("--partition=" + partition)
    if nodelist:
        command.append("--nodelist=" + nodelist)
    if direct:
        return [*command, *application]
    return [*command, python, str(helper / "rank_launch.py"), str(placement), *application]


def execute_with_tee(command, directory, env, output):
    """Equivalent to 2>&1 | tee, preserving srun's exit status."""
    with subprocess.Popen(command, cwd=directory, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, errors="replace") as process:
        for line in process.stdout:
            print(line, end="", flush=True)
            output.write(line)
            output.flush()
        return process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--block", type=Path, required=True)
    parser.add_argument("--legacy-build-log", type=Path, required=True)
    parser.add_argument("--legacy-identity", type=Path, required=True, help="legacy-build.json from build_legacy.py")
    parser.add_argument("--block-build-log", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True, help="Small frozen restart input directory")
    parser.add_argument("--out", type=Path, required=True, help="Must not exist")
    parser.add_argument("--ranks", type=int, default=83)
    parser.add_argument("--pairs", type=int, default=6)
    parser.add_argument("--detail-pairs", type=int, default=2)
    parser.add_argument("--sample", choices=("record", "stat"))
    parser.add_argument("--sample-ranks", default="0,30,52,82")
    parser.add_argument("--independent-srun", action="store_true",
                        help="Copy climateJ5 and launch each standalone srun directly, without rank wrappers")
    parser.add_argument("--partition", help="Optional partition for each srun; default uses normal Slurm selection")
    parser.add_argument("--nodelist", help="Optional node list for each srun; no nodes are guessed")
    args = parser.parse_args()
    try:
        validate_launch_mode(args.independent_srun, os.environ)
    except ValueError as error:
        parser.error(str(error))
    if args.ranks < 1 or args.pairs < 1 or args.detail_pairs < 0:
        parser.error("Invalid ranks/pair counts")
    if args.sample and any(int(r) not in range(args.ranks) for r in args.sample_ranks.split(",")):
        parser.error("Sample rank outside allocation")
    args.fixture = args.fixture.resolve(strict=True)
    if not (args.fixture / "simple.in").is_file():
        parser.error("Fixture must contain simple.in and its restart inputs")
    reserved = ("run.log", "RK4_profile.log", "RK4_production.log", "rank-placement", "samples")
    if args.independent_srun:
        reserved += ("climateJ5",)
    if any((args.fixture / name).exists() or (args.fixture / name).is_symlink() for name in reserved):
        parser.error("Use an input-only seed: remove executable/logs/rank-placement/samples from the fixture")
    args.out = args.out.resolve()
    if args.out.is_relative_to(args.fixture) or args.fixture.is_relative_to(args.out):
        parser.error("Output and fixture trees must be disjoint")
    binaries = {label: path.resolve(strict=True) for label, path in
                (("legacy", args.legacy), ("block", args.block))}
    identities = {label: digest(path) for label, path in binaries.items()}
    legacy_identity = json.loads(args.legacy_identity.read_text())
    if (legacy_identity.get("revision") != LEGACY_REF or not legacy_identity.get("tracked_contents_unchanged")
            or legacy_identity.get("binary_sha256") != identities["legacy"]):
        parser.error("Legacy identity does not match the requested unchanged main revision and executable")
    inputs = manifest(args.fixture)
    args.out.mkdir(parents=True, exist_ok=False)
    for label, path in (("legacy", args.legacy_build_log), ("block", args.block_build_log)):
        shutil.copy2(path, args.out / f"{label}-build.log")
    metadata = {"schema": 1, "legacy_reference": LEGACY_REF, "legacy_build_identity": legacy_identity,
                "launch_mode": "independent-srun" if args.independent_srun else "shared-allocation",
                "binary_paths": {k: str(v) for k, v in binaries.items()}, "binary_sha256": identities,
                "fixture": inputs, "ranks": args.ranks,
                "allocation": {k: v for k, v in os.environ.items() if k in
                    ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK")},
                "runs": []}
    helper = Path(__file__).resolve().parent
    # The identity wrapper needs a shared Python path on all allocated nodes.
    python = shutil.which("python3")
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env.update({name: "0" for name in SWITCHES})
    thread_settings = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                             "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    env.update(thread_settings)
    metadata["thread_settings"] = thread_settings
    metadata["mpi_settings"] = {k: env[k] for k in ("MPICH_ASYNC_PROGRESS", "I_MPI_ASYNC_PROGRESS",
        "OMPI_MCA_mpi_yield_when_idle", "UCX_TLS", "FI_PROVIDER") if k in env}
    # Refuse unknown inherited WAVETRISK switches; do not silently run an oracle.
    unknown = [k for k in env if k.startswith("WAVETRISK_") and k not in SWITCHES
               and not k.startswith("WAVETRISK_SAMPLE_")]
    if unknown:
        raise SystemExit(f"Unset unclassified WAVETRISK switches before comparison: {unknown}")
    reference = None
    for index, (group, pair, label) in enumerate(schedule(args.pairs, args.detail_pairs, args.sample)):
        binary_label = "legacy" if label == "legacy" else "block"
        binary = binaries[binary_label]
        if digest(binary) != identities[binary_label] or manifest(args.fixture) != inputs:
            raise SystemExit("Binary or input fixture changed during experiment")
        directory = args.out / f"{index:02d}-{group}-{pair}-{label}"
        shutil.copytree(args.fixture, directory, symlinks=True)
        # Do not permit checkpoint output to overwrite a symlinked seed file.
        for path in directory.rglob("*"):
            if path.is_symlink() and path.is_file():
                content = path.resolve()
                path.unlink()
                shutil.copy2(content, path)
        placement = directory / "rank-placement"
        placement.mkdir()
        run_env = env.copy()
        if label == "detail":
            run_env["WAVETRISK_PROFILE_PARALLEL_BLOCKS"] = "1"
            run_env["WAVETRISK_PROFILE_BLOCK_DETAIL"] = "1"
        application = [str(binary), "simple.in"]
        if args.independent_srun:
            local_binary = directory / "climateJ5"
            if local_binary.exists() or local_binary.is_symlink():
                raise SystemExit("Remove climateJ5 from the frozen seed; the driver copies the correct binary for each run")
            shutil.copy2(binary, local_binary)
            if digest(local_binary) != identities[binary_label]:
                raise SystemExit("Copied executable hash differs")
            application = ["./climateJ5", "simple.in"]
        if args.sample:
            run_env.update(WAVETRISK_SAMPLE_RANKS=args.sample_ranks, WAVETRISK_SAMPLE_MODE=args.sample,
                           WAVETRISK_SAMPLE_DIR=str(directory / "samples"))
            application = ["/bin/sh", str(helper / "sample_rank.sh"), *application]
        command = launch_command(args.ranks, python, helper, placement, application,
                                 args.partition, args.nodelist, direct=args.independent_srun)
        record = {"group": group, "pair": pair, "label": label, "directory": directory.name,
                  "command": command, "switches": {k: run_env[k] for k in SWITCHES}, "start_unix": time.time()}
        print(f"Starting {directory.name}", flush=True)
        print(" ".join(command), flush=True)
        with (directory / "run.log").open("x") as output:
            returncode = execute_with_tee(command, directory, run_env, output)
        record.update(returncode=returncode, end_unix=time.time())
        log_name = "RK4_profile.log" if label == "detail" else "RK4_production.log"
        (directory / log_name).symlink_to("run.log")
        record["allocation_job_ids"] = sorted({item["slurm"]["SLURM_JOB_ID"]
            for path in placement.glob("rank-*.json") for item in [json.loads(path.read_text())]
            if item.get("slurm", {}).get("SLURM_JOB_ID")})
        report = analyze(directory / "run.log")
        state = [" ".join(s.split()) for s in report["printed_states"]]
        if reference is None:
            reference = state
        record["printed_states_match"] = bool(state) and state == reference
        record["completed"] = report["completed"]
        metadata["runs"].append(record)
        # Owned output only. Preserve a usable manifest even when a run fails.
        (args.out / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if returncode or not report["completed"] or not record["printed_states_match"]:
            raise SystemExit(f"Stopped after failed/inconsistent run: {directory}; inspect before comparing timings")
        if label == "detail" and not any(w.get("region_count") == 65 for w in report["detail_windows"]):
            raise SystemExit("Missing Stage 179b detail schema; check that the revised executable was copied")
    print(f"Completed. Analyze with: python3 {helper / 'compare.py'} {args.out}")


if __name__ == "__main__":
    main()
