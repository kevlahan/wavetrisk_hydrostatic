#!/usr/bin/env python3
"""Compile an archived main revision, never editing tracked legacy source.

Only generates the physics package's untracked build-dependency file. Uses
the archived Makefile as-is, and verifies tracked contents after compilation.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess

from experiment import LEGACY_REF, digest


def tracked_entries(repo, ref):
    data = subprocess.check_output(["git", "-C", str(repo), "ls-tree", "-rz", ref])
    result = []
    for item in data.split(b"\0"):
        if item:
            header, name = item.split(b"\t", 1)
            mode, kind, sha = header.decode().split()
            if kind != "blob":
                raise ValueError("Archive contains a submodule; resolve build inputs explicitly")
            result.append((mode, sha, name.decode()))
    return result


def verify(root, entries):
    errors = []
    for mode, expected, name in entries:
        path = root / name
        data = str(path.readlink()).encode() if mode == "120000" else path.read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if actual != expected:
            errors.append(name)
    if errors:
        raise ValueError(f"Tracked legacy contents changed: {errors}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New directory, outside the repository")
    parser.add_argument("--ref", default=LEGACY_REF)
    parser.add_argument("--mpif90", default="mpif90")
    args = parser.parse_args()
    repo, root = args.repo.resolve(strict=True), args.out.resolve()
    if root.is_relative_to(repo) or repo.is_relative_to(root):
        parser.error("Build snapshot must be outside the repository tree")
    ref = subprocess.check_output(["git", "-C", str(repo), "rev-parse", args.ref + "^{commit}"], text=True).strip()
    entries = tracked_entries(repo, ref)
    root.mkdir(parents=True, exist_ok=False)
    archive = subprocess.Popen(["git", "-C", str(repo), "archive", ref], stdout=subprocess.PIPE)
    try:
        subprocess.run(["tar", "-x", "-C", str(root)], stdin=archive.stdout, check=True)
    finally:
        archive.stdout.close()
    if archive.wait():
        raise SystemExit("git archive failed")
    verify(root, entries)
    physics = root / "src/physics/simple_physics/phyparam"
    dep = physics / "Makefile.inc"
    if not dep.exists():
        names = sorted(p.stem for p in (physics / "physics").glob("*.F90"))
        combined = subprocess.check_output(["bash", "bash/concatenate_all_code.sh", *names], cwd=physics)
        dependencies = subprocess.check_output(["bash", "bash/makedeps.sh", "/dev/stdin"], cwd=physics, input=combined)
        with dep.open("xb") as output:
            output.write(dependencies)
    command = ["make", "-j1", "DEBUG=false", "PARAM=param_J5", "TEST_CASE=climate", "MPIF90=" + args.mpif90]
    with (root / "build.log").open("x") as output:
        output.write(f"Unmodified legacy revision: {ref}\nCommand: {command}\n")
        output.flush()
        result = subprocess.run(command, cwd=root, stdout=output, stderr=subprocess.STDOUT)
    verify(root, entries)
    if result.returncode:
        raise SystemExit(f"Build failed; tracked sources verified unchanged. See {root / 'build.log'}")
    identity = {"revision": ref, "tracked_contents_unchanged": True,
                "binary_sha256": digest(root / "bin/climate"), "command": command,
                "build_host": platform.uname()._asdict()}
    (root / "legacy-build.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
