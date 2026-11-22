"""Compare adaptive multilayer Drake checkpoints across a fresh MPI restart."""
import argparse
import filecmp
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def run(executable, grids, template, directory, resume):
    (directory / "grids").symlink_to(grids)
    settings = dict(run_id="drake", max_level=6, zlevels=6, dt_write=0.004,
                    time_end=0.010, resume=resume)
    for key, value in settings.items():
        template, count = re.subn(rf"(?m)^({key}\s+)\S+", rf"\g<1>{value}", template)
        if count != 1:
            raise ValueError(f"Expected one {key} entry in test_case.in")
    (directory / "test_case.in").write_text(template)
    result = subprocess.run(
        ["mpirun", "-n", "2", executable, "test_case.in"], cwd=directory,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300,
    )
    if result.returncode:
        raise RuntimeError(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path, help="J5 Drake executable")
    parser.add_argument("grids", type=Path, help="directory containing the optimized grids")
    args = parser.parse_args()
    executable, grids = args.executable.resolve(), args.grids.resolve()
    template = Path(__file__).with_name("test_case.in").read_text()
    with tempfile.TemporaryDirectory() as tmp:
        warm, fresh = Path(tmp) / "continuous", Path(tmp) / "fresh"
        warm.mkdir()
        fresh.mkdir()
        run(executable, grids, template, warm, -1)
        shutil.copy2(warm / "drake_checkpoint_0001.bin.zst", fresh)
        run(executable, grids, template, fresh, 1)
        checkpoint = "drake_checkpoint_0002.bin.zst"
        expected = subprocess.check_output(["zstd", "-dc", warm / checkpoint])
        actual = subprocess.check_output(["zstd", "-dc", fresh / checkpoint])
        if actual != expected:
            raise AssertionError("Fresh restart changed Drake checkpoint 2")
        if not filecmp.cmp(warm / checkpoint, fresh / checkpoint, shallow=False):
            raise AssertionError("Fresh restart changed the compressed Drake archive")
        print(f"Drake checkpoint 2 matched byte for byte ({len(actual)} bytes).")


if __name__ == "__main__":
    main()
