"""Exercise the J5 upwelling checkpoint writer with two MPI ranks and 160 domains."""
import argparse
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile


def run(executable, cwd, env, expected_error=None):
    result = subprocess.run(
        ["mpirun", "-n", "2", executable, "test_case.in"], cwd=cwd, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300,
    )
    if expected_error is None:
        passed = result.returncode == 0
    else:
        passed = result.returncode != 0 and expected_error in result.stdout
    if not passed:
        raise RuntimeError(result.stdout)


def configure(directory, grids, template, resume, **overrides):
    directory.mkdir()
    (directory / "grids").symlink_to(grids)
    values = dict(max_level=5, zlevels=4, dt_write=0.01, time_end=0.03, resume=resume)
    values.update(overrides)
    for key, value in values.items():
        template, count = re.subn(rf"(?m)^({key}\s+)\S+", rf"\g<1>{value}", template)
        assert count == 1, key
    (directory / "test_case.in").write_text(template)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path, help="J5 upwelling executable")
    parser.add_argument("grids", type=Path, help="directory containing the optimized grids")
    args = parser.parse_args()
    executable, grids = args.executable.resolve(), args.grids.resolve()
    env = os.environ.copy()
    zstd = shutil.which("zstd")
    if not zstd:
        parser.error("zstd is required")
    template = Path(__file__).with_name("test_case.in").read_text()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        writer, reader = root / "writer", root / "reader"
        configure(writer, grids, template, -1)
        configure(reader, grids, template, 0)
        run(executable, writer, env)
        name = "upwelling_checkpoint_0000.bin.zst"
        archive = (writer / name).read_bytes()
        (reader / name).write_bytes(archive)
        run(executable, reader, env)

        # Flip one domain's phase while preserving a valid compression frame.
        raw = bytearray(subprocess.check_output([zstd, "-dc", writer / name]))
        magic, version, domains = struct.unpack_from("=3q", raw)
        assert (magic, version, domains) == (0x5741564554524953, 2, 160)
        offset = struct.unpack_from("=q", raw, 24 + 4 * domains + 8 * 54)[0]
        phase = struct.unpack_from("=i", raw, offset + 24)[0]
        struct.pack_into("=i", raw, offset + 24, phase + 1)
        subprocess.run([zstd, "-qf", "-o", reader / name], input=raw, check=True)
        run(executable, reader, env, "checkpoint domain clocks or remap phases disagree")

        # Damage the verification read-back before the rename. The previous
        # published archive must remain byte-identical and usable.
        wrapper = root / "zstd"
        wrapper.write_text(
            '#!/usr/bin/env python3\n'
            'import os, subprocess, sys\n'
            'args = sys.argv[1:]\n'
            'subprocess.run([os.environ["CHECKPOINT_TEST_ZSTD"], *args], check=True)\n'
            'if "-df" in args:\n'
            '    with open(args[args.index("-o") + 1], "r+b") as f:\n'
            '        f.seek(4000)\n'
            '        value = f.read(1)[0]\n'
            '        f.seek(4000)\n'
            '        f.write(bytes([value ^ 1]))\n'
        )
        wrapper.chmod(0o755)
        injected = dict(env, PATH=f"{root}{os.pathsep}{env['PATH']}", CHECKPOINT_TEST_ZSTD=zstd)
        run(executable, writer, injected, "checkpoint archive failed byte verification")
        assert (writer / name).read_bytes() == archive
        # Evolve a J5-J6 grid, then compare automatic and fresh continuation
        # from the same nonzero checkpoint, including the adaptive tree.
        warm, cold = root / "adaptive", root / "adaptive_restart"
        settings = dict(max_level=6, dt_write=0.04, time_end=0.09)
        configure(warm, grids, template, -1, **settings)
        configure(cold, grids, template, 1, **settings)
        run(executable, warm, env)
        shutil.copy2(warm / "upwelling_checkpoint_0001.bin.zst", cold)
        run(executable, cold, env)
        evolved = "upwelling_checkpoint_0002.bin.zst"
        expected = subprocess.check_output([zstd, "-dc", warm / evolved])
        actual = subprocess.check_output([zstd, "-dc", cold / evolved])
        assert actual == expected, "adaptive restart changed the next checkpoint"
    print("Passed checkpoint save, fresh restart, domain phase rejection, failed publication, "
          "and adaptive continuation tests.")


if __name__ == "__main__":
    main()
