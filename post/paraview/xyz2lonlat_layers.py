# xyz2lonlat_layers.py
#
# Converts all vertical layers of a *.vtk.tgz file to *.vtp longitude-latitude files
# Can be run in parallel using mpi (where ntasks <= number of layers) or in serial.
#
# To run in serial:
#   python3 xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 365 1 10
#
# To run using mpi (number of tasks <= nz)
#   srun -n 30 python -m mpi4py xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 1 10
#
# To use as a module:
#   from xyz2lonlat_layers import run_xyz_layers
#
# Calling:
#   run_xyz_layers(run, Jmin, Jmax, nz, t1, t2)
#
# Input parameters:
#   run = prefix for run without _tri (e.g. SimpleJ8Z60)
#   Jmin = minimum level
#   Jmax = maximum level
#   nz   = number of vertical layers
#   t1   = first time (save index)
#   t2   = last time  (save index)
#
# module load gcc mvapich python py-numpy py-scipy py-mpi4py vtk

from pathlib import Path
import subprocess
from mpi4py import MPI
import  argparse, os, tarfile, time, sys, textwrap

class _Fmt(argparse.ArgumentDefaultsHelpFormatter,
           argparse.RawDescriptionHelpFormatter):
    pass

def parse_args():
    epilog = textwrap.dedent("""\
      Examples
        mpirun -n 4 python -m mpi4py xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 1 10
        python xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 1 10
    """)
    p = argparse.ArgumentParser(
        prog="xyz2lonlat_layers.py",
        description="Convert all vertical layers of spherical <run>_tri_<t>.vtk.tgz to longitude-latitude data <run>_lonlat_<z>_<t>.vtp.",
        formatter_class=_Fmt,
        epilog=epilog,
    )
    p.add_argument("run")
    p.add_argument("Jmin", type=int)
    p.add_argument("Jmax", type=int)
    p.add_argument("nz",   type=int)
    p.add_argument("t1",   type=int)
    p.add_argument("t2",   type=int)
    p.add_argument("--wdir", type=Path, default=None, help="Working/data directory")
    p.add_argument(
        "--prep-quiet",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Suppress prep/untar messages",
    )

    if len(sys.argv) == 1:
        print()
        p.print_help()
        sys.exit(0)

    return p.parse_args()

def prep_once_rank0_then_barrier(comm, tgz_path: str, dest_dir: str,
                                 timeout_s: int = 300) -> str:
    """
    Rank 0 extracts tgz_path into dest_dir and writes dest/.ready when finished.
    All ranks wait (Barrier). Non-root ranks poll for .ready to cover NFS lag.
    Returns absolute dest_dir.
    """
    tgz  = Path(tgz_path).expanduser().resolve()
    dest = Path(dest_dir).expanduser().resolve()
    ready = dest / ".ready"

    if comm is None:
        # Serial: just extract (no barrier)
        dest.mkdir(parents=True, exist_ok=True)
        if not tgz.exists():
            raise FileNotFoundError(f"tar not found: {tgz}")
        with tarfile.open(str(tgz), "r:gz") as tar:
            tar.extractall(dest, filter="data")
        return str(dest)

    rank = comm.Get_rank()

    if rank == 0:
        dest.mkdir(parents=True, exist_ok=True)
        if not tgz.exists():
            print(f"[prep] ERROR: tar not found: {tgz}", file=sys.stderr, flush=True)
        else:
            with tarfile.open(str(tgz), "r:gz") as tar:
                tar.extractall(dest, filter="data")
            tmp = dest / ".ready.tmp"
            tmp.write_text("ok")
            os.replace(tmp, ready)

    comm.Barrier()

    if rank != 0:
        t0 = time.time()
        while not ready.exists():
            if time.time() - t0 > timeout_s:
                raise TimeoutError(f"Timeout waiting for {ready}")
            time.sleep(0.1)

    return str(dest)


for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")

def run_one(z, run, Jmin, Jmax, t, use_flag="y"):
    cmd = [sys.executable, "xyz2lonlat.py", run, str(Jmin), str(Jmax),
           str(z), str(z), str(t), str(t), use_flag]
    return subprocess.run(cmd).returncode

def try_mpi():
    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD
    except Exception:
        return None

def run_xyz_layers_serial(run, Jmin, Jmax, nz, t1, t2, prep_quiet):
    failures = []
    for t in range(t1, t2 + 1):
        if not prep_quiet:
            print(f"Processing file: {run}_tri_{t:04d}.vtk.tgz", flush=True)
        for z in range(1, nz + 1):
            rc = run_one(z, run, Jmin, Jmax, t)
            if rc != 0:
                failures.append((t, z, rc))
    return failures


def run_xyz_layers(run, Jmin, Jmax, nz, t1, t2, prep_quiet):
    comm = try_mpi()
    rank = comm.Get_rank() if comm else 0
    size = comm.Get_size() if comm else 1

    wdir = str(Path.cwd().resolve()) if rank == 0 else None
    if comm:
        wdir = comm.bcast(wdir, root=0)
    os.chdir(wdir)

    failures = []

    for t in range(t1, t2 + 1):
        tar_file = f"{run}_tri_{t:04d}.vtk.tgz"
        if rank == 0 and not prep_quiet:
            print(f"Processing file: {tar_file}", flush=True)

        tgz_path = str(Path(wdir) / tar_file)
        dest_dir = wdir

        # Untar once on rank 0, barrier for others
        try:
            _ = prep_once_rank0_then_barrier(comm, tgz_path, dest_dir)
        except Exception as e:
            if rank == 0:
                print(f"[prep] ERROR for t={t}: {e}", file=sys.stderr, flush=True)
            failures.append((t, 0, 1))
            continue

        # Remove surface data *.vtk file on rank 0
        if rank == 0:
            surface_vtk = f"{run}_tri_000_{t:04d}.vtk"
            p = Path(surface_vtk)
            if p.exists():
                p.unlink()

        # Work-share over z
        for z in range(1 + rank, nz + 1, size):
            rc = run_one(z, run, Jmin, Jmax, t)
            if rc != 0:
                print(f"[rank {rank}] run_one failed: t={t}, z={z}, rc={rc}", flush=True)
                failures.append((t, z, rc))

    if comm:
        all_fail = comm.gather(failures, root=0)
        if rank == 0:
            failures = [x for sub in all_fail for x in sub]
        else:
            failures = []

    return failures
if __name__ == "__main__":
    a = parse_args()

    if a.wdir is not None:
        os.chdir(a.wdir.resolve())

    comm = try_mpi()
    rank = comm.Get_rank() if comm else 0

    if rank == 0 and not a.prep_quiet:
        print(f"Converting layers for run={a.run}, J=[{a.Jmin},{a.Jmax}], nz={a.nz}, t1={a.t1}, t2={a.t2}")

    failures = run_xyz_layers(a.run, a.Jmin, a.Jmax, a.nz, a.t1, a.t2, a.prep_quiet)

    if (comm is None) or (rank == 0):
        if failures:
            print("Failures:", failures)
