# xyz2lonlat_layers.py
#
# Converts all vertical layers of a *.vtk.tgz file to *.vtp longitude-latitude files
# Can be run in parallel using mpi (where ntasks <= number of layers) or in serial.
#
# To run in serial:
#   python3 xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 365
#
# To run using mpi:
#   mpirun -n 4 python xyz2lonlat_layers.py SimpleJ5Z30 5 5 30 365
#
# To use as a module:
#   from xyz2lonlat_layers import run_xyz_layers
#
# Calling:
#   run_xyz_layers(run, Jmin, Jmax, nz, t)
#
# Input parameters:
#   run = prefix for run without _tri (e.g. SimpleJ8Z60)
#   Jmin = minimum level
#   Jmax = maximum level
#   nz   = number of vertical layers
#   t    = time (save index)

from pathlib import Path
import subprocess
from mpi4py import MPI
import os, tarfile, time, sys

def prep_once_rank0_then_barrier(comm: MPI.Comm, tgz_path: str, dest_dir: str,
                                 timeout_s: int = 300, verbose: bool = True) -> str:
    """
    Rank 0 extracts tgz_path into dest_dir (same dir as tar is fine).
    Writes dest/.ready when finished. All ranks then pass a Barrier.
    Non-root ranks poll for .ready after the Barrier to cover NFS lag.
    Returns absolute dest_dir.
    """
    rank = comm.Get_rank()
    tgz  = Path(tgz_path).expanduser().resolve()
    dest = Path(dest_dir).expanduser().resolve()
    ready = dest / ".ready"

    if rank == 0:
        dest.mkdir(parents=True, exist_ok=True)
        if not tgz.exists():
            print(f"[prep] ERROR: tar not found: {tgz}", file=sys.stderr, flush=True)
            # still enter barrier so others can receive the error broadcast below
        else:
            # Extract once into the destination directory
            with tarfile.open(tgz, "r:gz") as tar:
                tar.extractall(dest)
            # Publish a .ready marker (atomic rename helps NFS)
            tmp = dest / ".ready.tmp"
            tmp.write_text("ok")
            os.replace(tmp, ready)

    # Everyone waits until rank 0 finishes the extraction
    comm.Barrier()

    # Handle NFS attribute caching on non-root ranks
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

def run_xyz_layers(run, Jmin, Jmax, nz, t):
    """
    Runs z=1..nz. If launched under MPI (mpirun/srun with mpi4py), distributes
    work across ranks; otherwise runs serial. Returns a list of (z, rc) failures.
    """
    comm = try_mpi()
    rank = comm.Get_rank() if comm else 0
    size = comm.Get_size() if comm else 1

    # Ensure all ranks agree on the working directory
    # (either run mpirun with --wdir /abs/path OR broadcast Path.cwd())
    wdir = str(Path.cwd().resolve()) if rank == 0 else None
    wdir = comm.bcast(wdir, root=0)
    os.chdir(wdir)

    # Build tar name and destination = SAME DIRECTORY AS THE TAR
    tar_file = f"{run}_tri_{t:04d}.vtk.tgz"
    tgz_path = str(Path(wdir) / tar_file)
    dest_dir = wdir  # same directory

    # Untar once on rank 0, wait for all
    data_root = prep_once_rank0_then_barrier(comm, tgz_path, dest_dir)

    # Remove surface data *.vtk file
    surface_vtk = f"{run}_tri_000_{t:04d}.vtk"
    p = Path(surface_vtk).expanduser()
    if p.exists():
        os.remove(surface_vtk)

    failures = []
    for z in range(1 + rank, nz + 1, size):
        rc = run_one(z, run, Jmin, Jmax, t)
        if rc != 0:
            failures.append((z, rc))

    if comm:
        all_fail = comm.gather(failures, root=0)
        if rank == 0:
            failures = [x for sub in all_fail for x in sub]

    # return failures on rank 0 (or always, if serial)
    if (comm and rank == 0) or (comm is None) or size == 1:
        return failures
    else:
        return []  # non-root ranks return an empty list

# CLI wrapper (still works from terminal)
if __name__ == "__main__":
    run, Jmin, Jmax, nz, t = sys.argv[1:6]
    fails = run_xyz_layers(run, int(Jmin), int(Jmax), int(nz), int(t))
    if fails:
        print("Failures:", fails, file=sys.stderr)
        sys.exit(1)
