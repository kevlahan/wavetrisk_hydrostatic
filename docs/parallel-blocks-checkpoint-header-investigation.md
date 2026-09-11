# Cluster checkpoint header investigation

Cluster job 54449 passed Stage 187 checked production (60.75 s) and checked
oracle-on (183.45 s), including their exact checkpoint comparison. The optimized
pre-187 baseline then completed sufficiently to enter its checkpoint gate, but
the independent reader rejected the checkpoint before a candidate comparison.
The user's direct inspection of decompressed checkpoint 4 reports:

| Run | Compressed bytes | Magic | Version | Domains |
| --- | ---: | --- | ---: | ---: |
| checked-off | 17435925 | 5741564554524953 | 1 | 160 |
| checked-on | 17435925 | 5741564554524953 | 1 | 160 |
| baseline-warmup | 18041373 | 0000000000000000 | 0 | 0 |

The zeroed on-disk header is invalid; no repair, tolerance relaxation or timing
acceptance is justified. The remaining payload has not yet been independently
validated. Optimized Stage 187 has not yet run in this campaign.

The checkpoint writer is shared by the pre-187 baseline and candidate. It calls
collective MPI_File_set_size(0), then rank zero immediately writes the header
and directory with independent writes. A barrier precedes collective payload
writes, but no explicit barrier separates truncation from the header write.
A delayed truncation is a hypothesis to test, not an established root cause.
The [MPI file-size semantics](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report/node394.htm)
describe resizing as a write operation when applying consistency semantics.

An immediate solver restart is not independent proof that the saved header is
valid: dump_adapt_mpi broadcasts and retains the directory in memory, and
load_adapt_mpi can reuse it without rereading the disk header. Its extent check
compares physical file size with the retained directory. The separate Python
comparison therefore supplies a check that this path can miss.

## Standalone probe

`test/parallel_block_profile/checkpoint_io_probe.f90` exercises the same I/O
ordering: collective truncate, four independent metadata writes on rank zero,
barrier, collective per-rank payload writes, close and readback using Fortran
stream I/O. It checks MPI transfer counts, header, directory and every payload.
These are synthetic data files, not valid climate checkpoints. Mode 0 preserves
the current order; mode 1 adds a post-truncation barrier. Eight trials per mode
delay a nonzero rank before truncation to expose non-synchronizing behavior if
the MPI implementation permits it. No corruption is artificially injected.

A four-rank local run with Open MPI 5.0.8, O2 and bounds checking passes both
modes. The smaller preliminary header-only probe also passed. Thus the local
runs have not reproduced the cluster failure or established that the extra
barrier fixes it. Multi-node cluster/filesystem behavior must be measured.

`checkpoint_io_probe.sbatch` requests 83 tasks in partition bb, compiles only the
small diagnostic on an allocated node, and runs in a fresh home directory.
It records MPI/compiler versions and header/payload results. Submit from shared
home after transferring both probe files to the kit directory:

```sh
sbatch "$kit/checkpoint_io_probe.sbatch" "$kit"
```

Return the complete `checkpoint-io-probe-JOBID.log`. If mode 0 fails and mode 1
passes, that supports the truncation-order hypothesis on this platform. If both
pass, continue investigating the actual writer and filesystem; do not treat the
probe as proof of climate checkpoint correctness. No production Fortran source
has been changed in this investigation.

## Cluster probe 54450 and expanded diagnostic

The user supplied the complete first-probe log: GCC 13.2.0, MVAPICH 4.0 with
ch4:ucx and ROMIO, 83 ranks on bb01 and bb02. Current ordering fails two headers
(trials 5 and 8), zero directories and one rank payload block. Barrier-only
ordering fails zero headers/directories and two rank payload blocks. Thus the
header signature is reproduced, but barrier-only is not an adequate fix.
These counts are captured in `parallel-blocks-checkpoint-io-54450-results.json`.

The revised probe has six modes, eight trials each:

| Mode | Write path |
| --- | --- |
| 0 | Original ordering |
| 1 | Barrier after truncation |
| 2 | Mode 1 plus collective file synchronization after metadata and payload |
| 3 | Mode 2 plus requests to disable ROMIO collective buffering/data sieving |
| 4 | Mode 2 with independent disjoint payload writes |
| 5 | Gather tiny synthetic payloads and use one Fortran stream writer |

Mode 3 reports the actual returned hints; a requested hint is not assumed to
be honored. [ROMIO documents hint queries and ignored unsupported hints](https://ftp.mcs.anl.gov/pub/romio/users-guide/node6.html).
Mode 2 uses synchronization based on the [MPI file-consistency semantics](https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node296.htm).
Mode 5 is a diagnostic control, not a proposal to gather large climate
checkpoints into one rank's memory.

Rank zero and the last rank each verify the files after close, reporting the
reader host, bad block, first bad absolute offset, wrong-byte count and expected
and observed byte. The job records mount type/options on both allocated nodes.
An independent Python checker rereads all files after MPI exits, writes
`postrun-verification.json`, and reports persistent failures. Passing after a
previous immediate failure is evidence relevant to visibility, not grounds to
discard the earlier failure. An optional second batch argument rechecks the
original probe directory before starting the new test:

```sh
sbatch "$kit/checkpoint_io_probe.sbatch" "$kit" "$HOME/checkpoint-io-probe-54450"
```

All six modes pass eight trials with both readers locally (four ranks, Open MPI
5.0.8, O2 and bounds checks). The independent checker also passes all 48 files.
The focused checker test rejects zero headers, a wrong payload byte and missing
trailing bytes. ROMIO hints are not reported by the local Open MPI backend, so
local mode 3 does not validate ROMIO's hint behavior. Cluster execution is pending.
