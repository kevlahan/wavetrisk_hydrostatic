# Checkpoint I/O workaround for bbserv

Probe job 54451 identifies a reproducible problem in the shared checkpoint I/O
path on bbserv's NFS4.2 home filesystem with MVAPICH 4.0/ROMIO. It does not
establish that one particular MPI or NFS component is defective. The original
probe failures persist when independently reread in a later job.

| Probe mode | Immediate verification | Independent post-run verification |
| --- | --- | --- |
| Current ordering | pass in this sample | pass in this sample; prior job had zero headers/payload holes |
| Post-truncation barrier | 3 bad rank-payload blocks on reader 0 | same 3 failures persist |
| Barrier plus file sync | 1 bad block on each reader | 1 failure persists |
| Sync plus ROMIO buffering/data-sieving disabled | both readers pass | all 8 files pass |
| Sync plus independent writes | both readers pass | all 8 files pass |
| Single stream writer | both readers pass | all 8 files pass |

The persistent holes start at offset 347288 and span 872 bytes at the beginning
of rank 42's payload. The second-reader-only failure starts at 344064 and spans
3224 bytes at the end of rank 41's payload. Those intervals partition the same
4096-byte region across the two-node rank boundary. This is evidence consistent
with interference at a shared page boundary; it is not a complete root-cause
proof. The two reader results also show why delayed checks alone are insufficient.

## Common writer patch

`test/parallel_block_profile/checkpoint_io_fix.patch` modifies only
`src/checkpoint.f90`. Apply this **same patch to both pre-187 and Stage 187**.
It does not change the integration scheme, field arithmetic, owner layout,
checkpoint format, compression settings or numerical comparison tolerance.

The writer now requests `romio_cb_write=disable` and `romio_ds_write=disable`
at file open, matching successful probe mode 3. It logs the hints actually
reported by MPI. If a reported value contradicts the request it stops. An MPI
implementation may not report unsupported ROMIO hints (the local Open MPI
backend does not), so a missing hint is logged and is not claimed to be honored.
[ROMIO documents these advisory settings and how to query them](https://ftp.mcs.anl.gov/pub/romio/users-guide/node6.html).

The patch also waits for every rank to finish truncation before writing the
header, synchronizes metadata and payload writes with barriers, checks all MPI
write transfer counts, and independently rereads the stored header/directory
after all handles close and before compression. That last check addresses the
immediate-restart path's reuse of directory data held in memory. It is a
metadata check, not a substitute for the runner's exact payload comparison.

These changes add work at checkpoint writes. Ordinary numerical timesteps do
not use the new functions. Both performance binaries must receive the patch;
old and new writer timings should not be mixed. No production-sized payload is
gathered onto one rank.

## Apply and rebuild on bbserv

Transfer `checkpoint_io_fix.patch` to the existing kit directory. Then:

```sh
kit="$HOME/wavetrisk-stage187-cluster-kit"
git -C "$HOME/wt-pre187" apply --check "$kit/checkpoint_io_fix.patch"
git -C "$HOME/wt-stage187" apply --check "$kit/checkpoint_io_fix.patch"
git -C "$HOME/wt-pre187" apply "$kit/checkpoint_io_fix.patch"
git -C "$HOME/wt-stage187" apply "$kit/checkpoint_io_fix.patch"
make -C "$HOME/wt-pre187" -j1 PARAM=param_J5 DEBUG=false BUILD_DIR=build-opt BIN_DIR=bin-opt
make -C "$HOME/wt-stage187" -j1 PARAM=param_J5 DEBUG=false BUILD_DIR=build-opt BIN_DIR=bin-opt
make -C "$HOME/wt-stage187" -j1 PARAM=param_J5 DEBUG=check BUILD_DIR=build-check BIN_DIR=bin-check
```

Use the same compute-compatible toolchain/flags as before. The corrected
physics module dependency includes should remain in place. The current runner
copies grids, expects checkpoint 4 and uses the fixture with time_end=0.315.
Resubmit the regular cluster campaign to a **new** output directory. The rebuilt
checked runs must pass before optimized warmups and timing begin. For each
checkpoint write, the bbserv log should report both hints as `disable`.
The exact checkpoint gate must still pass. Any new failure remains a blocker.

## Validation

The optimized and checked local solver builds pass. A four-rank MPI test extracts
and compiles the actual production output helpers. It passes metadata/payload
round-trip and fails as intended after deliberate header and directory damage.
`prepare_checkpoint_output_probe.py` preserves that test harness.

The two-step J4/J6 short RK4 restart fixture passes exact optimized comparison
against the unmodified writer (both about 25.74 seconds), and exact checked
comparison against the previous checked reference (77.37 seconds). Each writes
and reloads checkpoint 3; this short fixture has no post-reload timestep/remap
coverage. These are numerical validation runs, not performance measurements.
The patch applies cleanly to the identical original checkpoint source in both
comparison trees and reproduces its recorded SHA256. Machine-readable results
are in `parallel-blocks-checkpoint-io-workaround-results.json`. Cluster solver validation subsequently passed in job 54452: all checked and
optimized campaign gates completed, including exact checkpoint comparisons.
See [the cluster result](parallel-blocks-cluster-187-54452.md). This validates
the workaround for that campaign, without claiming universal I/O correctness. The known
legacy and checked-versus-optimized numerical discrepancies remain separate.
