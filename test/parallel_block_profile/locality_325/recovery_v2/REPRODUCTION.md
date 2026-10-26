# Reproduce the recovery-v2 package

Commit this directory as `test/parallel_block_profile/locality_325/recovery_v2`.
Its recipe references 13 shared Stage325 helpers from the parent directory,
so those helpers do not need duplicate copies. The 16 unique payload files
contain the recovery runner/controller, shell wrappers, site allocation,
historical input bindings, recovery tests, and retained validation evidence.

The sealed package used a derived-map boundary test fixture; that test is
retained here as an explicit recovery override. `verify322.py` differs from
the shared helper only by one final newline, which the recipe restores.
Every recovered payload byte must match the deployed sealed manifest.

From a checkout, with a nonexistent output directory:

```bash
python3 test/parallel_block_profile/locality_325/recovery_v2/build_recovery325.py \
  --output /tmp/rk325-recovery-reproduction \
  --run-tests
```

The builder verifies all hashes, writes a fresh package, runs its sealed
checker and optionally the 91 offline tests, then creates a tar.gz with
normalized ownership, mode, timestamps, and gzip header. Repeated builds
in the same Python/zlib environment produce identical archives. The
normalized archive checksum differs from the original deployment archive;
the packaged files and sealed manifest are identical. Both checksums are
recorded in the reproduction report. There is no scheduler submission,
compilation, solver execution, or remote access in this builder.

The v2 verifier is frozen here at its deployed hash. Later checker repairs in
the parent directory do not change reproduction of the failed v2 replay.
