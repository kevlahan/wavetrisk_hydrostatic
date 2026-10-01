# Stage248 bb13 diagnostic job 55185

The diagnostic archive `/Users/kevlahan/Desktop/rk248-diagnostic-55184-55185.tar.gz` has SHA-256 `6d2e4468c1ec820050212dbddaf5c00f0737029838739076e4ee9dc02a3e5853` and 43 members. It reuses the exact failed job 55184 J4 feature-off input, 16-rank placement, and optimized baseline binary. The stack limit was raised from 8192 KiB to unlimited, but the optimized run again segfaulted at startup: ranks 13–14 exited with signal 11, the remaining ranks were terminated, and `srun` returned 143.

A new `DEBUG=check` build of the same baseline compile-time selectors succeeded without build warnings. Its 16-rank run also segfaulted at the same startup point: ranks 3–4 exited with signal 11, the remaining ranks were terminated, and `srun` returned 143. Neither run recorded an ordinary step, and neither Fortran signal backtrace supplied a useful source frame. The simple stack-limit explanation is therefore ruled out. The exact invalid access remains unlocalized.

`test/parallel_block_profile/performance_248/diagnose_asan.sbatch` is a further bounded diagnostic: a fresh ASan/UBSan build and one 16-rank rerun of the same baseline. It packages logs and rank records without using timings for performance acceptance. Keep the Stage248 performance campaign stopped until the crash is located and corrected.
