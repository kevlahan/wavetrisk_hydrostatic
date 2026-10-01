# Stage251: localize the early-feature J5 regression before Domain removal

Stage250 job 55240 showed that disabling the mirror-free candidate's
storage-diagnostic write gives no measured speedup. Its quiet candidate still
took 1.664006 s/ordinary step at 83 ranks on bb13, versus 0.541002 s for
read-only `main`: a 1.123003 s residual gap. The Stage248 feature-off build
was already about 1.629 s/step, so the large cost likely exists before the
later RK, physics, lifecycle and mirror-free features are enabled. That
cross-campaign comparison is a search lead, not a causal measurement.

The Stage251 kit uses the same sealed Stage248 candidate source for seven
serial J5 builds, all with `MIRROR_FREE=off`, `BLOCK_RK=off` and all later
native features off. The dependency ladder is:

| Build | Additional early feature enabled |
| --- | --- |
| `c0` | None; Stage202/203 block pathway only |
| `c1` | `BLOCK_COLUMNS` |
| `c2` | `BLOCK_PRIMITIVES` |
| `c3` | `BLOCK_DIFFUSION` |
| `c4` | `BLOCK_MASS_OWNER` |
| `c5` | `BLOCK_VELOCITY_OWNER` |
| `c6` | `BLOCK_TENDENCY_INITIAL` |

The `c6` selector set exactly matches the Stage248 feature-off baseline.
All builds share the same frozen source, compiler, FP contraction setting and
23 external Simple-physics objects; only the listed compile-time flags differ.
Read-only `main` is built on the same bb13 host as a numerical and timing
anchor. Each build must first pass a full 83-rank J5 restart: twelve steps,
two remaps, both VTK events, matching adaptive trajectory, physical-core
placement and the established main-relative checkpoint numerical gate.
The single-run medians screen adjacent stages. The kit then selects the
largest positive adjacent jump and confirms it with two A-B-B-A batches,
using full runs as independent samples. Steps 1, 8 and 9 are excluded from
ordinary-step medians. If a validation gate fails, the job stops with a
collectible diagnostic result instead of treating an invalid build as speed
evidence. No code is promoted automatically.

This is a localization experiment, not another Domain removal. Once the
costly transition is confirmed, inspect its loops and MPI waits and make a
narrow optimization while preserving the numerical gates. Keep `main`
unchanged. The [kit README](../test/parallel_block_profile/feature_ladder_251/README.md)
has no-edit transfer, launch and collection commands. Cluster timing is
pending.

The packaged kit was unpacked and verified locally. All eight production
executables (`main` and `c0`–`c6`) compiled in separate build directories
with `make -j1`, FP contraction disabled, and the same 23 external physics
objects. Local compilation validates the selectors and build workflow; only
the bb13 runs can establish numerical equivalence and performance.
