# Stage 194: split inverse kernel timing and reuse inner-vector geometry

Stage 194 combines the measurement change with an opt-in computational change.
The previous-stage baseline is the exact Stage 193 source patch validated by
cluster job 54646. Every block case keeps Stage 193 scalar columns enabled and
Stage 192 sparse scatter disabled. Legacy remains a direct third comparison.

## Separate the measured work

The historical detail region 36, labelled "scalar/inner inverse kernels", also
encloses the outer-vector kernel. Consequently the Stage 193 reduction of
52.33% remains valid, but its residual cannot be attributed to inner-vector
work alone. Earlier review prose has been corrected; raw evidence is preserved.

Two appended detail regions retain all existing IDs:

- 66: scalar inverse kernels, covering both scalar phases.
- 67: inner-vector inverse kernel.
- Existing 21: outer-vector inverse kernel.
- Existing 36: inclusive combined parent, with its exclusive overhead reported
  separately. Do not add the inclusive parent to its children.

The analyzer checks that child times plus parent overhead conserve the parent,
that call counts cover all three children, and that complete ordinary-step
self times conserve the transaction. It reports whole-window rank means and
separate ordinary-step diagnostics. Timers surround producer calls rather than
individual cells or layers.

## Prepare horizontal geometry once per child

The inner-vector inverse formerly resolved the same 32 edge-length and
triangle-area values for each physical layer of each spatial cell. With
`WAVETRISK_STAGE194_INNER_GEOMETRY=1`, it resolves those values once per child
patch and reuses them over its physical column. The default remains off.

The cache is a small local `geometry(32,PATCH_SIZE/2,PATCH_SIZE/2)` array, 1 KiB
for the current 4-by-4 patch. It is rebuilt on every child and invocation, so it
does not persist through adaptation, remapping or restart. Explicit fixed slots
avoid dependence on expression evaluation order. Existing geometry accessors
still validate each prepared value; empty physical columns are skipped.

The field traversal, reconstructed samples, arithmetic expressions and division
order are preserved. This removes repeated geometry addressing and validation;
it does not remove physical reconstruction or alter inverse MPI transfers.
All comparisons use FP_CONTRACT=off. The selector is validated collectively,
including disagreement between ranks.

## Local validation

Fresh O2 and O0 checked binaries were built from the current numerical source.
The four-rank J4/J6 case has 30 atmospheric levels and 10 soil levels, resumes
checkpoint 2, writes and reloads checkpoint 3, and finishes 12 steps at t=0.3308
with two vertical remaps. Optimized cache-on, optimized cache-off and checked
cache-on runs all passed exact checkpoint comparisons against the validated
Stage 193 run and saved oracle, matched the printed timestep workload, and
passed the fixed legacy screen. Checked and optimized checkpoints also matched
exactly. Full cluster phase-oracle validation is still required.

Measured work in the optimized controls:

| Counter | Cache off | Cache on |
|---|---:|---:|
| Physical cell updates | 12,866,880 | 12,866,880 |
| Geometry reads | 411,740,160 | 13,724,672 |
| Prepared spatial cells | 0 | 428,896 |

Geometry reads fall **96.67%**, as expected for reuse across 30 layers. Actual
accessor counts in the disabled control match the reference count. Inverse
transfer calls, bytes and messages remain unchanged.

Single local instrumented runs gave these whole-window rank-mean times:

| Scope | Cache off | Cache on |
|---|---:|---:|
| Scalar | 0.593964 s | 0.628125 s |
| Inner-vector | 3.188420 s | 1.332318 s |
| Outer-vector | 1.603088 s | 1.683686 s |
| Combined parent, inclusive | 5.387000 s | 3.645750 s |

These are diagnostics on a workstation, not an accepted timestep speedup.
The targeted inner-vector time decreases, while unrelated scopes vary between
runs. Cluster acceptance depends on repeated uninstrumented timing.

The suite passed 128 tests with five optional MPI tests skipped. Separately,
both focused Stage 194 tests passed with MPI enabled. Focused Fortran coverage
executes the actual kernel at O0/O2, checks exact cached/uncached results,
physical/partial/empty columns, untouched soil, patch/boundary/ghost geometry,
cache invalidation between geometry epochs, invalid geometry, work counts,
and collective selector validation. Split-profile tests reject missing scopes,
unequal calls and failed timing conservation.

The machine-readable local record is `parallel-blocks-stage194-local-results.json`.
Exact checkpoints cover the save/restart boundary. Post-restart printed
workload and remapping are checked, but this short fixture has no final-state
checkpoint. External physics objects are shared; solver check flags do not
imply the external package was built with those checks.

## Cluster acceptance

The kit reconstructs the committed Stage 191 bootstrap, applies the exact frozen
Stage 193 baseline patch, builds that baseline, then applies only the Stage 194
MPI/profile delta. It verifies source hashes at each transition. Five fresh J5
builds include uninstrumented legacy and a separate legacy profile executable;
the latter retains its established 65-region schema.

There are 28 runs: two checked modes, three warmups, six rotating timing triples,
one disabled candidate control, and four profiles. The extra cache-disabled
Stage 194 profile provides identical split instrumentation for the cache-on
comparison. The original Stage 193 profile and uninstrumented baseline remain
unchanged. The campaign also checks identical cell work and inverse traffic.

Build and run on exclusive nk/bb13 with the established 83-rank J5/J7 fixture.
Follow `parallel-blocks-cluster-194.md`; collect one review archive before
transferring it to murier. Kit identities are in
`parallel-blocks-stage194-kit.json`. No Stage 194 cluster result is claimed yet.
