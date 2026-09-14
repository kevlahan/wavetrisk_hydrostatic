# Stage 193: scalar inverse computation shared across fields

The Stage 192 experiment reduced inverse scatter bytes by 7.03% on bb13 but
changed the ordinary-step median by only 0.37%. Its complete wavelet transaction
remained 4.37 times legacy. This stage targets repeated computation, using the
accepted Stage 191 configuration as its production baseline. Stage 192 sparse
scatter remains available explicitly, but its default is now off.

## Implemented change

`reconstruct_block_scalar_columns` implements both scalar inverse phases:
undoing scalar lifting and reconstructing the three detail nodes per coarse
cell. The reference routine performs stencil location, storage classification,
area/overlap access and frozen-mask lookup inside the scalar/vertical-layer
loops. Those quantities do not vary across physical fields.

The new routine prepares each spatial stencil once, then reads and updates
physical columns. Patch, boundary and ghost field addresses are resolved once
per source, and immutable geometry is read once per stencil. The per-layer
weighted sum still visits its terms in the same order; no weights are combined,
no division or multiplication is reassociated, and no arithmetic tolerance is
relaxed. With two scalars and 30 physical layers, this removes 59 of every 60
repeated address and geometry preparations in this kernel. It does not remove
59 of 60 arithmetic updates or claim a 60-fold solver speedup.

The routine retains finite-input/output checks, storage bounds, frozen-node
handling, atmospheric-only updates and the original coverage counters. It uses
small column scratch arrays allocated once per block callback and fixed stencil
arrays, with no persistent geometry/address cache to
invalidate at restart or remapping. Independent physical fields can be reordered:
lifting reads the coarse parent solution and fixed wavelets; detail reconstruction
reads the even/even nodes established by lifting. Detail outputs do not feed the
other detail stencils.

`WAVETRISK_STAGE193_SCALAR_COLUMNS=1` selects the candidate. The default remains
0 pending cluster performance acceptance. The original scalar kernel is retained
as a control and executes only when the new kernel is disabled. Option agreement
is checked collectively at inverse entry, outside per-block callbacks.
`WAVETRISK_STAGE192_SPARSE_OUTER=0` is explicit throughout the Stage 193 campaign.

## Duplication audit and scope

Production native `WT_after_step` already skips the full legacy scalar/vector
wavelet and inverse computations. Full legacy tendency calls remain gated by
the validation oracle. Removing those oracle routines from the source would not
save production work. The duplicate work removed here is repeated preparation
inside the native scalar inverse, not a second complete legacy solve.

Other possible work remains, including inner-vector geometry repeated across
layers, full Domain wavelet zeroing, and Domain/block field staging. They are
not removed speculatively in this patch. The new measurements determine how much
scalar inverse preparation contributes before selecting the next change.

## Validation and measurement

The focused test extracts the actual old/new Fortran kernels and compares their
results at O0 and O2 with contraction off and bounds checking. It covers all three
storage classes, multiple scalars, frozen nodes, negative/nonphysical layer
offsets, empty atmospheric ranges and nonfinite-input rejection. It independently
counts accessor calls and verifies the work counters against the reference.
MPI tests cover unset/explicit options, malformed settings and rank disagreement.

The full local protocol uses J4/J6, four MPI ranks, 30 physical layers, resume=2,
time_end=0.33 and dt_write=0.02: checkpoint 3 is written and reloaded, followed by
further integration and vertical remapping. Results and precise evidence are
recorded in `parallel-blocks-stage193-local-results.json`. Both final O2 and O0/check builds passed all 12 steps, checkpoint 3
write/reload, and two vertical remaps. Seven exact/fixed-screen comparison gates
passed, including checked-versus-optimized and both builds against the accepted
baseline, saved oracle and legacy screen. Printed workloads and inverse transfer
counts matched. The final source bytes match the compiled build tree.

Across the complete profiled window, 102,935,040 scalar field samples were
produced. Resolved field addresses fell from a reference-equivalent 694,811,520
to 11,580,192; geometry reads fell from 720,545,280 to 12,009,088. Both preparation
counts fell 98.33%. These are operation counts in this kernel, not a claim of
98.33% lower runtime. Local timings are not accepted as performance evidence
under workstation memory pressure. The
initial optimized run with no Stage 193 selector was rejected as a stale build;
it is excluded from all Stage 193 comparisons.

The cluster kit uses the unchanged J5/J7 fixture and 83 ranks on exclusive nk.
It compares accepted Stage 191, Stage 193 and pinned legacy Domain with common
checkpoint I/O hardening. It keeps six rotating timing triples, checked runs
with oracles off/on, exact block/checked gates, fixed legacy screening, a disabled
candidate control, and separate profile runs. The profile report requires
identical inverse transfer calls/messages/bytes and reports:

- Complete matched wavelet transaction and inverse parent time.
- Combined inverse kernel time with equal call coverage. Stage 194 inspection
  established that this historical scalar/inner timer also encloses the
  outer-vector child; its label does not isolate scalar and inner-vector work.
- Compiled scalar stencils, produced field samples, field address resolutions
  and geometry reads, plus reference-equivalent counts for that same workload.

Profile counts include the checkpoint step. Nested times are not added together.
The legacy comparison remains a fixed numerical screen, not bitwise equivalence.
No cluster gain is claimed until that campaign passes and shows repeatable
ordinary-step improvement. See `parallel-blocks-cluster-193.md` for the complete
build/run/single-archive collection protocol.
