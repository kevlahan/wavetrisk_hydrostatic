# Proposed Stage 197: reduce intermediate ownership transitions

Current direction: the user has required diagnosis of the main legacy/block cost causes before any further optimization. The [diagnosis objective and evidence requirements](parallel-blocks-performance-gap-diagnosis.md) supersede optimization recommendations below. Existing results are retained as historical evidence.

Status: original proposal. The implemented first ownership change and its
remaining scope limits are documented in parallel-blocks-stage197.md. Stage 196 is committed as
`9bc3a19556dab79b3cba3944cf5df53517389da3` on `codex/stage187-validated`.
Its numerical source is the previous production implementation; the stage
adds validated measurements and records the missing-counter limitation.

## Primary target: the complete WT ownership lifecycle

The four dynamics WT transactions contribute 0.280642 s of the measured mean
excess; 0.205488 s is in boundary/ownership transfer routines. Inverse transfer
planning/packing/install costs 0.055320 s and inverse posting/self-copy/wait
0.091748 s on blocks. This makes the full producer/consumer ownership chain
a better primary target than another isolated arithmetic expression or a
message-byte reduction. These times are locations, not recoverable-speedup
predictions; required dependencies and late peer computation can cause waiting.

The proposed production change is to retain native intermediate coefficients
through transform, compression and inverse reconstruction, with publication
only where a Domain consumer actually requires it. Start with the explicit
uncompressed writeback -> full import -> compression -> compressed writeback
sequence. Extend the lifetime of the inverse staging buffers so compatible
consumers reuse installed data rather than routing it again.

Relevant implementation points are in `src/parallel_block_mpi.f90`:
`activate_complete_block_wavelet_output`,
`prepare_block_native_wavelet_compression`,
`activate_block_native_wavelet_compression`, and
`prepare_block_native_inverse_boundaries`; inverse routes and buffers are in
`src/parallel_block_inverse.f90`.

Before removing a transfer, map its producer, authoritative owner, field
family, level, boundary refresh and next consumer. The source currently places
scalar coefficients in compact blocks but stages vector coefficients on Domain
owners. A shortcut must supply the complete required family, not assume both
families already reside together. Preserve the checkpoint snapshot, ADJZONE
mask application, geometric boundary aliases, provisional RK state and final
adaptation consumers. Invalidate retained data at grid-generation, remap,
restart and relevant field updates. Keep per-value arithmetic order unchanged.

Use one reversible candidate selector and compare against untouched Stage 196
and pinned legacy. Keep the production change limited to this ownership path;
measure the full WT and step to detect costs moved to another routine. The
selected transfers' fraction of the 0.205488 s envelope must be measured;
eliminating one small round trip alone is not success at closing the entire gap.

## Supporting measurement changes

Repair counters 4, 5, 12–15 at their call sites so they operate independently
of the disabled old profiler. Add coverage checks conditioned on actual work,
not a universal requirement that every rank has a positive count.

Record actual request-bearing exchanges by dependency phase and level, values
packed/installed, field family and buffer reuse. Retain rank distributions and
local post-to-wait intervals. Do not add barriers or assume independent rank
clock epochs are synchronized. Separate no-request completion calls from MPI
waits. This tests whether removing ownership transitions reduces traversal and
serialized exchange phases, rather than merely shrinking payloads.

Also separate matching scalar tasks on both implementations: direct flux,
restriction, divergence and temperature-boundary work. Count visited versus
computed cells/edges/fields and distinguish setup/lookup from numerical loops.
Stage 196's 0.082273 s scalar-group difference has unequal scope boundaries and
cannot be treated as an arithmetic-only slowdown. Do not combine a separate
scalar algorithm rewrite with the ownership change in the first candidate.

Mass exchanges by vertical layer are a subsequent independent experiment.
Their 0.074860 s wait cost warrants attention, but the present layer-major loop
shares flux/scratch storage with velocity and diffusion consumers. Batching
layers requires an explicit storage/dependency redesign, not simply posting
all existing requests earlier.

## Short cluster campaign: 14 executions initially

First run local J4/J6 correctness and checked restart tests. Keep cluster
J5/J7, 83 ranks, the same nk placement, compiler flags and physics objects.
Retain the existing short checkpoint-3 fixture and end time: it exercises
checkpoint 4, reload and post-restart work. Shorten repetition rather than
remove the restart/adaptation coverage.

| Purpose | Runs |
|---|---:|
| Checked candidate, validation oracles off and on | 2 |
| Untimed warmup/reference: Stage 196, candidate, legacy | 3 |
| Two ordinary timing rounds, all three implementations | 6 |
| One focused profile per implementation | 3 |
| Total | 14 |

The 14-run schedule requires a new runner; it is not an option already wired
into the frozen Stage 196 kit. Run untouched references before comparisons
that require them. Reverse the implementation order in the second timing
round. Baseline and legacy profile executables use the same new observer as
the candidate; each must match its untouched reference. Keep the independent
reference builds and common build provenance even though there are fewer runs.

Every run still checks completion, matching workload, restart/remap coverage,
source/binary identities and affinity. Require candidate/Stage 196 exact
checkpoint equality, checked/optimized agreement and the existing fixed legacy
screen. Profiled executions must match their reference, have valid operation
coverage and conserve the common timing envelopes. Keep swap/RSS records.

Use ordinary-step wall time, complete WT/outside-WT wall time and targeted work
counts to decide whether the mechanism helps. Report both rounds and their
spread, not a precise confidence interval from two samples. A decrease in MPI
wait alone is insufficient if computation or another phase rises correspondingly.
Inspect profile-versus-bare overhead; one profile per implementation is a
screen, not a calibrated statistical upper bound. Add focused observer checks
if the new hooks perturb the result materially.

For a clear repeatable gain and passing correctness, this is enough to choose
the next development direction. If the gain is small, inconsistent or close to
the observed run spread, add pairs before claiming an improvement. A final
performance claim can receive three or more pairs and broader workload testing
when warranted. There is no need to repeat all four observer modes six times
for each small development iteration.

Process CPU time is not the principal performance criterion: MPI spinning can
consume it. Wall-clock time and causal work counts answer the development
question more directly. Keep collecting results and build/run/profile records
as one inventory-bearing `stage197-review-JOB.tar.gz`, transferred in one scp
command to `murier.math.mcmaster.ca:~/Desktop/` after the new runner is prepared.
