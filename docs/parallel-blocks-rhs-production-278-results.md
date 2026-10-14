Current reviewed result: [kit-v4 job55501](parallel-blocks-stage278-job55501-results.md) passes science/ownership and improves traffic/time, but still fails absolute performance, memory and scaling acceptance. Stage266 remains accepted.

# Stage278 connected RHS/RK: ready for the bbserv campaign

Current decision2026-10-13: [job55500 review](parallel-blocks-stage278-results.md)
verifies science and real zero-Domain RHS/RK ownership, but finds46–56% slower
ordinary steps than266. Candidate remains unpromoted;266 stays accepted. The
preparation and delivery record below is historical.

Updated 2026-10-09 after job 55499. The corrected [kit-v3 archive](/Users/kevlahan/Desktop/rk278-connected-rhs-kit-v3.tar.gz)
is sealed, uploaded and extracted on bbserv; remote source/fixture preflight
passed. The replacement has not been submitted. Accepted
solver remains Stage266 `ff740744608731a05d53771c32a3055309f3a586`.

Job55499 failed during the checked build after 3:02 (Slurm start 17:58:13,
end 18:01:15), following three successful optimized builds and before any solver
run. GCC13.2 treated automatic allocation at `lengths=transpose(resolved)` as
possibly uninitialized under `-Werror`. The local GCC15 build had not reported it.
The correction explicitly allocates the transposed extent and assigns to the
array section. All checked flags, numerical operations and acceptance gates stay
in force. Original job evidence is preserved in a 17-file verified archive.

The same GCC13 compiler reproduces the original warning and compiles the fixed
module; the remaining 21 changed Fortran files also compile to objects with the
same strict flags, using existing module interfaces. Complete local checked,
optimized and oracle J4 replays all pass again (12 steps, exact checkpoint, 62 exact VTK members each).
The full J5 checked build/run remains a cluster gate. Controls and campaign
helpers are unchanged; the candidate differs from kit-v2 in exactly one source
file. Four fresh serial builds retain the established campaign workflow.

Archive SHA256: `e72cd9e60334c709e249f3569dbe20ecd8b2af436cde06c2fdb3a1fdcfbb1091`. Size: 19,987,701 bytes.
The seal verifies 50 outer kit files and freshly extracts all three source
families, including 251 candidate build-source files and the explicit 23-file
overlay. Stage273 partition changes and the Stage277 cache remain excluded.

Recovery evidence: [evaluation](../.stage278-validation/gcc13-recovery-v3/evaluation.json),
[GCC13 checks](../.stage278-validation/gcc13-recovery-v3/compiler-check/results.json),
[GCC13 object builds](../.stage278-validation/gcc13-recovery-v3/compiler-codegen/results.json),
[failed job archive](../.stage278-validation/gcc13-recovery-v3/rk278-review-55499.tar.gz),
[raw replay manifest](../.stage278-validation/gcc13-recovery-v3/files-sha256.json),
[helper tests](../.stage278-validation/gcc13-recovery-v3/campaign-tests.json),
[seal](../.stage278-validation/gcc13-recovery-v3/seal.json).

## Implemented ownership and paired retirement

The opt-in selector now reads catalogue interiors from resident block state,
executes the connected pressure/Bernoulli/Exner, shared flux/PV/KE/vorticity,
scalar/velocity diffusion, restriction, source and gradient chain on target block
owners, and retains authoritative tendencies in the real resident RK accumulator.
It reuses the direct kernels and ownership/alias evidence developed after266.
There is no ordinary-stage full primitive/reference return stream.

The selected caller bypasses the old Domain compatibility RHS, scalar tendency
capture and numerical record transport, velocity tendency/remainder export,
full temperature-stage publication, and displaced interior mass/velocity RK.
The retired RHS/capture entries fail if reached; existing legacy-call and
writeback-count checks remain active. The actual RK consumer records48 retained
transactions per rank, not just producer assignments.

Required boundary/scaffold state and named halo/coarse RK readers remain bounded
adapters. The separate split Simple-physics transaction, inverse/lifecycle
canonical publication, old native geometry compilers and allocated caches remain.
No savings are claimed for those retained paths or allocations. The geometry-only
wavelet/inverse adapter supplies metadata without reinstating scalar RHS capture.
The accepted `(raw_rhs + Exner_term) - Exner_term` velocity rounding remains on
resident owners; raw mass/velocity RK aliases keep the accepted boundary epoch.

This advances independent RHS/RK ownership and removes displaced work in the same
slice. Whole-timestep independence, remaining publication/physics migration,
finer subdivision and useful P>=4D scaling remain separate scorecard milestones.
The delivered immutable fixtures exercise diffusion orders(2,2,2); orders0/1 have
component tests, not complete active-production qualification. Scope remains
restarted compressible non-mode-split RK4 climate, Simple physics, no NCAR/SSO.

## Qualification and preserved failures

All three final source-bound four-rank J4 replays pass: checked, optimized, and
optimized with independent dynamics/adaptation/remap oracles. Each completes12
steps,48 RHS/RK stages,30 physical layers,12 topology/mask plan epochs, both reload
and remap events, an exact accepted266 checkpoint and all 62 exact VTK members.
No numerical tolerance was relaxed. The first full active replay failed at about
1.47e-9 temperature-wavelet error; preserving accepted velocity cancellation and
raw RK-reader rounding resolved it. The1138-file development archive preserves
original failures, the deliberately stopped diagnostic run, and corrected proof.

Final exchange tests pass in optimized and checked modes:34 MPI invocations in
total, one/four ranks, identity-copy pruning and in-place sign changes, provider
demand, large rendezvous messages, empty participants and14 expected protocol
rejections per compiler mode. Twelve campaign-helper tests, shell syntax, complete
packaged analysis imports, immutable grid/source checks and fresh extraction pass.
Local limits remained serial builds and at most four MPI ranks.

Original kit-v2 evidence (preserved): [production evaluation](../.stage278-validation/production-evaluation.json),
[raw replay manifest](../.stage278-validation/production-replay/files-sha256.json),
[exchange results](../.stage278-validation/exchange-final/results.json),
[helper tests](../.stage278-validation/campaign-final-tests.json),
[seal](../.stage278-validation/production-seal.json).

## Quantitative performance case

The historical kit-v2 bare local pair used the same fixture, selectors,23 identical external
physics objects and contraction-disabled optimized builds. Its ordinary-step
medians are 4.823805s for266 and 4.842877s for278: candidate+0.3954%, speedup0.9961.
This single local pair shows near parity, not a speedup or a calibrated cluster
forecast. It does not justify repeating local timing to decide J5 cluster scaling.

The kit-v2 replacement ledger records4,472,284,800 remote payload bytes plus
6,266,112 epoch-header bytes in783,264 stage messages across the entire four-rank
J4 optimized replay. These totals include input, derived and RK streams; topology
metadata is additional. Resident input gathers comprise386,709,120 values and
explicit uncatalogued adapters132,021,960 values. Bounded raw RK readers account
for416,232,000 values; resident RK deliveries119,433,600 values. These are complete
local trajectory counts and must not be compared with Stage275 sampled per-layer
budgets as if the denominators matched.

The sum of per-component allocation maxima ranges358.7–437.3MiB per rank for
numerical buffers, routing payloads and rounding columns. This is not peak RSS
and excludes descriptors, compilation peaks, retained old workspaces, geometry
adapters and MPI storage. Primitive stencil geometry remains per target/event;
no reduced geometry-memory or cache-speed claim is made. The largest rank records
18.856s backend preparation,14.621s RHS evaluation and1.679s RK-reader preparation
across the complete optimized replay. These narrow timers omit other necessary
work; full trajectory timing and actual child peak RSS are collected on bbserv.

The [prior quantitative case](parallel-blocks-rhs-performance-case-278.md) remains
conditional: `saving_ms = 49.196 R - 193.987 d - N`. R is the removed fraction of
the measured mixed capture/remainder pool, d the producer cost change, N additional
costs. R=0.5,d=0,N=10 gives14.60ms; d=0.1 instead gives-4.80ms. Ownership redistribution
could reduce producer time, but its seconds are uncalibrated. The local ledger and
near-parity pair justify a bounded matched-control campaign; they do not establish
R,d,N at83/160/256 or promise closure of the39.80%/54.72% historical legacy gaps.
The new campaign must demonstrate that retirement/placement gains exceed all
replacement costs. Repeatable regression blocks promotion and further stacking.

## Campaign and dual-goal decision

The fixed initial nk/bb13 campaign uses one exclusive256-task node, mem=0, a60-minute
wall cap, four serial builds and sixteen serial runs: eight validations, eight bare
timings. It compares main/266/278 at83 and160, then266/278 at256. Validation includes
16-rank checked/oracle runs and a256-rank oracle/ownership run. All96 zero-Domain
ranks must execute authoritative RHS/RK work when they own blocks; local tests do
not establish that cluster gate.

The independent review recomputes science, ownership and traffic evidence, ordinary
step medians, complete trajectory time, rank peak RSS/faults, the remaining main gap,
160-to-256 speedup and relative efficiency. No automatic promotion occurs. The
initial kit does not repeat automatically; decision-changing ambiguity permits a
separately reviewed reversed group, at most eight extra bare runs. The conditional
four-node bb readiness leg follows a resolved Leg A; it is intentionally outside
this first kit. Do not merge nk and bb measurements into one scaling curve.

256 is intermediate for J5 D=160. Useful >=640-rank execution, finer granularity
and a user-agreed efficiency target remain open. The campaign tests both ownership
and performance rather than treating a zero-Domain correctness pass as scaling.
Planning allowance is25–40 minutes, provisional from the earlier18:45 schedule and
new local build/replay work;60 minutes is a resource cap, not a runtime guarantee.
Real Slurm start/end/elapsed and the exact run schedule are collected automatically.

## Copyable user commands

Kit-v3 is already uploaded and extracted on bbserv, and a copy is on the Mac
Desktop. The failed job archive has already been collected. In the established
bbserv compiler/MPI/Python environment, launch with:

```sh
bash "$HOME/rk278-connected-rhs-kit-v3/launch.sh"
```

Check progress:

```sh
bash "$HOME/rk278-connected-rhs-kit-v3/status.sh"
```

After the job ends, including a failed job:

```sh
bash "$HOME/rk278-connected-rhs-kit-v3/collect.sh"
```

Collection verifies the review archive and transfers it to the established murier
Desktop destination. No new user guidance is required to complete preparation.
Earlier stops were incomplete implementation checkpoints, not a need for further
authorization. Preparation is now complete; cluster results remain to be measured.
