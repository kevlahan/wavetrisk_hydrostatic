# Stage272: explain scaling beyond the Domain count

Stage271 measured accepted Stage266 at 549.270 ms/ordinary step on 160 ranks
and 584.620 ms on 256: 60% more ranks produced a 6.436% slowdown. The fused
proof-reduction experiment was rejected and solver source restored exactly to
Stage266 `ff740744608731a05d53771c32a3055309f3a586`. The diagnostic starts there.

The question is whether useful shared producer work remains concentrated on
geometric-Domain owners while extra block owners wait, or whether replicated
preparation, representation transfers or communication growth dominate instead.
Source inspection confirms retained Domain traversal; it does not quantify
which mechanism accounts for the observed slowdown.

## Instrumentation and evidence

A sealed diagnostic copy retains Stage270's bounded nested phase observer and
adds seven slots: ownership census, patch inventory, Domain basic operators,
Domain primitive preparation, native-mass MPI wait, Domain boundary MPI wait,
and the local Domain producer patch batch. Existing producer, capture, restriction, kernel,
preparation, inverse and migration scopes provide context without another main
profiling campaign. The unchanged production build contains no observer calls.

Before each of four RK producer stages, record `size(grid)`, actual final-owner
block count, Domain patch inventory over active levels, and block patch inventory.
The census has its own measured scope. Domain ownership is stable within the RK
transaction; analysis rejects fractional stage-averaged ownership rather than
misclassifying mixed records. Global Domain totals must equal 160 at every
ordinary step. Groups are formed from actual rank/step records. Group membership
is retained, including an empty Domain-free group at 160 if appropriate.

Time `basic_operators` at its existing level/layer transaction boundary; count
producer patch-layer visits and Domain-level-layer visits using patch-list lengths
once per Domain. Separate the first local producer patch batch from the enclosing
basic-operator transaction, and measure native-mass and Domain-boundary MPI waits
explicitly so those waits are not labeled useful producer computation. Record native scalar restriction's existing parent-patch,
restricted-edge and divergence-replay counters. Inventories and executed work
have distinct labels and must not be equated. Capture publication records actual
remote logical doubles selected from the executed full/deduplicated/compact
transport branch, excluding self copies. Capture, restriction, temperature-edge
and inverse send/receive totals must balance globally in each ordinary step.
Existing migration byte counts remain available.

The observer adds no MPI call, barrier or timed-step file output. It records
bounded transactions and integer batch work, not clocks inside per-node physics
callbacks. Per-rank files are written after the solver completes. Raw files are
hashed, collected and independently reparsed locally; reported nested wall/CPU
self-time conservation and scientific results are checked.

Each owner group reports minimum, median, p90, maximum and mean phase duration,
call count and work count, retaining each observer run separately. Statistics
pool rank/step samples within a run, not independent repetitions. Mean-rank
self-time changes and the existing longest-local-step rank composition provide
context; independent region maxima must not be added into an invented critical
path. MPI scopes combine late-arrival waiting and data transfer. CPU includes
spinning and is not evidence of useful numerical work. Inclusive times overlap.

## Campaign and calibration

One exclusive nk/bb13 job, physical-core placement verified before and during
every run, one hardware thread per physical core. Four serial builds share the
same 23 Simple-physics objects: accepted Stage266, legacy main `6a245db7` with
established checkpoint I/O hardening for correctness only, optimized observer,
and checked observer.

Six validations: main, Stage266, observer, checked observer and observer oracle
at 16 ranks, then observer oracle at 256. All retain established checkpoint,
trajectory, restart/remap, storage and 62-member same-count VTK gates. Block
checkpoints remain exact across counts; both oracles retain 48 inverse proofs
and 84 installations. There is no newly optimized production candidate.

Eight timing runs use this fixed sequence:

| Run | Ranks | Version |
|---|---:|---|
| 1 | 160 | untouched Stage266 |
| 2 | 160 | observer |
| 3 | 256 | untouched Stage266 |
| 4 | 256 | observer |
| 5 | 256 | observer |
| 6 | 256 | untouched Stage266 |
| 7 | 160 | observer |
| 8 | 160 | untouched Stage266 |

This yields A-B-B-A independently at each count and brackets the 256 group with
160 runs. It reduces rank-group order bias but does not eliminate hardware or
run-to-run variation. Use ordinary steps 2–7 and 10–12 from the unchanged J5
fixture; oracle and checked timings are excluded. No automatic repeats.

Apply the established observer calibration separately at each count: absolute
mean shift <= max(3 ms, 1% of control), each paired shift <= max(5 ms, 2%), and
within-version spread <= 2%. These are diagnostic resolution gates, not a
statistical confidence interval or production performance acceptance threshold.
If either fails, retain all results but mark cross-count phase attribution
unqualified. Do not transfer the prior 83-rank calibration to these counts.

## How findings guide the next implementation

If producer computation remains concentrated on Domain owners while Domain-free
ranks wait at capture/restriction, select one complete producer-to-block-consumer
chain to relocate to block owners, removing its corresponding transfer/duplicate
work while preserving Stage267 input epochs, boundaries and exact oracles.
If preparation or communication growth dominates instead, address that specific
cost first. Do not label all producer time removable or promise a recoverable
speedup from phase time alone. Hardware bandwidth/NUMA/frequency remain alternate
contributors requiring further evidence if the software pattern is inconclusive.

Any resulting implementation must pass an uninstrumented comparison against
accepted Stage266 and supported legacy main anchors before acceptance, including
83/160/256 ranks. The goal remains efficient execution on at least 640 ranks for
160 Domains; this is an intermediate single-node diagnosis, not completion.

## Local validation and runtime estimate

Serial optimized/checked builds and sequential four-rank J4 optimized, checked
and independent-oracle replays qualify the observer overlay. Clock-stack tests
cover disabled/invalid selectors, mismatched scopes and capacity failures at
O0/O2. Analysis tests exercise nontrivial owner rank identities, missing census,
impossible Domain-free work and remote payload corruption. Workflow tests cover
all four builds/fourteen runs, each calibration, launch/duplicate protection,
purged/failed collection, independent review and 256-core placement.

Stage271 took 21:47 for four builds/twenty runs; Stage270 took 15:38 for five
builds/fifteen runs. Estimate **16–22 minutes after allocation** for this four-build,
fourteen-run campaign, subject to the 256-rank oracle and filesystem conditions.
Queue and collection time are additional. Both runner monotonic elapsed time
and raw/parsed Slurm Start/End/Elapsed are collected automatically.

Use the [no-edit commands](../test/parallel_block_profile/scaling_272/README.md).

Local qualification is complete: both final serial builds and all three final
four-rank replays pass exact checkpoint, step-state and 62-member VTK comparisons.
The oracle preserves 48 exact inverse-return proofs and 84 installations. Every
rank has all twelve steps and 67 observer slots, with conserved nested timing,
actual ownership totals and selected remote payload totals. Thirteen workflow,
analysis and observer tests pass, including twelve compiled clock-stack cases.
Fresh delivery extraction verifies all 51 sealed files, source manifests, full
helper imports, shell syntax, selectors and absence of macOS xattr headers; all
thirteen delivered tests pass again. Production source still matches Stage266.

An initial local MPI launch was blocked by sandbox socket restrictions before
solver execution. Completed builds were reused with local socket access. The
initial observer draft passed three replays; after separating native-mass and
Domain-boundary waits from the producer, both builds and all three replays were
repeated for the final overlay. Only the final overlay is delivered. No build
or second MPI replay ran alongside an active replay.

The [sealed preparation record](parallel-blocks-stage272-preparation.json)
records final evidence and archive identity. Cluster calibration/scaling remains
pending; local timing is not used to claim a performance improvement.
