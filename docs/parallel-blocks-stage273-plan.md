# Stage273: isolated contiguous partition correction

There is good evidence for an assignment defect, not proof of a speedup.
Stage272's exact counters show that at 256 ranks, 160 Domain owners retain all
producer patch visits and 77.21% of block patch/divergence work. The other 96
ranks have 37.5% of the cores but only 22.79% of that measured block work. Both
observer timing calibrations failed, so phase durations do not supply a
quantitative speedup estimate. Untouched timing nevertheless worsened from
547.507 to 594.543 ms going from 160 to 256 ranks.

## One implementation change

Replace fixed-average, increasing-rank next-fit with a **minimax contiguous
partition**, followed by balanced feasible cuts. Catalogue order, leaf-weight
formula, block identities/subdivision, Domain placement, numerical kernels,
communication protocol and all payload proof reductions remain unchanged.
Only `src/parallel_block_mpi.f90` differs from accepted Stage266 `ff740744`.

The helper constructs 64-bit prefix sums and binary-searches the smallest load
capacity allowing at most `min(blocks,ranks)` contiguous groups. A greedy
feasibility scan is exact for nonnegative weights. A suffix feasibility table
then restricts each cut to choices that satisfy this capacity and leave enough
blocks/groups for the remaining ranks. Among feasible choices it selects the
cut closest to the cumulative equal-work target, preferring the earlier cut
on ties. Comparisons use integer quotient/remainder arithmetic, avoiding
floating threshold ambiguities and full-total-times-rank overflow.

Each active rank gets at least one block; excess ranks are empty. Empty input
returns immediately; all-zero weights balance block count. Negative weights
and invalid extents are rejected. Diagnostics use 64-bit load sums as well.
For positive production weights, feasibility search costs O(B log W), followed
by linear suffix and cut passes; workspace is O(B), where B is block count and
W is total weight. This is replicated catalogue work and belongs in the
complete-step benchmark; its added cost must not be ignored.

The peak-load constraint matters: unconstrained nearest cumulative cuts can
improve the distribution while making the heaviest rank worse. The delivered
algorithm instead guarantees the smallest possible maximum load for the
current contiguous weight model. That cannot exceed the maximum of any valid
contiguous legacy assignment on the same weights. It does not minimize runtime,
communication, variance, or retained Domain work, and is not a proof of speedup.

For 640 equal-weight blocks on 256 ranks, the old rule assigns three blocks
to 192 ranks and one to 64. The candidate assigns two to 128 and three to 128,
with the same optimal maximum of three and no underfilled one-block tail.
This is an algorithmic counterexample, not a reconstruction of actual J5 weights.

## Validation evidence and real catalogue checks

The test harness extracts and compiles the actual candidate helper and the
actual baseline assignment loop from the sealed source trees at O0/check and
O2. Candidate cases include nonuniform/heavy indivisible blocks, zero/empty
inputs, more ranks than blocks, equal blocks at 83/160/256/640 ranks, randomized
cases and totals above 32-bit range. An independent dynamic-programming and
exact-rational reference validates optimal capacity and cut choices. The old
Fortran loop is checked against the offline legacy reconstruction as well.

Candidate scientific validation runs enable the existing detailed-diagnostics
selector. Rank zero prints each actual catalogue's block ID, root Domain/patch,
leaf weight, final owner and geometric owner. Independent review checks complete
records, identities, split counts, candidate cuts and maximum load. It computes
legacy ownership on these same weights and reports per-rank loads/block counts,
Domain-owner/free weight shares, zero-block ranks and reassignment counts.

These records establish whether the real J5 weighted distribution improves.
They do not establish equal execution cost or measure communication overhead;
weights are still leaf-patch proxies. The offline old assignment is clearly
labeled a counterfactual. No source-level phase timers are added. Every timing
run explicitly disables catalogue diagnostics, and review rejects such output
in a timing run.

Local qualification uses two serial builds and sequential four-rank optimized,
checked and independent-oracle J4 replays against source/binary/physics-sealed
Stage266 controls. Local timing is not a cluster performance result. Retain all
exact checkpoint/VTK, state, restart/remap, boundary and oracle gates; do not
weaken scientific tolerances merely because ownership changed.

## Cluster campaign and acceptance

The kit is self-contained except for the existing sealed J5 fixture. One
exclusive nk/bb13 allocation makes four serial builds sharing 23 physics objects:
Stage266, legacy main `6a245db7` with established checkpoint I/O hardening,
optimized candidate and checked candidate. Placement uses distinct physical
cores, alternating sockets for counts above 16, and is checked for every run.

| Purpose | Ranks | Versions/order |
|---|---:|---|
| Validation | 16 | main, Stage266, candidate, checked, candidate oracle |
| Validation beyond Domain count | 256 | candidate oracle |
| Timing | 83 | main, Stage266, candidate, candidate, Stage266 |
| Timing | 160 | main, Stage266, candidate, candidate, Stage266 |
| Timing beyond Domain count | 256 | Stage266, candidate, candidate, Stage266 |

This is twenty solver runs. Two block runs per version/count support paired
comparisons for a potentially modest change; main has one supported-count anchor.
Use ordinary steps 2–7 and 10–12; checked/oracle runs are excluded from timing.
Report both pairs, spreads, full-step medians and scaling. Timing rows must
contain neither catalogue diagnostics nor added phase profiler output.

Accept only if science passes and uninstrumented performance supports a useful
benefit without repeatable regression at another count. Better proxy balance
alone is insufficient. Communication/locality changes, fixed per-block costs,
replicated partition overhead, hardware effects and retained Domain computation
can offset gains. No automatic repeat/promotion and no promised speedup percent.
Reject a regression before stacking another change.

Even a successful partition correction leaves the Domain producer ceiling.
The subsequent architectural work must move a complete producer/consumer chain
to block owners, preserving input epochs and retiring its corresponding transfers.
The goal remains efficient execution on at least 640 ranks for 160 J5 Domains,
with enough independent block granularity and eventual multi-node verification.

Budget **20–25 minutes after allocation**, based on Stage271's measured 21:47
for four builds/twenty runs. Validation catalogue logging adds some I/O; the
three-hour Slurm limit is a ceiling. Queue and transfer time are additional.
Runner UTC/monotonic elapsed time and Slurm accounting are collected automatically.
See the [no-edit commands](../test/parallel_block_profile/partition_273/README.md).

## Local qualification completed

Both serial builds and all three sequential four-rank J4 replays pass. Optimized,
checked and oracle cases each match all twelve step states, numerical records
and 62 VTK members exactly. The oracle passes 48 exact inverse transactions and
84 oracle installations. Each replay's fourteen actual catalogues matches the
independent reference and has no worse peak weight than the reconstructed old
assignment on those same weights.

The actual Fortran assignment harness passes 792 valid cases plus four invalid
input checks across O0/check and O2. All ten workflow tests pass. Fresh extraction
of the final archive verifies 46 sealed files, all three source families, imports,
shell syntax and disabled timing diagnostics; the delivered workflow and Fortran
tests pass again. The archive has no macOS extended-attribute headers.

The candidate remains unpromoted pending cluster performance. Machine-readable
hashes and evidence are in [the preparation record](parallel-blocks-stage273-preparation.json).
