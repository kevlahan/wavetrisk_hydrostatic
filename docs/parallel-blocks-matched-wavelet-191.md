# Resolve the Stage 191 inverse/boundary profiling ambiguity

We can now identify a major cost center using the existing cluster data. The
complete `WT_after_step` transaction is 4.33x slower in Stage 191 than legacy.
Its measured excess accounts for 48.7% of the whole-window profiled timestep
excess. Inside that transaction, the native inverse spends more time gathering,
scattering and exchanging aliases between ownership layouts than in arithmetic.
This is stronger evidence than the previous comparison of exclusive timer rows.
It does not identify every remaining cost or prove that all transfer time can
be removed.

## Use the complete transaction, with the same call count

Both implementations bracket `WT_after_step` from before initialization of its
wavelet work through completion of compression and inverse reconstruction. The
native and legacy callbacks differ internally, but the entry/exit contract is
the same: compute wavelets and reconstruct the adaptive-grid stage state.
The inclusive timer includes all nested boundary and MPI work regardless of
which region name accounts for its exclusive time.

All three profiles record 60 calls per rank across the same 12 timesteps, with
identical printed workload and passing checkpoint gates. That establishes a
matched parent scope. By contrast, bare native inverse and legacy inverse timers
record 48 and 73 calls/rank respectively; they are not comparable transactions.

| Complete WT_after_step | Stage 190 | Stage 191 | Legacy |
| --- | ---: | ---: | ---: |
| Calls per rank | 60 | 60 | 60 |
| Rank-mean inclusive seconds, whole window | 5.61379 | 5.60228 | 1.29358 |
| Seconds per step, including checkpoint step | 0.46782 | 0.46686 | 0.10780 |

Stage 191's wavelet-transaction excess is 4.30870 s over the window, or 0.35906 s
per step including the checkpoint step. The corresponding whole-profile step
excess is 8.84604 s, so this transaction accounts for 48.7% of that excess.
These are consistently scoped profiled means. Do not divide 0.35906 by the
unprofiled ordinary-step excess of 0.783 s or claim exactly the same fraction
of the production critical path. Existing records lack per-step inclusive
parent times, so the checkpoint step cannot be removed from this comparison.
The installation change barely affects this transaction, as the Stage 190 and
191 columns demonstrate.

## What the native inverse actually spends time doing

Within the Stage 191 WT parent, the native inverse timer is 3.59347 s for 48
calls per rank. Its dedicated gather/scatter/alias helpers report:

| Native inverse helper | Rank-mean seconds | Aggregate remote payload |
| --- | ---: | ---: |
| Gather | 0.43163 | 6.36 GB |
| Alias exchange | 0.83682 | 10.56 GB |
| Scatter | 0.77543 | 25.35 GB |
| Total | 2.04388 | 42.26 GB |

The helper total is about 56.9% of native inverse inclusive time and covers
496,004 aggregate messages. It includes packing, communication waiting, local
copies and installation. These helper times are contained within the inverse
and WT parent times: never add them to either parent. The payload counters count
remote sends, excluding self copies; they are not total memory traffic or a
measurement of NIC traffic on this single-node allocation.

The existing inner/outer kernel inclusive timer is 0.76038 s. Explicit inverse
wait self time is 1.41249 s, with process CPU time 1.40342 s; MPI spinning can
explain that CPU use, so it must not be described as useful computation.
Individual phase maxima are not summed into an invented critical path.

Source inspection explains the transfers. Regular scalar and inner-vector
kernels execute on final-owner blocks, while sparse outer-edge production and
geometric aliases retain their established route owners. `transfer_nodes`
therefore gathers block values into the inverse dependency workspace and scatters
results or boundaries back. `synchronize_block_inverse_scalars` and
`synchronize_block_inverse_vectors` repeat the appropriate transfers at levels;
provisional transforms also complete an all-level boundary contract and publish
state for Domain consumers. This is retained ownership/layout interoperability,
not evidence that the complete legacy inverse is being recomputed alongside
every native inverse.

In particular, `WT_after_step` skips the legacy forward-wavelet loops when
`native_output` is true and takes the native inverse branch. The number of full
WT calls is equal to legacy. The words “validate” in several callback names do
not imply duplicate production evaluation; their production branches generate
the native result. The separate complete legacy tendency remains oracle-only.

## Rank imbalance still matters

Over the 11 ordinary candidate steps, rank 52 averages 1.099 s outside the named
MPI regions and 0.235 s in them; rank 82 averages 0.472 s outside and 0.862 s
inside. Their total step times are both about 1.334 s. This is strong evidence
of unequal arrival/work distribution. “Outside named MPI” is not synonymous
with pure arithmetic, and the remaining waits cannot all be blamed on network
latency. Ownership reconciliation can add local work as well as communication.

The scalar/tendency ambiguity is not fully resolved by a single existing row:
legacy's shared-tendency parent includes its scalar work, while the block code
performs native tendency assembly after the retained compatibility producer
returns. Comparing only those shared-parent rows would omit block work. A future
profile should bracket the complete RK-stage input-to-RHS transaction in both
implementations and record per-step parent times and bytes, without inserting
new barriers. That remaining instrumentation need does not weaken the matched
WT result above.

## Recommended next experiment

Target inverse ownership reconciliation before another small arithmetic change:
retain needed boundary state on final block owners across levels, or move the
sparse outer-edge/alias producer-consumer chain to those owners. First narrow
redundant gather/scatter routes while preserving geometric write order, consumer
aliases and exact numerical gates. Keep the old route as a control and compare
the complete WT transaction against the accepted installation-only baseline and
legacy. Count remote/local bytes and route visits, and verify that the reduction
survives the whole timestep. The 42 GB payload total is evidence of work, not a
promise that all of it is avoidable.

Storage reuse remains worthwhile, but the remaining construction/migration group
is only about 0.101 s per ordinary step; it is no longer the strongest first
candidate for closing this gap.

## Reproduction and commit boundaries

- `37ac5364`: exact original cluster-tested Stage 191 implementation.
- `149b519f`: cluster evidence, test/profiling infrastructure and log collection.
- `352e803c`: installation on / experimental mass overlap off by default, with
  selector tests and a fresh local exact restart run. Three installation-only
  cluster ablations already exercised the selected configuration.

Run `analyze_matched_wavelet_191.py /path/to/extracted/stage191-review-54561`.
It audits the campaign first against the pinned tested source commit, checks
matched parent call counts and complete rank coverage, and reconciles transfer
call totals with the detailed wait counter. Output is preserved in
`parallel-blocks-matched-wavelet-191-analysis.json`. No new performance gain is
claimed by this diagnosis, and no solver instrumentation was changed for it.
