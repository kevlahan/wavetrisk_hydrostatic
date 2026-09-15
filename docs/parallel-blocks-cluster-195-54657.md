# Stage 195 cluster review — job 54657

Stage 195 establishes a repeatable, matched accounting of the legacy/block
performance difference. The unchanged Stage 194 block implementation remains
2.262× slower by the ordinary-step median: 1.260 s versus 0.557 s.
The coarse observers have negligible observed effect at this campaign's
resolution. Almost all measured excess is in dynamics and adaptation/reconstruction.
Its division into redundant mathematics, representation work, and waiting is
still incomplete. No new production optimization or speedup is claimed.

## Evidence and validation

The supplied archive is `stage195-review-54657.tar.gz`, SHA256
`77e5b7f7f03d79121fd702df63030f396d137db12f7915af28dba94a04270767`.
The audit verifies all 14,102 inventoried files, five build-log digests and source
identities against the frozen Stage 195 kit. The baseline is Stage 194 commit
`adb2e91e0d1c16dd21a9bf1597816e42891da719`; the common observer source hashes
match in the block and legacy diagnostic builds. Numerical source in the local
repository still matches the untouched block baseline.

All 58 runs completed. The 146 supplied gates comprise 60 exact checkpoint
comparisons, 58 workload checks and 28 fixed-tolerance legacy comparisons.
The exact checks include observer-mode/reference equivalence and
checked-versus-optimized comparisons. Legacy equivalence uses the existing
fixed limits, not bitwise equality. Maximum reported absolute differences are
3.724e-12 velocity, 2.502e-12 mass and 7.859e-10 temperature (rounded upwards).
Checkpoint payloads and binaries are absent from the log archive: the audit
revalidates their recorded reports and identities, not the original arrays.

The common workflow was reparsed from all 3,652 supplied rank files. Phase
coverage, conservation, nested totals and detailed event intervals pass. All
runs have 12 steps, restart 3, save/reload 4, and post-restart remapping.
Ordinary accounting excludes checkpoint step 8 and includes remap steps 5 and 10.
All optimized runs have validation oracles disabled; checked-on enables them.

Builds use GNU Fortran 13.2.0, J5, RK4, contraction off, identical recorded
external physics objects; optimized builds use O2 and march=native, checked O0
with runtime checks. All ranks ran on bb13 with the same 83 distinct core
bindings (0–41 and 128–168), single-threaded numerical libraries and no recorded
swap-in/out change. Summed per-rank peak RSS medians are 21.693 GiB for blocks
and 11.029 GiB for legacy. These are not concurrent node RSS measurements and
do not establish a bandwidth bottleneck.

## Profiling calibration and authoritative performance

Six interleaved repeats per implementation and observation mode give:

| Mode | Blocks, s/step | Legacy, s/step | Statistic |
|---|---:|---:|---|
| Untouched executables | 1.260 | 0.557 | Median of printed ordinary-step medians |
| Minimal observer | 1.266636 | 0.557489 | Median of common maximum-rank wall medians |
| Coarse observer | 1.265504 | 0.558089 | Same common wall statistic |
| Detailed observer | 1.312418 | 0.626061 | Same common wall statistic |

Coarse versus minimal differences are -0.089% blocks and +0.108% legacy.
These are observed differences across repeated runs, not a formal upper bound
on observer overhead. They support using coarse phase accounting here.
Detailed observation changes those medians by +3.614% and +12.300%, respectively.
It remains unsuitable for an uncorrected legacy/block performance ratio or a
single uniform overhead correction.

Bare elapsed medians are 30.980945 s blocks and 22.554070 s legacy. All six bare
pairs favor legacy; ordinary-step median excess ranges from 0.693 to 0.706 s.
This is an attribution experiment on unchanged numerical implementations.

## Matched accounting of the ordinary-step excess

Use arithmetic means to make the phase ledger additive: mean over the same 11
ordinary steps, 83 ranks and six repeats. Blocks take 1.207112 s and legacy
0.535930 s, giving 0.671183 s excess. These means differ from the 0.703 s
headline median excess because the statistics and observation modes differ;
the two faster post-restart steps affect the mean more strongly.

| Matched part of the workflow | Blocks minus legacy, s/ordinary step | Share of mean excess |
|---|---:|---:|
| Dynamics: complete WT_after_step transactions | +0.281022 | 41.87% |
| Dynamics outside those transactions | +0.276661 | 41.22% |
| Adaptation and reconstruction | +0.115829 | 17.26% |
| All other outer phases, signed net | -0.002330 | -0.35% |
| Total | +0.671183 | 100% |

The six paired mean excesses range only from 0.668677 to 0.672195 s.
Phase differences conserve each paired excess to floating-point rounding.
Using the actual slowest rank once per step gives 0.671841 s excess and the
same leading outer phases. This does not exclude waiting: ranks can finish
close together because they synchronize. Independently selected phase maxima
are never added.

The full advance phase is 1.009242 s blocks versus 0.451559 s legacy, excess
0.557683 s. Adaptation/reconstruction is 0.178052 versus 0.062222 s, excess
0.115829 s. Physics contributes only +0.000257 s. Remap and ordinary output
contribute -0.002213 and -0.000374 s net. With-remap and without-remap steps
are separated in the JSON: the remap phase itself has positive excess on
actual remap steps, despite its negative all-ordinary average.

Checkpoint steps are accounted separately: mean excess 0.833569 s, including
only 0.070234 s in the output/restart phase. Ordinary slowdown is not caused by
checkpoint I/O.

## What the matched wavelet transactions establish

Both implementations make four dynamics WT_after_step calls plus one
adaptation-phase call per rank per step: 4,980 total rank-summed calls over
12 steps. Detailed local event traces validate these parent phases and the
non-overlap of wavelet and tendency transactions.

The four dynamics calls together take 0.353759 s blocks versus 0.072736 s
legacy, a 4.864× ratio. Their 0.281022 s excess accounts for roughly half the
dynamics difference. These envelopes contain the full transform, compression,
reconstruction, boundary and representation work; they are not pure arithmetic
timers. The adaptation-phase WT calls take 0.026461 versus 0.028403 s, slightly
favoring blocks. This distinction matters: the excess is associated with the
native dynamics transaction, not every invocation of the WT routine.

The separate tendency tags take 0.616206 versus 0.373093 s. Their boundaries
are implementation-specific, so their 0.243113 s difference is not a matched
estimate of redundant tendency computation. Keep it as a diagnostic clue
inside the broader matched dynamics envelope.

## Production-source evidence and unresolved causes

The source rules out treating the block version as a full duplicate legacy
step. `prepare_block_multistage_tendency` calls the complete legacy tendency
consumer only under `validate_oracle`. The compatibility producer retains
shared Domain primitive/physics inputs, but legacy scalar tendency and
velocity source/gradient comparisons are conditional on validation. Native
WT output skips the legacy atmospheric wavelet loops and legacy inverse
branch. Recorded production oracle fields and oracle self time are zero.
A routine or timer containing the word "replay" does not prove the same
numerical result is also computed by legacy production code.

There is concrete evidence of extra representation work:

* `src/main.f90`, time_step: state is published for Domain consumers and then
  catalogued, extracted, migrated and retained again after adaptation.
* `src/parallel_block_mpi.f90`, validate_candidate_block_scalar_wavelets:
  the refreshed Domain boundary is imported into blocks and block SOL ghosts
  are refreshed before the native scalar transform.
* The same file, prepare_block_native_wavelet_compression: the wavelet field
  is imported from Domain ownership into compact block storage. Compression
  runs there; activate_block_native_wavelet_compression writes it back to
  Domains before native inverse boundary preparation and reconstruction.
* The scalar tendency chain captures shared inputs and performs block-native
  direct flux, restriction, divergence and repeated ghost refreshes. This
  introduces setup and data traversal beyond the stencil expressions.

These are actual producer/consumer and ownership transitions, involving local
traversal, address work, copying and possibly MPI. They are not all synonymous
with message volume, nor does removing a call preserve required consumers.
The current timers do not establish what fraction can be eliminated safely.

The detailed profiles remain diagnostic: their labels split work differently
between implementations, and they have substantial unequal observer effects.
Even named MPI waits cannot distinguish network service, late peer arrival
from computation, packing delays or progress. Rank-local event epochs are not
synchronized absolute timestamps. Hardware instruction, cache and bandwidth
measurements were not collected. A causal explanation of the entire gap is
therefore not yet complete.

## Recommended continuation

Preserve Stage 195 and this review as a diagnostic milestone before changing
the production path. Continue the direct comparison; do not select another
small arithmetic kernel or a message-volume reduction from these totals alone.

The next diagnostic should partition the complete dynamics WT transaction
on both sides into numerical transform/reconstruction work, geometry/topology
preparation, allocation/initialization, representation/packing/install work,
and MPI posting/wait. Keep all subparts inside the conserved common envelope,
with clocks around batches and operation counts rather than per-cell timers.
Track the same field's producer, ownership, generation and next consumer to
identify repeated preparation and unnecessary round trips concretely.

Apply the same method to the dynamics intervals between WT calls and the
post-adaptation lifecycle. Compare complete RK-stage production/consumption
boundaries rather than subtracting the two current tendency tags. Keep
uninstrumented references and repeat coarse-versus-minimal calibration.
Only after that evidence separates dominant costs should a controlled change
to persistence, representation or ownership be used to measure recoverable
excess with unchanged numerical work and the existing correctness gates.

## Reproduction

The source tree is unchanged numerically. The review adds an archive auditor,
four focused tests and these evidence records. Four audit tests, three common
workflow tests and six previous complete-gap tests pass locally.

```sh
python3 test/parallel_block_profile/analyze_stage195_campaign.py \
  /private/tmp/stage195-review-54657 \
  --kit /private/tmp/wavetrisk-stage195-cluster-kit \
  --archive /Users/kevlahan/Desktop/stage195-review-54657.tar.gz
```

`parallel-blocks-cluster-195-54657-campaign.json` preserves the submitted
campaign. `parallel-blocks-cluster-195-54657-analysis.json` contains the audit,
paired phase ledgers, restart/remap subsets, transaction counts, detailed
regional means and resource summaries. No commit was made during this review.
