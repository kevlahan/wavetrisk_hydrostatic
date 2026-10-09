# Stage275: input reuse passes through 256 ranks

Review date: 2026-10-09. Job **55485** completed all eight required cases.
Independent review verified **2,334 sealed files**, reproduced the minimal-input
budgets and joined the raw resident/Domain input records. The accepted numerical
baseline remains Stage266; Stage275 is a diagnostic overlay and has no qualified
performance comparison.

Archive: `rk275-review-55485.tar.gz`, SHA256
`eadae59d76427e6d090b1a726e160e5cf6195d9b16cdeae3e46502905e1ba912`.
Exact values and provenance are in `parallel-blocks-stage275-results.json`.

## Scientific and execution checks

All runs complete 12 steps, two remaps and 62 VTK outputs. Main-relative
scientific tolerance gates and exact Stage266-family checkpoint comparisons
pass, including the 83/160/256-rank cases. Checked execution passes. The 16- and
256-rank oracle cases each retain 48 exact inverse transactions and 84 oracle
installations, with storage and publication policies unchanged.

Each of the six instrumented runs verifies **672,000 source calls** and
**230,400 gradient calls**, spanning all 48 RHS transactions and 30 atmospheric
layers. Supplying only the 19 used components of each declared 3-by-9 flux/PV
stencil, while poisoning unused input locations, preserves every checked output
bit. Gradient shadow calls likewise preserve outputs using only active stencil
nodes. These checks execute on Domain owners and do not demonstrate numerical
execution on the additional block ranks.

Binary/checkpoint/VTK payloads are omitted from the archive. Their recorded
hashes and numerical reports were verified; the cluster solvers were not rerun
locally during review.

## Resident scalar state is correct in the sampled interiors

Each configuration supplies **424,580 matched mapped input records**, with zero
missing records and zero differences in scalar state, mean or their sums.
These are repeated observations across 16 phase/layer combinations, not 424,580
distinct mesh nodes. The phases are RHS calls 1, 2, 3, 9, 17, 33, 41 and 48,
at layers 1 and 30. All three 16-rank build/oracle variants pass.

| MPI ranks | Matched mapped records | Records with different Domain and block owner ranks |
|---:|---:|---:|
| 16 | 424,580 | 214,108 |
| 83 | 424,580 | 344,510 |
| 160 | 424,580 | 410,882 |
| 256 | 424,580 | 417,340 |

A separate reconstruction of the 256-rank raw files identifies **96,618** of
those matching records on the **96 ranks owning no Domains**. Their status was
derived from each rank's gradient-check counter increment at the same phase:
one increment per local Domain. Every sampled phase has 160 Domain owners and
96 zero-Domain ranks. No assumption about rank numbering was used to classify
them. The observed zero-Domain ranks happen to be 160–255.

This establishes sampled scalar correctness on block owners beyond the Domain
count. It does not establish balanced work or speedup: those ranks hold 22.756%
of the matched input observations, while making up 37.5% of ranks. Input counts
are not time-weighted numerical work and cannot quantify a load-balancing gain.

The scope remains **owner interiors**. Destination halos, unmapped boundary and
scaffold aliases, intervening calls/layers, and derived flux/PV/physics/Bernoulli/
Exner fields are not proved current by this comparison. Those checks must
accompany the executable producer prototype.

## Input trimming is insufficient for a velocity-only relocation

| MPI ranks | Reduction from full-argument to minimal Domain-fed input budget | Minimal input / current tendency export | Hypothetical upstream-owned input / export |
|---:|---:|---:|---:|
| 16 | 4.05–4.13% | 2.98–3.22 | 0.50–0.52 |
| 83 | 4.52–4.56% | 3.27–3.46 | 0.64–0.65 |
| 160 | 4.98–5.01% | 3.12–3.38 | 0.66–0.69 |
| 256 | 5.45–5.47% | 3.21–3.48 | 0.76–0.79 |

Each range covers the three rank-local structural sample groups. Ratios use the
same observations in numerator and denominator, per layer. Groups are not
synchronized mesh censuses or call-weighted trajectory totals. The hypothetical
upstream scenario assumes unimplemented field ownership and excludes upstream
production and communication costs. None of these percentages is a speedup.

At 256 ranks, flux and potential vorticity together account for **52.84–52.91%**
of the remaining minimal Domain-fed input budget; density and density-temperature
account for **15.41–15.45%**. Even an unrealistically free removal of the latter
would leave substantial Domain-fed traffic. The new scalar evidence supports
local reconstruction but cannot be credited as eliminating every destination's
halo transfer.

## Recommended next change

Proceed to a **bounded executable fused block primitive-producer prototype**,
targeting flux and PV while preserving the useful shared work in `step1`.
Stage270 independently measured about 194 ms of self time in the shared producer
region; Stage275 now identifies flux/PV as the largest proposed input-volume
pair. Together these justify prioritizing the producer/consumer boundary. They
do not identify a recoverable 194 ms or prove that flux/PV dominate runtime.

The prototype should compile compact patch/stencil addresses at valid topology,
mask and ownership epochs, form scalar sums directly from resident state/mean,
and validate scalar/vector halo and alias inputs at their actual consumption
phase. Use direct fused loops and retain the exact flux restriction, physics,
Bernoulli and Exner sequencing. Keep temporary duplicate calculations within
validation; do not accept them as a production migration.

The deliverable is executable block-owned primitive work with exact comparisons
at the producer/consumer boundary and an explicit cost account for preparation,
input refresh, local computation and publication. Integrate it with direct local
consumers before retiring the corresponding Domain work and tendency export.
A production candidate must pass uninstrumented full-step comparisons against
Stage266 and supported-rank main, including 160/256 scaling. Stop a repeatable
regression within that stage. No extra generic timing survey is needed to choose
this target; the unresolved halo/alias checks belong with its prototype.

Stage273 remains separate and unpromoted. Efficient execution at >=640 ranks for
J5's 160 Domains remains unqualified and will require multi-node testing.

## Runtime

Slurm reports **00:14:35**, from **2026-10-08 18:05:40** to **18:20:15** on bb13,
partition nk. Runner elapsed is 874.363 seconds. This includes four serial builds
and eight validation runs, but excludes queue wait and collection/transfer.
It is below the conservative 20–30-minute estimate. For a comparable diagnostic
schedule, use roughly 15–20 minutes; do not apply that estimate automatically to
a new numerical prototype or multi-node campaign. The exact record is saved in
`parallel-blocks-campaign-runtimes.json`.
