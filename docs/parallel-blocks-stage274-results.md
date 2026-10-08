# Stage274: ownership contract review

Job **55484** successfully recovered the diagnostic campaign. Independent review
verified all 1,790 sealed files and reproduced the ownership budgets from raw
plans. All eight required cases pass: main, Stage266, candidate, checked and
oracle at 16 ranks, candidate at 83/160, and candidate with oracle at 256.
Each completes 12 steps, two remaps and 62 VTK outputs. Exact Stage266-family
checkpoint checks, established legacy tolerance gates, storage and publication
checks pass. Each oracle case records 48 exact inverse transactions and 84 oracle
installations. This verifies diagnostic neutrality; it does not qualify performance.
Binary/checkpoint/VTK payloads are omitted from the review: their cluster-computed
hashes and numerical reports were checked, not independently rerun locally.

Archive: `rk274-review-55484.tar.gz`, SHA256
`8c96d310ebeed2eea9944fbd1fe50aa12dad9fb2d29b4b55da348443d17e5b9c`.
Compact evidence and exact budget values are in
`parallel-blocks-stage274-results.json`.

## Findings that guide the replacement

Across **3,228 sampled plans**, there are zero detected cross-owner restriction
operand transfers, zero source-to-gradient transfers, zero same-level cross-owner
source dependencies, and zero uninitialized source reads. All sampled numerical
targets map to blocks. Thus these samples support executing ordered source,
restriction and gradient work locally within the proposed block ownership.
They do not justify an execution-depth MPI interpreter. This is a sampled
structural result, not a proof for every future mesh, splitting scheme or epoch.
Unmapped boundary/scaffold **inputs** still exist; their canonical identity and
refresh phase need an explicit adapter even though target work is mapped.

The proposed full-argument input packing is substantially larger than the
completed velocity tendency exported today:

| MPI ranks | Input/export ratio with Domain producers retained | Input/export ratio with hypothetical upstream block producers |
|---:|---:|---:|
| 16 | 3.11–3.36 | 0.60–0.62 |
| 83 | 3.43–3.63 | 0.76–0.78 |
| 160 | 3.28–3.56 | 0.80–0.83 |
| 256 | 3.40–3.68 | 0.92–0.95 |

Ranges cover three rank-local ordinal groups. Each ratio divides summed proposed
source-plus-gradient input bytes by summed current tendency export bytes for
the **same observations**, per vertical layer. These are conservative protocol
budgets, not measured traffic, lower bounds, synchronized mesh censuses or
whole-step totals. The upstream scenario assumes field placement that has not
been implemented and excludes the costs of producing and communicating those
fields. It cannot support a speedup estimate.

At 256 ranks, the first group budgets 3,794,248 input bytes versus 1,031,808 export
bytes per layer. Hypothetical upstream ownership reduces that input budget to
980,448 bytes, only about 5% below the export budget; later groups are about 8%
below it. This provides little communication headroom without better input reuse,
although moving useful producer computation could still improve scaling.

The corrected reader is exercised on actual ownership changes: four Domains at
16 ranks and 104 at 83 ranks appear on multiple observing ranks. Each 16-rank
diagnostic has 492 plans; 83 has 792; 160 and 256 each have 480. Optimized,
checked and oracle 16-rank contracts agree. Ordinals and generations remain
rank-local; even groups with 160 observations are not claimed to be synchronized.

## Recommended next deliverable

Prepare a compact block-local source/restriction/gradient prototype on accepted
Stage266, with an **input residency and value-epoch gate before numerical
migration**. Establish which inputs can be reused or cheaply reconstructed on
block owners; for every remaining input, identify its producer, boundary refresh,
restriction phase, consumers and communication cost. In particular:

- Density and density-temperature are currently formed from `sol_mean + q` in
  `compute_native_velocity_gradient`. Prove equivalent block values and halo
  validity at that exact call before crediting local reconstruction.
- Source consumes native mass flux at each level's restriction phase. Gradient
  consumes final restricted Bernoulli and Exner. Neither can use an arbitrary
  earlier cached value. Preserve PV and physics callback dependencies explicitly.
- Preserve the direct arithmetic and ordered restrictions in compact local
  arrays. Keep upstream fusion and batch only the actual required boundary data.
  Do not move velocity alone using the conservative bulk input protocol above.

The prototype should report actual preparation, packing, communication and
execution cost alongside exact numerical comparisons. A production candidate
must retire the displaced Domain computation and tendency export in the same
qualified slice, then pass uninstrumented full-step gates against Stage266 and
supported-rank main, including 160/256 scaling. If required inputs force costly
bulk staging, revise the producer/consumer cut before promotion. This is a
bounded prototype recommendation, not approval to migrate every upstream
producer at once. No repeatable regression may be carried into a later stage.

Stage266 remains accepted. Stage273 remains separate and unpromoted. The eventual
target is efficient execution at >=640 ranks for 160 Domains; these samples do
not establish performance there or prove load imbalance is the sole bottleneck.

## Runtime and recovery

Recovery job 55484 ran **00:08:21**, from 17:20:28 to 17:28:49 on 2026-10-08
as reported by Slurm. It reused four verified builds and three successful solver
outputs from job 55483, revalidated those outputs, and executed the remaining five
cases. Its runner elapsed was 500.149 seconds.

Job 55483's runner recorded 349.662 seconds before the reader failure. Its Slurm
accounting query failed with a database connection refusal; exact batch elapsed
is unavailable. Combined measured runner time is 849.811 seconds, approximately
**14:10**, not an exact combined Slurm duration. Both records are preserved in
`parallel-blocks-campaign-runtimes.json`; queue and collection/transfer are excluded.
The 8:21 recovery is within the 6–12 minute recovery estimate and must not be used
as the duration of a fresh four-build/eight-run campaign.
