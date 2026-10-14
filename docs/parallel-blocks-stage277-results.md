# Stage277 review: correctness passes; no qualified performance benefit

Reviewed 2026-10-09. Job55497 completed all sixteen runs in **18:45**.
The independent local reviewer verified all **4,323 sealed files**, science,
placement, assignment records and uninstrumented timing summaries. The archived
identity, kit manifest and source manifest match the delivered kit exactly.

Archive: `rk277-review-55497.tar.gz`, SHA256
`aa268d0d5e239ea4bd945ade6397e7921bcc323d99d4824cb96dcd9a5cb1fbb8`.
Machine-readable review: `parallel-blocks-stage277-results.json`.

## Decision

**Do not promote Stage277. Keep accepted Stage266 `ff740744`.** The arithmetic
and assignment evidence are useful for the complete ownership migration, but
this geometry-owner cache has no demonstrated full-step benefit and changes no
execution ownership. Preserve its direct kernels and tests as reusable evidence;
do not automatically inherit the cache, temporary storage or installation scheme.

One timing run per configuration does not prove a repeatable small regression.
Nevertheless, an additional campaign is unnecessary for the decision to leave
this candidate unpromoted: the observed changes are small and mixed, while the
strategic ownership benefit is absent. Do not spend further stages tuning this
standalone optimization before completing the connected ownership change.

## Scientific and assignment evidence

All eight validations and eight timing runs pass the required science checks,
including twelve steps, two reloads, two remaps, 62 VTK members, block-family
exact comparisons and supported main-relative comparisons. Checked execution
and the independent 16/256-rank oracles pass. Binary/checkpoint/VTK payloads are
omitted from the review; hashes and comparison reports were verified, not the
absent payloads or a new local solver replay.

Each of the six instrumented configurations checks **61,115,520 node observations
and 1,161,194,880 assignments**, covering all 48 RHS transactions and all 30
layers with zero differences. Cache rebuild counts are twelve on every rank.
At 256 ranks, all 96 zero-Domain ranks report **zero primitive nodes**. This
confirms the unchanged Domain-owner execution limit; the block-owner prototype
success from Stage276 must not be attributed to Stage277 production.

Timing runs have assignment checks disabled, the expected production selector,
and no added profiling or assignment output. Four serial builds share the same
23 physics objects. The separate Stage273 partition candidate is excluded.

## Same-allocation bare controls

| Ranks | Main ms/step | Stage266 ms/step | Stage277 ms/step | Stage277 minus Stage266 |
|---:|---:|---:|---:|---:|
| 83 | 544.225 | 760.845 | 757.689 | -3.156 ms / -0.415% |
| 160 | 356.235 | 551.163 | 557.192 | +6.028 ms / +1.094% |
| 256 | unsupported | 590.552 | 593.222 | +2.670 ms / +0.452% |

Stage277 remains 213.464 ms (39.22%) behind main at 83 and 200.957 ms (56.41%)
behind at 160. These are ordinary-step run medians under the established
exclusions, one independent run per configuration. Steps within a run are not
independent repetitions. No uncertainty interval or speedup qualification is
claimed for the small candidate/baseline differences.

Stage266 160-to-256 speedup is 0.9333; Stage277 is 0.9393. Both slow down with
more ranks. Relative efficiencies are approximately 58.33% and 58.70% against
their own 160-rank references, not measured serial efficiencies. There is no
evidence that the Domain-owner cache resolves the scaling problem.

## Runtime and next objective

Slurm reports 2026-10-09 10:39:48–10:58:33, **00:18:45 (1,125 seconds)**.
Runner time is 1,123.923 seconds. Four serial builds, eight validations and eight
single timing runs; queue and collection/transfer are excluded. This was slightly
below the 20–30-minute planning estimate. A comparable unchanged schedule can be
budgeted at 18–25 minutes; a substantial ownership replacement needs its own
estimate. Exact accounting is retained in the campaign runtime ledger.

The next objective is the connected block-owned production RHS chain, including
the physics/pressure/restriction consumers that would otherwise force a bulk
return to Domain owners. Its acceptance must demonstrate retired production work
and transfers as well as science and full-step performance. See the revised
[objectives](parallel-blocks-objectives.md) and maintained
[scorecard](parallel-blocks-scorecard.md).
