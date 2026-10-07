# Stage265 proposal: direct inverse-result publication

Status: source change and campaign proposed, not implemented or packaged.
Baseline: exact Stage230 commit `4ea640128a71a857f9b8c6b9411dc34fc0aeb641`.
Legacy anchor: main commit `6a245db7f1632b0205c2bebe95e5f0cd535c2582`,
with the same established checkpoint I/O hardening as Stage264.

Stage264 job 55374 independently verified the baseline. At 83 ranks on bb13,
Stage230 averaged 0.762952294 s/ordinary step versus main's 0.542412091 s:
a 220.540203 ms gap, or 40.659%. This proposal targets one redundant transport;
it does not assume that this transport explains the whole gap.

## Source change

Replace only the complete SOL block-to-Domain return in
`activate_block_native_inverse_transform`, after the first
`install_inverse203`. The native inverse result already lives on the geometry
rank that owns the required consumer. Publish from that local workspace rather
than exporting to final block owners and fetching it back.

Add a narrow canonical-interior publication routine in
`parallel_block_inverse.f90`. Reuse the current integer node/address maps.
For an interior covered by `restore203_index`, use `restore203`; otherwise use
the canonical workspace value. Copying `value` indiscriminately is incorrect
because boundary dependency updates can overwrite interior aliases.
`publish_native_inverse` currently skips interiors and cannot simply replace
the return call.

Retain both existing `install_inverse203` calls, the provisional all-level
boundary refresh, and accepted-step physics consumer alias restoration.
Keep their distinct value epochs. Keep all still-required collectives on ranks
without geometric Domains: these ranks can own final blocks.

Audit every destination of `write_domain_patch_prognostic`, including untouched
scaffold nodes and the global `sol` fields outside atmospheric levels. The old
writer reaches soil/non-atmospheric state even when passed an atmospheric
`domain_scaling` slice. Prove that omitted writes are exact no-ops at every
transaction or supply an explicit bounded replacement. This proof is required
before removing the old call.

Validate map generation, coverage and destination extents before publication.
Avoid a new per-value owner lookup, a new full production numerical image, or
per-step transport-plan allocation. Update the production writeback counters
and their assertions to describe the actual calls. Add diagnostic counts of
replaced returns and direct publications; enable their logging only for
validation. Verify the expected four inverse exits per ordinary RK4 step
rather than assuming the count also applies to restart/remap steps.

In the independent validation path, reconstruct the old return result from
final blocks and compare every affected destination bitwise at the transaction
boundary. Include canonical interiors, scaffold, boundary aliases and
non-atmospheric/soil state. Test-only scratch is permitted here; none is
retained in production. Preserve existing dynamics/adaptation/remap oracles.

No velocity, diffusion, primitive, adaptation algorithm, wavelet formula,
logging-policy or MPI-placement changes belong in this candidate. Coalescing
the duplicate installation is a subsequent, separately measured change.

## Local correctness gate

Use serial builds, four MPI ranks maximum and one thread per numerical library.
Run the full established restarted J4 fixture with exact Stage230 as reference:
candidate optimized, candidate checked, and independent oracles. Require exact
checkpoint and same-count VTK equivalence to Stage230, all twelve states,
two reloads and two remaps. Exercise the transaction comparator and publication
counts across provisional and accepted inverses. Local speed is not evidence.

## Proposed bbserv campaign

Prepare a self-contained Stage265 kit using the working Stage264 submission,
placement, status and collection structure. Bundle all source snapshots and
helpers, seal their hashes, and require only the established J5 fixture.
Do not depend on previous extracted kits, job state or mutable build trees.
Preserve fresh result directories, duplicate-submit protection and collection
of early failures. Test the transferred archive with the independent reviewer,
including JSON census serialization, before delivery.

Use one exclusive `nk`/bb13 allocation; no bb04 exclusion. Build optimized main,
optimized Stage230, optimized candidate and checked candidate with `make -j1`,
sharing the same 23 physics objects. Use Stage264's physical-core placement,
runtime selectors, inputs and arithmetic controls.

First run seven correctness cases:

| Ranks | Version | Purpose |
|---|---|---|
| 16 | main | Legacy numerical anchor |
| 16 | Stage230 | Immediate exact reference |
| 16 | candidate optimized | Production equivalence and removed-return evidence |
| 16 | candidate checked | Bounds/runtime checks and equivalence |
| 16 | candidate with oracles | Independent transaction and scientific checks |
| 167 | candidate optimized | Exact cross-count result and zero-Domain ranks |
| 167 | candidate with oracles | Zero-Domain collective and transaction checks |

Stop on any correctness failure. Require Stage264's fixed checkpoint hashes,
unchanged main-relative numerical tolerance, exact block comparisons, 62 VTK
members and same-count hashes, twelve states, two reloads/remaps, and absent
production Domain trend storage. Candidate oracle storage follows the existing
oracle contract. Validate physical placement and binary/source identity.

Then run two independent unprofiled timing batches at 83 ranks:

`main, Stage230, candidate, candidate, Stage230, main`

This is twelve timing runs, four per version, and nineteen solver runs total.
Every timing run also passes the scientific gates. Exclude steps 1, 8 and 9 as
in Stage264; use the precise `test_log`, compute each run's ordinary-step median,
then the mean of four independent medians. Do not time checked/oracle modes.
Report four order-balanced candidate-minus-Stage230 gaps, candidate-minus-main,
the main-relative ratio, all raw run medians and measured variability.

## Decision and next stage

A repeatable candidate slowdown blocks advancement. A clear improvement can
qualify this removal only after all correctness gates pass. Near-noise or
mixed-sign evidence receives at most one additional predeclared batch of the
same six-run sequence; call unresolved evidence inconclusive. Neutrality needs
evidence within the observed variability, not a generic percentage allowance.
There is no automatic promotion or claim of legacy parity.

If slower, diagnose the replacement cost within Stage265: direct publication,
preflight, remaining block install and boundary work. Any diagnostic profiling
is a separate matched follow-up; it is not included in acceptance timings.
Do not introduce another migration while this result remains unresolved.

After acceptance, consider removing the first duplicate block installation
only if its intervening consumers can be proven not to require it. Repeat the
candidate-versus-immediate-baseline gate and retain main as an anchor. Continue
until the full production path removes non-essential Domain dependencies and
matches or beats legacy main; a neutral pathway removal alone does not meet
that performance goal.
