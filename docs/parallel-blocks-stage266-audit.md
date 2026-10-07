# Stage266: defer the provisional first block installation

Status: audit implemented; exact optimized, checked and independent oracle
four-rank J4 replays passed. The existing canonical/boundary epoch regression
test passed at both `-O0` and `-O2`. Cluster job 55398 supports accepting the
change at 83 ranks with a 10.033 ms gain against the immediate baseline; no
repeat is recommended. Legacy parity remains unachieved.

Baseline: accepted tiled Stage265 commit `fab72486`. Job 55382 verified thirteen
runs and improved the 83-rank ordinary-step median by 14.189 ms against Stage230.
Its main-relative gap remained 210.321 ms. Preserve this immediate baseline and
main in the next campaign; do not combine other removals or layout changes.

## Consumer audit

The audited region is `activate_block_native_inverse_transform`, from its first
`install_inverse203` through the existing provisional boundary refresh and
second installation. An installation saves alias-overwritten workspace values,
substitutes canonical `restore203` values, exports scalar/vector interiors to
final block owners, restores the saved workspace, then exports scalar/vector
boundary views. It has four transfer calls and leaves the workspace unchanged.
Its side effect is the final-owner block image, not the canonical producer.

| Intervening operation | Reads / writes | Need the first installation? |
|---|---|---|
| Profile enter/leave | Profiling counters only | No |
| Early `publish_native_inverse(domain_scaling)` | Copies workspace boundary values to Domain fields | Skipped already for provisional production; retained for final/oracle |
| `publish_native_inverse(domain_wavelet, WAV_COEFF)` | Copies the wavelet workspace to Domain boundary fields | No; it reads the separate wavelet family |
| `publish_canonical_inverse203` | Checks and copies canonical `restore203`/workspace interiors plus existing scaffold to Domain scaling | No; integer patch descriptors and the tile do not read final-owner blocks |
| `assert_canonical_inverse_return` | Reconstructs the old return from final-owner blocks and compares Domain/scaffold/soil bits | Yes; oracle first installation must remain |
| Provisional `native_inverse_gather` for scalars/vectors | Resets scaffold and copies canonical restore lanes into resident workspace; production returns before transfer | No under Stage203 production; oracle gather comparison reads blocks |
| All-level `native_inverse_boundary` | Ordered local aliases and remote snapshots from workspace | No; preserve scalar/vector order and exact epochs |
| Oracle Domain boundary refresh/comparison | Independent Domain update and phase comparison | Retain unchanged, with both installations |
| Production boundary publication | Workspace to Domain boundary fields | No |
| `native_inverse_scatter(...,.false.)` | Under resident production, bookkeeping then return; oracle scatters existing boundary plan | No production dependency |
| Final `install_inverse203` | Exports canonical interiors and refreshed boundaries to final block owners | Required before the caller consumes the result |
| Final oracle block/Domain assertion | Reads installed block image | Retain after the final installation |

No intervening production path reads `saved203`: it is scratch private to the
installation. No new payload, transport allocation or collectives are needed.
The canonical tile, complete finite/extent preflight, scalar/vector kernels,
publication and writeback counters remain unchanged. The retained second
installation substitutes the same canonical interiors and the current boundary
workspace as before. Accepted exits have no later installation, so keep theirs.
Fallback/non-Stage203 paths retain the original sequence.

## Narrow implementation and actual-call check

Defer the first call only when all three conditions hold: provisional exit,
Stage203 inverse enabled, dynamics oracle disabled. Retain every boundary refresh
and the later installation unconditionally within the provisional branch.
Do not guard any retained collective by geometric Domain count: ranks with no
Domains can still own final blocks and must participate.

Add a private monotonic counter inside `install_inverse203`, incremented only
when the resident installation executes. A read-only accessor lets the caller
assert actual call deltas on every Stage203 exit and every rank. Require one
for all production exits and accepted oracle exits; require two for provisional
oracle exits. This would detect an extra installation introduced by a helper.
The counter does not alter `work203`, transport, numerical data or allocation.
One startup marker identifies the new production contract. Per-exit actual call
counts are logged only with the existing independent oracle, which also retains
all 48 exact return proofs. With the established twelve-step RK4 fixture the
production sequence has 48 full installations rather than 84, while oracle
execution retains 84. The actual per-exit assertions do not assume that count
for arbitrary runs. Finite/extent and scientific checks remain in timed runs.

The source delta consists of the two inverse files: the guarded call and audit
assertion in `parallel_block_mpi.f90`, and the counter/accessor in
`parallel_block_inverse.f90`. No numerical formula or layout changes belong in
this test. The source seal binds exact accepted baseline, main and candidate.

## Separate correctness and timing campaign

Use serial local builds, four MPI ranks maximum and one library thread. Require
optimized and checked production replays to match the accepted tiled reference
exactly, exercising the deferred branch. The independent oracle replay retains
both installations, checks each actual call count, and proves all 48 return
transactions including twelve accepted exits. Require all states, restart/remap,
checkpoint and 62 VTK outputs to match. Local timing is not performance evidence.

On one exclusive nk/bb13 allocation, build main, accepted Stage265t, candidate
optimized and candidate checked serially with the same 23 physics objects.
Retain seven validation runs: main/baseline/candidate optimized/checked/oracle
at 16 ranks, then candidate production/oracle at 167 ranks to exercise ranks
with zero geometric Domains. Retain fixed checkpoint hashes, numerical policy,
exact block comparisons, VTK members, placement, restart/remap and storage checks.

Run six unprofiled 83-rank timing cases:
`main, Stage265t, candidate, candidate, Stage265t, main`.
Compare precise ordinary-step medians excluding steps 1, 8 and 9, the two
order-balanced candidate/baseline differences, and the same-allocation main gap.
Use one optional extra batch only for small/noisy/mixed evidence, with verified
binary/build/toolchain/placement reuse and no repeated correctness runs. No
automatic repeat or promotion. Any repeatable regression blocks further removal.
The delivered kit is self-contained except for the established J5 fixture and
uses the previously validated launch/status/collection structure.

## Reviewed result: job 55398

Independent review verified 1,935 archive files and all thirteen solver runs.
Transferred source, kit and local-validation manifests match the delivered kit.
The 16-rank checked case passed. Both 16/167-rank independent oracle cases
proved all 48 inverse returns, including twelve accepted exits, and retained
the original 84-installation sequence. Boundary/scientific, checkpoint, VTK,
reload/remap, storage, physical placement and zero-Domain-rank gates passed.

The six timing runs averaged main 0.541753610 s, accepted Stage265t
0.758481670 s and Stage266 0.748448779 s per ordinary step. The two paired
candidate/baseline gains were 11.031772 and 9.034010 ms; the mean gain is
10.032891 ms (1.323%). Both candidate medians are below both baseline medians.
Recommend accepting this one omission at the measured 83-rank configuration
without the optional repeat and using candidate `ff740744` as the next immediate
baseline. Full evidence is in `parallel-blocks-stage266-results.json`.

The remaining main-relative gap is 206.695169 ms (+38.153%). This change closes
4.629% of the same-allocation baseline/main gap; the inherited regression is
still substantial. Do not infer speed at 16/167 ranks from correctness runs or
subtract results from different allocations to claim a cumulative gain.

Next audit the remaining provisional canonical workspace-to-Domain interior
publication and every subsequent reader. The native provisional refresh avoids
re-importing resident interiors, which makes this a plausible next audit target,
but it does not prove every other consumer is independent of Domain scaling.
Retain required boundary views and accepted-step physics/output adapters. Only
propose a separate removal after the complete reader/epoch audit, with the same
science gates, immediate-baseline/main comparison and regression stop.
