# Stage265 refinement: copy canonical results by patch

Status: local and cluster correctness passed. Job 55382 supports acceptance of
this refinement at the measured 83-rank configuration; no repeat recommended.
Legacy performance parity remains unachieved. Result details are in
`parallel-blocks-stage265-tiled-results.json`.

This remains within Stage265. Its first candidate (job 55381) passed science
but averaged 3.818 ms slower than Stage230, with mixed-sign paired differences.
It has not qualified for performance. No further Domain removal is proposed.

The new candidate changes only `src/parallel_block_inverse.f90` relative to
Stage265 commit `ea3d4b2d`. The original two-file return-removal delta relative
to accepted Stage230 remains. Both current block installations, all boundary
epochs, the independent old-return comparator, and the finite/extent preflight
are retained unchanged. Solver arithmetic is unchanged.

## Source hypothesis and implementation

The first publisher visits one workspace point, writes every variable and
vertical level in separate Domain arrays, then moves to the next point. This
touches many distant destination ranges before completing any adjacent field
range. It also selects restored-versus-current canonical data inside every
field copy. This is a concrete candidate for avoidable memory traffic and
descriptor/branch cost; the previous timings do not attribute the penalty to it.

The refinement compiles a small integer descriptor for each catalogued patch
with the inverse topology generation: Domain index, contiguous destination
start, and signed canonical source addresses for its points. Positive addresses
refer to the current workspace; negative addresses refer to `restore203`.
No payload pointers or numerical geometry are cached. Compilation checks full
patch coverage against the existing interior transfer inventory.

One reusable numerical tile holds a single physical patch:
`PATCH_SIZE**2 * (nscalar+EDGE) * zlevels` doubles. It is allocated only with a
new plan generation, shared by all local Domains, and fully overwritten for
each patch. No full Domain image or per-step allocation is added. Canonical
selection occurs once per point, followed by complete adjacent scalar-field
copies and adjacent edge writes. All original source finite-value and
destination-extent checks still precede authoritative writes.

## Gates and campaign

Local builds are serial, with at most four MPI ranks and one library thread.
Optimized, checked and independent oracle replays compare exactly against the
original Stage265 candidate, whose Stage230 equivalence is already verified.
Require all twelve states, reloads/remaps, checkpoint and 62 VTK members.
The existing oracle verifies all 48 returns, including scaffold and soil.
Local elapsed times do not establish performance.

The new kit bundles exact Stage230, the established main snapshot, the refined
candidate, all helpers and source seals. The previous kits are needed only
during local preparation; runtime on bbserv requires only the established J5
fixture. Seal checks enforce that only the inverse source differs from the
original Stage265 candidate. Preserve historical Stage264/265 kits and reviews.

Use exclusive nk/bb13, serial compilation of four solvers sharing the same 23
physics objects, and unchanged 16/83/167-rank physical placement. Seven
correctness cases remain: main, Stage230, candidate optimized/checked/oracle at
16, then candidate production/oracle at 167, including zero-Domain ranks.

The accepted baseline is still Stage230. Comparing only against the unqualified
Stage265 candidate could carry its existing penalty into later stages. The
initial unprofiled 83-rank sequence is therefore:

`main, Stage230, candidate, candidate, Stage230, main`

Thirteen solver runs total: seven correctness cases and six timing cases.
Every timing case retains scientific checks. Use ordinary steps excluding
1, 8 and 9, each run's precise median, two order-balanced candidate/Stage230
differences and the same-allocation main anchor. Any gain attributed to tiling
specifically would need a matched original-Stage265 comparison; this campaign
qualifies the complete refined return removal against the accepted baseline.

If the result is small/noisy/mixed, `repeat.sh` submits one more identical
six-run timing batch. It verifies and reuses the original binaries, source,
build logs, physics objects, toolchain and placement without repeating builds
or correctness runs. The original review stays immutable; initial metadata is
also preserved before the repeat. Each extra timing run passes science gates.
The launcher and runner reject a second repeat. No repeat is automatic.

A repeatable regression blocks further migration. Gain or supported neutrality
must be established against Stage230, with no accumulated main-relative
penalty. A mixed result remains inconclusive; no generic percentage allowance
or automatic promotion applies. Remaining regressions are resolved within this
same return-publication scope before considering the duplicate installation.

## Reviewed result: job 55382

Independent review verified 1,935 archive files, the exact delivered kit
identity and all thirteen solver runs. Both 16/167-rank oracle cases proved all
48 inverse returns, including twelve accepted exits. The checked case and all
scientific, storage, placement and zero-Domain-rank gates passed.

The six unprofiled timing runs averaged main 0.542783628 s, Stage230
0.767293997 s and candidate 0.753104779 s per ordinary step. Candidate/Stage230
pairs improved by 15.272575 and 13.105863 ms; the mean gain is 14.189219 ms
(1.849%). Every candidate median is below every Stage230 median. These results
support accepting the complete refined Stage265 change at 83 ranks without
the optional repeat. They do not establish speed at other counts or isolate
tiling from the original return-removal change.

The remaining main-relative gap is 210.321151 ms (+38.749%). This refinement
closes 6.320% of the matched Stage230/main gap; the inherited regression remains
the larger problem. Recommend using candidate commit `fab72486` as the next
immediate comparison baseline while retaining main in every timing batch.

The next focused proposal is to omit the first full block installation only
on provisional Stage203 production exits. The subsequent resident gather reads
canonical workspace values; the all-level boundary refresh then installs the
current block view. Audit every intervening consumer before implementation.
Retain the first installation for accepted exits, non-Stage203 execution and
the independent oracle, and retain all boundary epochs and the final install.
Measure this one omission separately with the same seven validation cases and
six timing runs; do not bundle another pathway removal or copy-order change.
