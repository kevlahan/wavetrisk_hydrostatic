# Stage 196 cluster review — job 54659

The unchanged production implementations remain 2.264× apart: blocks 1.260 s
versus legacy 0.5565 s by median ordinary-step time over six bare repeats.
Stage 196 substantially narrows the location of that gap. Boundary/ownership
transfer routines account for 73.2% of the excess inside dynamics WT calls,
but only 44.1% outside those calls. Across the full step, their signed net
contribution is about 49% of the measured mean excess. This is evidence for
examining ownership transitions, packing and dependency scheduling together;
it is not evidence that reducing message bytes alone will recover that time.

## Evidence

Archive: `stage196-review-54659.tar.gz`, SHA256
`b03f38ab77f2421926d7780b7bd3d8ebc73d134107d549a473ec7c4cc1225537`.
All 17,173 inventoried files pass size/hash checks. Five build logs and source
identities match the frozen Stage 196 kit and Stage 195 baseline `5dfe0db1`.
The observer source is identical in both diagnostic implementations. Builds
use GNU Fortran 13.2.0, J5/RK4, contraction off, O2/march=native for optimized
and O0/runtime checks for checked, with 23 identical recorded physics objects.

All 58 runs and 146 recorded gates pass: 60 exact checkpoint comparisons,
58 workload comparisons and 28 fixed-tolerance legacy screens. The exact gates
include checked/optimized and observer/reference comparisons. Legacy remains
within the existing numerical screen, not bitwise equal: maximum reported
absolute differences are 3.723244e-12 velocity, 2.501111e-12 mass and
7.858035e-10 temperature (rounded upwards). The archive excludes checkpoint
payloads and executables: reports and identities were audited, not original
checkpoint arrays recomputed.

Every run completed 12 steps, restarting checkpoint 3, saving/reloading 4 and
remapping after restart. All 3,652 common workflow files and 3,652 operation
files were independently reparsed. Required scope coverage, flags, dense
records and per-rank/per-context conservation pass. Recomputed operation
summaries agree with the supplied results, allowing only floating-point
summation roundoff. Production timing cases have all validation oracles off.

All ranks ran on bb13 with the same 83 distinct core bindings, 0–41 and
128–168. No recorded swap-in/out changes occurred. Median summed rank-peak RSS
is 21.683 GiB blocks versus 11.016 GiB legacy; this is not concurrent RSS or
proof of a memory-bandwidth bottleneck.

## Observer effects

| Mode | Blocks, s/ordinary step | Legacy, s/ordinary step |
|---|---:|---:|
| Bare, printed median | 1.260 | 0.5565 |
| Minimal, common wall median | 1.260792 | 0.569607 |
| Boundary operations, common wall median | 1.268540 | 0.569582 |
| Full operations, common wall median | 1.264064 | 0.569648 |

Full operation observation versus minimal changes whole-step medians by
+0.260% blocks and +0.007% legacy. Within the matched WT context, observed
changes are +0.077% and +0.009%; outside WT, +0.357% and +0.296%. Boundary-only
and full-operation attribution also agree closely: WT boundary excess is
0.204776 versus 0.205488 s; outside-WT boundary excess 0.118232 versus 0.117894 s.
These observations support the broad attribution, unlike the substantially
unequal Stage 195 detailed-profiler effects. They are not formal upper bounds.

The minimal legacy executable itself is about 2.34% slower than bare by the
printed-step statistic. Disabled observer hooks/common clocks and executable
layout can still perturb execution. Thus the following is an instrumented
mean ledger, not an exact decomposition of the bare median gap; no uniform
overhead correction is applied.

## Additive full-operation accounting

Mean over the same 11 ordinary steps, 83 ranks and six repeats: blocks
1.207550 s, legacy 0.548173 s, excess 0.659377 s. Paired excess ranges from
0.654795 to 0.663286 s. Checkpoint steps are separate. The ordinary mean differs
from the bare median excess of 0.7035 s because the statistic and observation
mode differ.

| Disjoint part | Excess, s/step | Share of mean excess |
|---|---:|---:|
| WT boundary/ownership transfer routines | +0.205488 | 31.16% |
| WT other kernel, staging and residual work | +0.075154 | 11.40% |
| Outside-WT boundary/ownership transfer routines | +0.117894 | 17.88% |
| Outside-WT other kernel, storage and residual work | +0.149248 | 22.63% |
| Adaptation/reconstruction | +0.114188 | 17.32% |
| Other outer phases, signed net | −0.002595 | −0.39% |
| Total | +0.659377 | 100% |

The two WT rows are one matched context, excess 0.280642 s. Boundary/transfer
work therefore explains 73.22% of its excess. Outside WT totals 0.267142 s,
with 44.13% in boundary/transfer work. Together these contexts are 83.08% of
the full mean gap. These are signed differences, including legacy savings.

Within the boundary category, approximately 0.200590 s net is in named MPI
scopes and 0.122792 s in preparation/install/publication scopes. The latter
include mixed routines and MPI posting in mass begin, so this is not a strict
CPU-versus-network decomposition. MPI waits can reflect late producers,
packing delays, progress and synchronization. Neither quantity is an estimate
of safely removable time.

## What dominates inside WT

Selected self costs, seconds per ordinary step averaged over ranks/repeats:

| Operation | Blocks | Legacy in that group | Signed excess |
|---|---:|---:|---:|
| Inverse transfer planning/packing/install | 0.055320 | 0 | +0.055320 |
| Inverse MPI post/self-copy/wait | 0.091748 | 0 | +0.091748 |
| Block ghost preparation/install plus MPI | 0.034282 | 0 | +0.034282 |
| Domain writeback MPI | 0.026777 | 0 | +0.026777 |
| General Domain boundary CPU plus MPI | 0.015826 | 0.054536 | −0.038709 |

These are selected rows, not a complete ledger. Inverse wait alone is 0.088324 s;
packing is 0.026674 s and install 0.019598 s. Inverse completion scopes average
177.27 per rank/step across all four WT calls. That count includes calls with
no requests; it is not 177 network synchronizations on every rank. Plans are
already reused: one build and three reuses per rank/step. Repeated plan rebuilding
at every RK call is therefore not the main inverse cost.

The measured full/boundary imports and publication add concrete representation
work. Per rank/ordinary step in WT: eight full-import calls, four boundary-import
calls and twelve publication calls. Full imports zero about 1.354 million
buffer doubles; boundary imports another 1.631 million. These logical counters
show real work but are not hardware DRAM-traffic measurements.

The common Domain boundary path does less work in blocks, not more: its recorded
sent doubles in WT fall from 1.217 million to 0.359 million per rank/step, and
posted requests from 337.31 to 64.28. The extra specialized paths more than
consume that saving. This rules out blaming inefficient repeated use of the
old generic boundary interface alone. Blocks use reusable specialized paths;
multiple storage/ownership representations are the broader issue to test.

The other 0.075154 s of WT excess includes native wavelet/compression self work,
inverse numerical loops/traversal, scaffold staging and residual driver work.
Legacy arithmetic lies partly in driver/inverse residual scopes, so a zero in
a block-specific kernel group does not mean legacy performs no corresponding
arithmetic. The full WT envelope remains the matched comparison.

## What dominates outside WT

The scalar restriction/divergence/traversal group is 0.096169 s blocks versus
0.013895 s legacy, a 0.082273 s difference. Scalar setup/storage contributes
another 0.025816 s on blocks; the broader input-storage/geometry/capture group
is 0.026273 s. These are substantial computational/data-layout candidates.
However, block scope 25 contains direct flux, restriction, divergence,
temperature-boundary work and associated traversal; legacy scope 25 brackets
individual cpt_or_restr_scalar/flux calls. Their ratio is not a matched
arithmetic-kernel slowdown and does not prove duplicated calculations.

Other selected block costs are native mass MPI_Waitall 0.074860 s, scalar
transport MPI 0.045544 s, scalar transport/install 0.028884 s and block ghost
preparation/install 0.035710 s. General Domain boundary CPU/MPI saves
0.105756 s relative to legacy here. Again, adding only positive block costs
would overstate their net contribution to the gap.

Mass waits average 660.72 nonempty Waitall calls per rank/ordinary step, with
8,402.37 posted mass requests. The source schedules exchanges by position,
horizontal level and vertical layer. This is a concrete reason to measure
batching and producer arrival, rather than assume byte volume is the cause.
Scalar replay labels do not imply that a complete legacy tendency is also
being calculated: the production full validation path remains disabled.

## Counter limitation discovered in this review

Counters 4, 5, 12, 13, 14 and 15 are unavailable in this campaign. Their
call sites remain guarded by detail_enabled, which is intentionally false:
local geometry writes/field samples, generic boundary pack/route/install
candidates and generic boundary buffer zeroing therefore record zero despite
real work. See comm_mpi.f90:1037, 1066, 1094, 1510 and
parallel_block_mpi.f90:20980. Existing timing scopes are independently enabled
and conserved; boundary send/receive/request counts and new import/mass
counters are present. The frozen kit is preserved; these zeros must not be
used to calculate per-cell or per-route efficiency. The audit records this
limitation explicitly. The next measurement revision must enable those
call-site counters independently and check their nonzero coverage.

## Recommended continuation

Commit Stage 196 as a validated diagnostic milestone with this counter
limitation. No production performance improvement is claimed, and no new
optimization has been committed by this review.

The next comparison should test two concrete mechanisms, retaining untouched
Stage 196 numerical and legacy references:

1. Ownership and dependency scheduling: record request-bearing exchange batches
   by level/field/layer, useful local work between post/wait, packing/install
   counts and per-rank work distribution. Compare a change that eliminates an
   ownership round trip or batches a legal dependency phase. Measure full
   WT/outside-WT envelopes as well as transfers; smaller byte counts alone are
   not the success criterion. Avoid adding synchronization to measure arrival.
2. Matched scalar computation: bracket the same direct-flux, restriction,
   divergence and temperature tasks on both sides, separating traversal,
   lookup, initialization and numerical loops. Collect cells/edges/fields
   visited and actually evaluated, with the missing counters repaired. Only
   then attribute differences to redundant work, address lookup or arithmetic.

The 0.114188 s adaptation/reconstruction excess remains outside this stage's
inner-operation coverage. Keep it in the accounting and investigate its
rebuild/migration/publication work separately; it is not unexplained dynamics.
There is still no evidence for a complete duplicate legacy step. The evidence
now favors an expensive hybrid ownership/data-layout workflow plus localized
scalar work, but the removable fraction and hardware cause remain unproven.

## Reproduction and records

Run from the repository, after extracting the review archive:

```sh
python3 test/parallel_block_profile/analyze_stage196_campaign.py \
  /private/tmp/stage196-review-54659 \
  --kit /private/tmp/wavetrisk-stage196-cluster-kit \
  --archive /Users/kevlahan/Desktop/stage196-review-54659.tar.gz
```

The saved analysis JSON contains signed group ledgers, individual region
means/calls, work counters, calibration, memory and phase accounting. The
campaign JSON preserves the supplied results semantically, compactly encoded;
its byte digest differs from the original, whose digest is in the audit.
Continue collecting future reviews as one inventory-bearing tarball.
