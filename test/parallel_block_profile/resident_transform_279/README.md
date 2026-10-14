# Stage 279 resident provisional transform campaign — kit v3

Prepared 2026-10-14 for local delivery. No upload or submission is implied.
Stage266 remains accepted. Previous is the exact job55504/v7 source; candidate
adds the connected resident provisional transform. Main is unchanged.

The repair moves all three provisional RK4 forward transforms, compression and
inverse operations to numerical owners and retains canonical results for the
next connected RHS. It bypasses the old provisional WT path, including its
vector staging, coefficient round trip, geometry-owned inverse and interior
installation/publication. It preserves the explicit raw-RK boundary and fixed
scaffold adapters, and the final transform/adaptation/remap/physics path. This
is a qualified slice toward whole-timestep independence, not completion of it.

Canonical inputs use prebound local reads, avoiding a general provider exchange.
Fixed all-column kernels preserve arithmetic order. Compact routes use one
borrowed communicator per epoch and reuse a single physical-owner directory
across phases. Their send/receive scratch is shared and sized for the largest
route, with all requests completed before reuse. Plans invalidate on topology/ownership or mask epochs. Runtime
uses prebound peer exchanges; it has no per-value owner lookup or dependency-
depth interpreter. Canonical Domain interiors are not published; the remaining compatibility and
block-boundary views are refreshed explicitly.

The exact source delta and local test hashes are in IDENTITY.json and
LOCAL_VALIDATION.json. ENTRY_GATE.json must pass before launch. The established
GCC13 physics dependency fix and identical 23-object Simple-physics build remain.
Never overlay an arbitrary working checkout onto this kit.

## Validation correction after job55506

The solver completed83-rank validation and passed science, but the evidence
checker incorrectly required every canonical owner to have scalar AND inner AND
outer records. Some owners legitimately have only outer-vector records. V3
requires nonzero total records, preserves the per-block census and adds exact
global class-by-class census agreement across all36 transactions and validation
rank counts. It retains all science/epoch/traffic/zero-Domain checks. Solver
source, flags, fixtures and the49-run schedule are byte-for-byte unchanged from
v2. The runner and reviewer now share one revision identity to prevent drift.

## Prior compiler portability correction

Job55505 stopped during the GCC13 checked build before any solver runs. GCC13
flagged allocatable descriptor self-reassignment in compile_fixed_inputs. V2
uses explicit allocation and stable ordered filling, including empty selections.
No warning or runtime check is disabled. Only that source file differs from v1;
the full failure record and exact-source GCC13 checked/optimized compiler checks
are included. Fresh local optimized/checked/oracle replays bind the correction.
Build failures now print their compiler-log tail in the Slurm output.

The prior v1 local timing comparison remains inconclusive: candidate point
estimates were +2.11% median/+4.53% mean versus v7, with overlapping ranges under
large host-load variation. Those timings are historical v1 evidence, not a new
v2 performance result. No speedup or cluster performance acceptance is claimed.
The49-run schedule, controls, fixture and performance gates are unchanged.

## Campaign and decisions

Use the established GCC13.2.0/MVAPICH4 environment and fixture
~/stage187-fixture-20260911-145307. One exclusive nk/bb13 allocation, all node memory, 256 CPUs, verified physical-
core/NUMA placement. Five serial builds and 49 serial solver runs:

- Nine validation runs, including checked16, oracle16/256, and audited
  production83/160/256. Oracles compare the old path and restore candidate SOL
  for real consumers; only production validations establish paired retirement.
- At83/160: two main and two266 controls bracketing the candidate comparisons.
  At256: two266 controls; main cannot run above160 Domains.
- Four bare runs each of v7 and candidate per count, in ABBABAAB order.
- One separate coarse-profile pair at each count; no fine profiling or adaptive
  additions. Profiles and audits never determine production speedup directly.

Primary: mean of independent run medians over nine ordinary steps, excluding
indices0,7,8 and retaining remaps4,9. Also report means/totals, remap/non-remap
subgroups, all run observations and spread. Production checkpoints are rare;
checkpoint-heavy validation costs do not replace ordinary-step performance.

Retain exact accepted-family science, main-relative tolerance, all12 states,
2 reloads,2 remaps and62 VTK members. Validate all36 provisional transformations
per rank and useful transform/RHS work on96 zero-Domain ranks at256. Record actual
per-block transform descriptors plus existing RHS work, traffic, named payload
storage, rank RSS/CPU/faults, placement and Slurm accounting. Descriptor counts
are not FLOPs or calibrated time weights; route calls are not global sync rounds.

Regions85–90 separate transform plan/input/forward/inverse/output/exchange costs.
Use self times on the same slowest whole-step rank, and compare complete-step
bare times against v7,266 and main. The previous binary has no regions85–90.
Do not interpret MPI wait as pure network time. Named buffer storage excludes
compiler peaks, kernel/geometry metadata, descriptors, the retained final-transform workspace and MPI;
all-rank RSS includes those costs, but summed individual peaks are not simultaneous.

Success requires defensible ordinary absolute-time improvement, preserved science
and ownership, acceptable memory and improved candidate160→256 scaling, assessed
separately. No automatic promotion or invented efficiency threshold. The broad
old transform budget contains necessary work and other callers, so it is not a
savings forecast. This slice cannot be promised to close the full140ms160-rank
main gap. Fixed640-block granularity remains unchanged; operation records inform
its separate repair. J6/J8 and shared-memory windows are outside this campaign.
J5 >=640 matched multi-node scaling remains an eventual required gate.

Budget approximately50–70 minutes,90-minute wall limit; queue and collection are
extra. Job55504 took42:33 with five builds/48 runs. New kernel compilation,
oracle checks and one additional validation justify extra allowance; this is a
planning estimate, not a measured cluster duration.

## Preserve failed job55506 before relaunch

On bbserv, run the original collection command once; it preserves the failed
build and satisfies the launcher's previous-job evidence gate:

```bash
bash "$HOME/rk279-resident-transform-kit-v2/collect.sh"
```

## Transfer from the Mac

```bash
scp /Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage279-validation/rk279-resident-transform-kit-v3.tar.gz /Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage279-validation/rk279-resident-transform-kit-v3.tar.gz.sha256 kevlahan@bbserv.math.mcmaster.ca:~/
ssh kevlahan@bbserv.math.mcmaster.ca
```

## On bbserv

```bash
set -e
cd "$HOME"
sha256sum -c rk279-resident-transform-kit-v3.tar.gz.sha256
test ! -e rk279-resident-transform-kit-v3
tar -xzf rk279-resident-transform-kit-v3.tar.gz
python3 ~/rk279-resident-transform-kit-v3/check_kit.py
bash ~/rk279-resident-transform-kit-v3/launch.sh
```

```bash
bash ~/rk279-resident-transform-kit-v3/status.sh
```

After the job ends, including if it fails:

```bash
bash ~/rk279-resident-transform-kit-v3/collect.sh
```

Collection verifies the archive and copies it to the established murier Desktop
destination. Inspect failures before another submission. The immutable kit and
all previous campaigns remain separate.
