# Stage226: production support for ranks without Domains

Stage225 completed the 41/40 and 167/160 capability matrix exactly. This stage
integrates that support into production and tests larger task/Domain ratios.
Main remains unchanged; work is on codex/performance-investigation.

## Production change

The climate driver defers the old startup limit until input and dynamics
configuration are available. Other drivers retain the original early limit.
Main accepts ranks without Domains only for restarted climate Simple physics,
compressible non-mode-split RK4, physics enabled, no NCAR topography and no SSO.
This deliberately qualifies the tested path instead of silently enabling
cold starts or other dynamics. There is no environment-variable bypass.

Startup rejects ranks above N_CHDRN*N_GLO_DOMAIN. After every block catalogue
construction, including adaptation and restart, a collective-derived actual
block count is checked before ownership balancing. A run that coarsens below
its rank count stops with a diagnostic. Above the Domain count, each timestep
also requires the resident multistage block path before any Domain fallback.
No production lifecycle files, test ownership overrides or new Domain
replication are introduced. The Stage224 checkpoint-position and Stage225
extent corrections remain in place.

## Current decomposition limit

build_parallel_block_catalog examines root patch 1 of each geometric Domain.
If all four children exist, the four child subtrees become blocks; otherwise
the entire Domain is one block. It does not recursively split large subtrees.
Thus B = D + 3*S, where D is the Domain count and S is the number of split
Domains. D <= B <= 4*D. A block can contain many descendant patches; those are
not independently distributable blocks. The legacy arch variable n_glo_block
still counts geometric Domains; size(block_catalog) is the actual block count.

The implemented support requires tasks <= B so every rank can own useful work.
It does not claim support for ranks with zero blocks. J4's 40 Domains can
provide at most 160 blocks; J5's 160 Domains at most 640 blocks. A single nk
node caps non-oversubscribed tests at 256 tasks, so the 4:1 boundary is tested
with J4. Stage225 raw installation records independently show 160 total J4
blocks and 640 total J5 blocks at every recorded rebuild, so both fixtures
currently attain that maximum. Any failure stops later gates; do not weaken oracles or capacity checks.

## Cluster protocol

Six fresh builds: for each J4/J5, the frozen Stage225 reference, optimized
production integration, and checked integration with a separate test-only
lifecycle/ownership observer. The production and oracle runs do not use
WAVETRISK_STAGE223_ALLOW or write capability223-rank files. Checked observers
retain those lifecycle filenames and a layout selector solely for regression
evidence, without overriding either production rank guard.

| Fixture | Ranks | Runs |
| --- | --- | --- |
| J4, 40 Domains | 40 | Frozen reference and production |
| J4 | 41 | Checked observer, production, production oracles, checked empty-root |
| J4 | 80 | Checked observer and production |
| J4 | 120 | Checked observer and production |
| J4 | 160 | Checked observer, production and production oracles |
| J5, 160 Domains | 160 | Frozen reference and production |
| J5 | 167 | Checked observer, production and production oracles |
| J5 | 256 | Checked observer, production and production oracles |

All 21 runs require 12 steps, two checkpoint reloads, two remaps and the same
adaptive trajectory as the fixture anchor. Same-count comparisons are exact;
cross-count comparisons retain the existing fixed numerical policy and report
whether they are also exact. Checked observer records require useful blocks
and completed dynamics/steps on every rank, complete Domain coverage and the
expected empty ranks. Full oracles run at both established boundaries and new
maximum ratios. The J4 empty-root case retains nonzero checkpoint-source coverage.

Use the repaired whole-node srun policy, distinct physical cores across both
128-core sockets, the same build/run node, contraction off, identical physics
objects, Stage202 enabled and Stage203 disabled. No local O1 warning overrides
are applied to cluster builds. The per-run timeout is 15 minutes and batch
limit four hours. Archive collection records all results including a failure,
verifies hashes, and transfers from bbserv to murier; retries reuse the archive.

## Follow-up requested by the user

Find and review the earlier proposal for recursively splitting sufficiently
large refined block subtrees beyond four blocks per Domain. Keep that design
review separate from this production integration and existing-split validation.
A block root has at most four immediate children; the relevant threshold would
refer to descendant work/subtree size, not unbounded immediate child count.

## Partition choice

Keep Stage226 on nk: J4 at 160 tasks already exercises the current 4:1 ratio
and one-block-per-rank boundary; J5 at 256 tasks adds the larger fixture at the
whole-node limit. The user reports that bb can provide 640 tasks. After these
gates pass, a separate bb campaign can exercise J5 at 640 tasks/640 blocks,
with physical topology/placement verified for that partition and explicit
multi-node coverage if its allocation spans nodes. Do not reuse nk's fixed
two-socket placement map on bb.

## Completed local validation

The fresh production checked build (O1 locally, all runtime/initialization
checks and FP traps retained) passed the real 41-rank/40-Domain fixture with
all three oracles enabled. It completed 12 steps, both reloads and remaps,
matched the validated adaptive trajectory and produced an exactly matching
checkpoint. No Stage223 bypass, layout variable or lifecycle file was present.
Fourteen compiled MPI guard cases passed, including unsupported configurations,
the theoretical upper bound, the actual block limit and equality at capacity.
All four checked-observer modules compiled separately with warnings as errors.
All 55 production compilation commands retain the required check/FP flags.
The same five previously reviewed O1 warnings remain visible/nonfatal locally;
cluster build policy remains checked O0 with no such overrides.
