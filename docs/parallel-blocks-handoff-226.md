# Stage226 handoff: production rank-count integration

Production changes are in src/arch.f90, src/main.f90,
src/parallel_block_mpi.f90 and test/climate/climate.f90. The standalone patch
against the frozen Stage225 reference is
`test/parallel_block_profile/stage226-production.patch`. No main-branch change,
commit or merge is part of this integration; earlier cleanup/fixes remain.

The supported restarted climate Simple/RK4 path no longer needs a campaign
rank-limit bypass. The legacy drivers and unvalidated configurations retain
the Domain-count restriction. Startup enforces the present four-blocks-per-
Domain maximum; every catalogue rebuild enforces its actual global block
capacity before balancing; each surplus-rank timestep requires resident block
dynamics before any Domain fallback. Production contains no Stage223 hooks.
See `parallel-blocks-stage226.md` for exact supported configuration and limits.

Local validation and source identities are recorded in
`parallel-blocks-stage226-local-results.json`. The cluster package, once
prepared, is recorded in `parallel-blocks-stage226-kit.json`. Follow
`parallel-blocks-cluster-226.md` for transfer/build/run/retrieval. The higher-
ratio campaign has 21 gates, through J4 160 tasks/40 Domains and J5 256/160.
It uses a separate checked observer for per-rank evidence; plain production
and oracle builds do not include the observer. Checkpoint saving is intentional.

Requested later follow-up: find and review the earlier proposal for refining
large block subtrees beyond the present one-level/four-child split. Repository
history contains recursive subtree extraction at 20a34227 and deep-extraction
validation at 428c4cf3; these are leads, not yet identified as the proposal.
No deeper-splitting implementation is included in Stage226. Examine work/patch
counts below the root (each node has at most four immediate children), interface
closure, stable ownership keys, coarsening and migration before proposing a
threshold-driven split policy. Preserve current exact oracle/restart gates.

The user also reports 640-task availability on partition bb. The recommendation
is to complete nk first (J4 already reaches 4:1), then validate bb placement and
run J5 up to 640 tasks as a separate capacity/multi-node extension. Stage226
retains nk and does not assume bb's CPU or node topology.

Local completion: the actual 41-rank production checked full-oracle run passed
all twelve steps, both reloads/remaps and exact checkpoint/trajectory comparison
without prototype hooks. Fourteen compiled MPI policy/capacity cases also pass,
and the separate checked observer compiles. Source/build evidence is persisted.

Packaging completion: the self-contained Stage226 kit is prepared, not submitted.
Re-extraction reproduced baseline, production and checked-observer inventories.
All 19 packaged helper tests and the complete isolated bootstrap/launch/collect/
repeat-transfer workflow pass. Python and shell syntax checks pass. The verified
archive is supplied on the Desktop; transfer/build/run/retrieve commands are in
parallel-blocks-cluster-226.md. The workflow simulation submits no remote jobs.

## Cluster review update: job 54731

All 21 runs and 33 exact numerical comparisons pass; 6,386 archive hashes and
raw lifecycle/placement/trajectory evidence rechecked. Production capability
now passes J4 160/40 and J5 256/160. Corrected ordinary-step medians show poor
relative scaling: J4 40→160 gives 1.275x speedup/31.9% efficiency; J5 160→256
is 9.5% slower/57.1% efficiency. Each count has only one independent bare run;
exclude startup, save/reload step and next full step (1,8,9). See the Stage226
review, JSON and CSV for raw selections and limitations.

New hard constraint: NEVER use bb04 (user reports unreliable). bb01–bb12 are
reported as 2 sockets, 64 Slurm CPUs/node; exclude bb04 leaves 11 nodes/704
configured slots. A 640-task allocation needs ten nodes, not nk's one-node
layout. Future bb build/probe/run scripts must use --exclude=bb04 and validate
hostnames, physical topology, CPU compatibility and per-host local-rank mapping.
Machine-readable requirements: parallel-blocks-bb-campaign-requirements.json.
No bb campaign was prepared/submitted in this review. Do not infer performance
benefit from reaching the block-count ceiling; establish a bb-local reference
and keep capacity validation distinct from performance confirmation.

## Agreed optimization direction after review

Pause the bb/640-task campaign and deeper splitting while closing the remaining
performance gap with legacy. Prioritize complete redundant Domain traversals,
reconstruction, communication and intermediate storage, with an explicit
producer/consumer replacement and bounded feasibility evidence before another
performance campaign. Necessary physics must be preserved when moving its
execution to blocks. Keep the passed above-Domain-count cases as regressions.
The current J5 comparison at tasks=Domains means 160 MPI tasks but 640 blocks;
a literal one-block-per-Domain decomposition is a separate diagnostic change.
Use the agreed ordinary-step screening and A–B–B–A confirmation policy and
complexity-dependent acceptance threshold. Main must remain untouched. bb04
remains excluded whenever a future bb campaign is resumed.
