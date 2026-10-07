# Block refactoring strategy

User guidance reaffirmed on 2026-10-07:

- The primary goal is efficient execution on **four or more MPI ranks per geometric Domain** (`P >= 4*D`). The block data structure is intended to remove the geometric-Domain limit on useful parallelism.
- Supporting ranks without Domains proves correctness only. Success requires measured reductions in total step time and useful parallel efficiency beyond the Domain count. Matching legacy main below that count is a necessary regression check, not the end goal.
- Prioritize measured bottlenecks that restrict useful work to Domain owners, duplicate Domain/block work, require representation transfers or global synchronization, or prevent adequate block load balancing. Remove non-essential Domain pathways progressively while preserving numerical correctness and required oracle checks.
- Qualify each substantial change against the accepted candidate and legacy main where supported, and include performance beyond the Domain count. Resolve regressions within the stage that introduces them; do not stack changes and defer performance recovery.
- J5 has 160 geometric Domains: 256 ranks is only an intermediate scaling check; the fourfold target is at least 640 ranks and needs multiple nodes on the available cluster. Current four-child-subtree splitting provides at most four blocks per Domain, so finer subdivision/load balance must be considered for efficient execution at or beyond that target.

See `docs/parallel-blocks-measured-restart.md` for campaign history and the current strategy. No numerical efficiency threshold has yet been agreed; report measured speedup and efficiency rather than inventing an acceptance percentage.

Operational constraints: the local machine has four usable cores; use serial builds and at most four MPI ranks locally, without concurrent build/replay work. Partition `nk` contains only bb13 and needs no bb04 exclusion; exclude bb04 on `bb` or `all`. Provide user shell commands in fenced code blocks without `#` comment lines.
