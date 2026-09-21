# Stage225: complete the above-Domain-count capability gates

Stage224 job 54727 passed J4 40-rank baseline/candidate and 41-rank checked/
optimized production runs, including both reloads and remapping. The 41-rank
oracle run then aborted on an extent check: a zero-Domain rank has zero logical
transform values but a one-element dummy allocation. Exit 137 came from Slurm
terminating peers after MPI_Abort, not an observed out-of-memory condition.

The correction compares allocated capacity with max(1, logical count), matching
the existing allocation contract. Both checked local four-rank regressions now
complete 12 steps, two reloads and two remaps with all three existing oracles
active. Source rank 0 or 1 owns exactly one Domain; receiver rank 3 or 0 owns
none. Restored metadata and final checkpoints match exactly. All ranks perform
block work. Ten predicate cases also pass. This is evidence to resume cluster
correctness checks, not evidence of a speedup.

The frozen reference is the Stage224 reference with only
src/parallel_block_mpi.f90 changed by the extent correction. The candidate uses
the same opt-in task-limit bypass and per-rank lifecycle recording as before.
Packaging reproduces the local regression instrumentation and requires all 136
Fortran source hashes to match the completed local tests. It also re-extracts
the source archive and reproduces candidate identities. Full source is included;
there are no Git bundle prerequisites.

The candidate deliberately retains WAVETRISK_STAGE223_ALLOW,
WAVETRISK_STAGE223_LAYOUT and capability223-rank-*.txt for compatibility.
Orchestration, source manifests, build records, state and archives use Stage225.
The production task-count guard remains unchanged; no branch is modified.

Repeat the quick J4 anchors before the previously blocked oracle gate, then
continue through empty-root and J5 gates. Do not reuse Stage224 binaries or skip
oracle checks. Local validation used checked O1 builds to bound runtime; cluster
builds retain the established DEBUG=check (O0) and optimized policies with no
local warning overrides. All six cluster builds are fresh.

## Sequential correctness matrix

| Order | Fixture / Domains | Ranks | Build / layout |
| --- | --- | --- | --- |
| 1 | J4 / 40 | 40 | Corrected reference, optimized |
| 2 | J4 / 40 | 40 | Candidate, optimized, native |
| 3 | J4 / 40 | 41 | Candidate, checked, native |
| 4 | J4 / 40 | 41 | Candidate, optimized, native |
| 5 | J4 / 40 | 41 | Candidate, optimized, existing numerical oracles |
| 6 | J4 / 40 | 41 | Candidate, checked, empty rank zero; source rank 1 |
| 7 | J5 / 160 | 160 | Corrected reference, optimized |
| 8 | J5 / 160 | 160 | Candidate, optimized, native |
| 9 | J5 / 160 | 167 | Candidate, checked, native |
| 10 | J5 / 160 | 167 | Candidate, optimized, native |
| 11 | J5 / 160 | 167 | Candidate, optimized, existing numerical oracles |

Stop on the first failure. Require useful block ownership and completed block
dynamics on every rank, agreement in global clocks/lifecycle, complete Domain
coverage, the expected zero-Domain participants, both checkpoint reloads,
remapping and matching adaptive trajectories. Same-count checkpoint comparisons
are exact; cross-count/layout comparisons use the unchanged established
numerical policy and report whether they are also exact.

Build all six binaries fresh on one nk node, with contraction disabled,
Stage202 enabled and Stage203 disabled. Physics objects must match across
builds. The run job selects that same node. nk has two 128-core sockets; reserve
one exclusive 256-core node. The verified mappings spread 40, 41, 160 and 167
ranks across both sockets, pinning one hardware thread per physical core.
Keep the repaired whole-node step policy: no step `--cpus-per-task` or `--exact`;
record/clear inherited step restrictions and retain probe stdout/stderr.

Per-step timeout remains 15 minutes; the batch limit is four hours. This is
restart/capability validation, so checkpoint saving is deliberately enabled.
There are no performance acceptance thresholds, ABBA repetitions or speedup
claims. Later performance work must exclude saving and first-restart steps.
