# Stage224: repeat the capability gates with corrected checkpoint metadata

Stage223 job 54725 passed the 40-rank baseline/candidate and all placement
probes, then failed at 41 ranks during restart. The source rank owned one
Domain and attempted to reread metadata at the consumed record cursor. The
Stage224 correction preserves the original metadata position.

The local checked MPI regression reproduced EOF with the original source in
both one-Domain-source layouts. Corrected layouts (source 0/empty rank 3 and
source 1/empty rank 0) completed 12 steps, two reloads and two remaps, matched
all restored metadata bytes and produced final checkpoints exactly equal to
the validated reference. This justifies repeating the bounded cluster gates;
it does not establish more-ranks-than-Domains capability yet.

The frozen reference is the original Stage223 reference with only
`src/checkpoint.f90` corrected. The separate prototype keeps the original
opt-in task-limit bypass and per-rank lifecycle evidence. Its four changed
files and corrected checkpoint source are checked against the local regression
by reproducing the test-only instrumentation and comparing source hashes.
The prototype deliberately retains `WAVETRISK_STAGE223_ALLOW`,
`WAVETRISK_STAGE223_LAYOUT` and `capability223-rank-*.txt`; the stage number of
the orchestration, source manifests, build identities and result archives is
224. No production guard change or automatic promotion is made.

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
