# Stage223 — cross the Domain/task limit

This is a bounded correctness and capability campaign. It freezes the locally
validated Stage222 cleanup as the baseline and changes the task-count startup
guard only in a separate prototype source tree. The Git `main` branch and the
production startup restriction are unchanged. Success does not automatically
promote the prototype or establish a performance gain.

## Question and evidence

Can ranks without geometric Domains own blocks, complete resident RK4 dynamics,
adapt, remap vertically, write a checkpoint, reload it and continue with the
same numerical solution and adaptive trajectory?

The prototype requires `WAVETRISK_STAGE223_ALLOW=1`. Native-layout runs retain
both existing distributors. One additional J4 stress case assigns Domains to
ranks 1–40, leaving rank zero empty; this changes only the Domain map, retaining
the real block distributor. Its checkpoint metadata source is therefore rank 1.

Each rank writes short lifecycle records at block installation, completed
multistage dynamics, completed steps and checkpoint load. The gate requires:

- Every rank owns nonzero blocks at every recorded installation/dynamics/step.
- Every rank completes every recorded RK4 dynamics transaction and timestep.
- Domain ownership counts sum to 40 or 160 at every lifecycle event.
- Rank clocks, output counters, lifecycle order and checkpoint metadata source
  agree; ranks without Domains remain full collective participants.
- Saved fields are finite; topology, thresholds and printed timestep/dt/mesh
  trajectories agree with the same fixture's baseline.
- Same-count comparisons are exact. Cross-count/layout comparisons use the
  existing fixed numerical policy, with actual exactness reported separately.
- The saved checkpoint is reloaded, followed by further evolution and vertical
  remapping. Merely reaching startup or writing a checkpoint is insufficient.

There are no broad profiles, cache variants or speed thresholds in this stage.
Lifecycle file writes and checked/oracle builds make these unsuitable for timing
claims. Saving is retained deliberately to exercise restart. Any later timing
campaign must exclude saving and first-restart steps and use run-level medians.

## Sequential run matrix

| Order | Fixture | Ranks | Build and layout |
|---:|---|---:|---|
| 1 | J4/J6, 40 Domains | 40 | Stage222 baseline, optimized |
| 2 | J4/J6 | 40 | Prototype, optimized, native layout |
| 3 | J4/J6 | 41 | Prototype, checked, native layout |
| 4 | J4/J6 | 41 | Prototype, optimized, native layout |
| 5 | J4/J6 | 41 | Prototype, optimized, existing oracles enabled |
| 6 | J4/J6 | 41 | Prototype, checked, rank zero has no Domains |
| 7 | J5/J7, 160 Domains | 160 | Stage222 baseline, optimized |
| 8 | J5/J7 | 160 | Prototype, optimized, native layout |
| 9 | J5/J7 | 167 | Prototype, checked, native layout |
| 10 | J5/J7 | 167 | Prototype, optimized, native layout |
| 11 | J5/J7 | 167 | Prototype, optimized, existing oracles enabled |

The first failure stops later runs and preserves partial evidence for review.
Each `srun` has a 15-minute limit to bound a collective stall; the batch budget
is four hours. There are no automatic repetitions or ABBA timing runs here.

## nk placement and builds

The user confirms nk has two sockets with 128 physical cores each: 256 total.
Reserve one exclusive node, 256 tasks, one CPU per task, no SMT. The build job
selects an nk node; the run job requests that same node for native compilation.

The launcher discovers socket/core/NUMA/L3 IDs, selects one hardware thread per
physical core and verifies that all 128 cores on each socket are available.
Ranks alternate sockets and are spread across the discovered LLC/NUMA groups.
The mappings use 20+20, 21+20, 80+80 and 84+83 physical cores for 40, 41, 160 and
167 ranks respectively. All mappings are probed before solver execution. Every
rank is pinned before MPI initialization; memory uses `numactl --localalloc`.

Six fresh solver builds cover J4/J5 baseline, optimized candidate and checked
candidate. Contraction is off, Stage202 is on, Stage203 is off. External physics
objects must be identical across all six. Runtime selectors retain the accepted
Stage215/222 choices. The kit contains the clean source snapshot, source hashes,
prototype patch, J4 restart fixture and J4/J5 coarse grids. J5 uses the existing
`stage187-fixture-20260911-145307` checkpoint/input on bbserv.

## Local preparation

The checked four-rank empty-root prototype completed the full 12-step J4/J6
trajectory, both checkpoint reloads and remapping, with an exact Stage215
checkpoint match. Rank zero owned blocks without Domains; rank 1 supplied
restart metadata. The larger-count runs remain cluster tests, not local claims.
See `parallel-blocks-stage223-local-results.json` for source/binary identities.

Review outcomes must distinguish capability from scaling: Domain-based physics
and other remaining producers can still constrain performance above 160 ranks.
If a gate fails, fix the specific empty-rank/ownership/collective assumption;
do not restore full Domain solution replicas on every block owner as a workaround.
