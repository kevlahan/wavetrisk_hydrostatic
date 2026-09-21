# Stage222 handoff

Read `parallel-blocks-stage222-cleanup.md` and the corresponding results JSON.
The local correctness gate has passed: all four trajectories match Stage215
checkpoints exactly. The checked no-Domain rank completes adaptation, both
reloads and remaps with 36–38 owned blocks. The next campaign is 41 ranks on
40 Domains, then 167 ranks on 160 Domains. Stage223 is now packaged; see
`parallel-blocks-handoff-223.md` and `parallel-blocks-cluster-223.md`.
The user prioritizes removing unnecessary Domain replication as architectural
preparation for more tasks than Domains, even without large efficiency gains.
Do not reintroduce a universal 100 ms acceptance gate for this work.

The working source removes early provisional inverse publication, excess local
temperature staging, and repeated final native SOL/WAV grid-change writebacks.
It also allows empty source-block stores on ranks without Domains and supplies
those ranks with global checkpoint metadata independently of Domain records.
The latter fixes a demonstrated divergence between checkpoint and timestep
collectives caused by an uninitialized restart clock/output counter. It retains
the distinct saved alias version needed by physics and keeps the startup task
limit until a real above-Domain-count capability test passes. The Git branch
`main` is untouched.

The checked probe then exposed a second empty-rank assumption: shared restart
cleanup returned immediately for `size(grid)==0`, leaving vertical arrays
allocated. That return is removed. This edits `src/main.f90` on the investigation
branch; the Git branch `main` remains untouched. Checkpoint comparisons also
support physical record permutations caused by changed Domain ownership, while
still rejecting overlapping, missing or truncated records.

The next above-limit campaign must test real owned-block work on ranks without
Domains, not just successful launch or idle MPI participation. Preserve exact
topology/threshold/workload checks, numerical field comparisons, adaptation,
remap and restart. A lower task-count performance screen is not a substitute
for this capability test. Any timing comparison must use ordinary-step medians,
exclude initial/restart and saving steps, and use independent run repetitions
only when needed to resolve uncertainty.

Start with 41 ranks on the 40-Domain fixture, then 161 or 167 on the 160-Domain
fixture. Change the startup guard only in the bounded prototype until this
passes. Keep the real distributors and record per-rank Domain/block counts
after redistribution. Test an empty rank zero as well as an empty last rank;
the source of checkpoint global metadata is the owner of Domain zero, which
is not necessarily MPI rank zero.

When preparing the cluster kit, provide complete transfer, build, run, archive
creation/verification and `scp` retrieval instructions **from bbserv to murier**.
Do not require the collector to recreate an existing valid archive. Keep cluster
accounting-disabled and expired `squeue` job-ID behavior in mind.
