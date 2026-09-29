# Stage230 handoff: Domain tendency storage removed locally

The user authorized pursuing complete removal of unnecessary Domain numerical
storage and pathways, then asked to resume the work after a usage-limit pause.
The active branch is `codex/performance-investigation`, based on `5c457a5d`.
The user separately advanced `main` to `6a245db7`; this work leaves it untouched.

Seven production source files have uncommitted changes. The qualified restarted
RK4/Simple climate path no longer allocates or extends `trend`. Its native
producers feed RK updates without entering the generic legacy RK consumer.
VTK has short-lived mass-divergence scratch instead of borrowing `trend`.
Independent oracles and other configurations retain their required reference
storage. All-rank initialization checks and timestep guards enforce the policy;
an explicit legacy fallback without storage fails immediately.

See `parallel-blocks-stage230.md` for scope, the remaining-storage ledger and
validation details. `parallel-blocks-stage230-local-results.json` has status
`passed`: optimized 4-rank production, full 4-rank oracles, 41/40 production,
checked 4-rank production, and deliberate forbidden-fallback rejection all pass.
The four positive runs each cover twelve steps, both reloads and remaps, exact
checkpoints and exact VTK output against a reference at the same rank count.
The reference checkpoint is also exact across rank counts.

An initial 41-rank VTK review used a four-rank reference. The unchanged baseline
reproduces that difference because exported vertex lists depend on MPI rank
count. A matched 41-rank reference resolves all 62 VTK members exactly. The
driver now rejects mismatched references before launch; no numerical tolerance
was relaxed and the candidate solver did not need a rerun for this correction.

Old `/private/tmp/stage230-local` artifacts were purged during the pause. Do not
use those empty directories or claim their logs were recovered. New raw evidence,
build trees, an independently regenerated evolved seed, a review patch and hash
manifests are retained in the ignored worktree directory
`.stage230-validation/20260929/`. This directory is outside system temporary
storage. The new fixture hash and pinned source/compiler identities are recorded
in the results. All 23 physics objects match. No timing gain is claimed.

Repeatable case reviewer/runner:
`test/parallel_block_profile/run_local_230.py`.
Real MPI negative test:
`test/parallel_block_profile/reject_fallback_230.py`.
The raw `check-state.json` is `passed`, and `summarize.py` revalidates source,
binary, input and log hashes before refreshing the results document.

The next architectural unit remains the full inverse → shared primitives →
tendencies → RK → forward/inverse cycle on block owners. Domain `sol`/`q1`,
diffusion scratch, geometry-indexed native mass/velocity workspaces, wavelet
fields and their current consumers still exist. Do not call this a complete
block-owned cycle, turn the old Stage218/219 snapshots into production by
renaming them, or resume the superseded narrow inverse-publication timing target.
Use the existing intermediate oracles and zero-Domain coverage as those consumers
move. J5 cluster qualification and any subsequent commit remain outstanding.


Stage230 cluster campaign is now sealed. See `parallel-blocks-cluster-230.md`
for complete transfer/build/run/retrieval instructions and
`parallel-blocks-stage230-kit.json` for identity and packaging tests. The durable
kit is in `.stage230-validation/campaign-kit/`; its archive is also delivered to
Desktop. It uses three J5 builds, correctness at16/167, bare16/83 pairs and83
ABBA confirmation with a3% regression guard and no minimum gain. Same-rank VTK
comparisons are mandatory. The positive path has15 solver runs, maximum19.
Source reproduction,14 packaged unit tests and build/shell workflow simulations
passed; those simulations did not run compilers, submit jobs or transfer files.
Four real local Stage230 logs independently passed the new storage/timing
reviewers. No bbserv jobs have been submitted and no source commit was made.


## Cluster result received: job55089

All15 planned runs passed on bb13. Reviewed archive:
`/Users/kevlahan/Desktop/stage230-review-55089.tar.gz`;2590 inventory hashes pass.
See `parallel-blocks-stage230-review.md/json` and the repeatable independent
reviewer `test/parallel_block_profile/review_cluster_230.py`. Raw extracted data
is retained under `.stage230-validation/review-55089/`.
Checked/oracle16/167 and production167 pass, with exact same-count checkpoints
and62 matched VTK hashes per run. Production has0/N Domain trend allocations;
oracles retainN/N. The83-task ABBA mean is0.7574031835s baseline versus
0.7548926995s candidate (2.51ms/0.33% observed reduction), but paired savings
are+6.60ms and-1.58ms. Accept the architectural cleanup within the3% regression
guard; claim no demonstrated speedup. No legacy rerun or new legacy-gap estimate.
The review changed documentation only plus its reviewer; source edits are still
uncommitted and main is untouched. Recommend commit next, then the complete
block-owned producer/consumer cycle already described above. Do not pursue
another timing campaign for this removed tendency storage.
