# Branch roles and Stage278 continuation

Checked2026-10-14 against actual Git refs/worktrees, following the SourceTree view.
The active continuation is **codex/stage278-structural-repair**, checked out at:

`/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic`

It branches directly from acceptedStage266 ff740744608731a05d53771c32a3055309f3a586.
Its first commit 1b9928a0ed76387adc2f381fa25acf60f87fd4b5 records the exact job55503 predecessor; the following
commit contains the six-file v7 solver repair and qualification/campaign material.
All251 candidate source paths are checked against the sealed55504 manifest.
Commit history records evidence; production acceptance remainsStage266.

| Branch | Role |
|---|---|
| main | Legacy main at6a245db7; branch, remote refs and primary checkout untouched |
| codex/stage278-structural-repair | Active v7 candidate and next repair work |
| codex/domain-migration-measured | Historical Stage267–277 experiments, ending atb5f03c8e; retained with existing uncommitted handoff work |
| codex/performance-investigation | Earlier investigation history with divergent commits; retained |
| codex/stage187-validated | Earlier validated-stage history; retained |
| parallel-blocks | Existing historical branch; retained |

Stage266 is an ancestor of domain-migration-measured. The latter includes
unpromoted partition/cache experiments, so it is not merged into this candidate.
performance-investigation and domain-migration-measured are divergent (53/31
unique commits at inspection); deleting either as though fully merged would be
misleading. stage187-validated is an ancestor of performance-investigation, but
its named historical reference and attached worktree are intentionally preserved.

The cleanup gives the formerly detached Stage278 work a named continuation branch
and keeps large .stage278-validation archives/builds/raw outputs out of Git.
No branch deletion, rebase, reset, force-push or main checkout switch is needed.
SourceTree can show the new branch with All Branches enabled; opening this worktree
shows its candidate checkout. The26d8 worktree stays on its historical branch.

Tracked structural_278 helpers include source identities and compact qualification
reports. Full raw evidence remains locally preserved outside Git. Commit scope is
explicit; unrelated files in the primary or historical worktrees are not staged.
