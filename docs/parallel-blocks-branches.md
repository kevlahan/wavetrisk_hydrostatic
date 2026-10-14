# Branch roles and Stage280 continuation

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
| codex/stage278-structural-repair | Reviewed Stage280 candidate, retaining v7 ancestry and controls |
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

## Stage280 reviewed candidate snapshot

The next commit after de558ed03798f721991f38eeecc00f00eafc8370 records the complete
Stage280 source tested by job55508, including the prerequisite provisional
transform implementation from Stage279 and the two-file construction/publication
repair. These had remained uncommitted together. The combined source changes
12 files versus v7; all255 packaged candidate source paths match local files and
the sealed kit. This is one reviewed end-state snapshot, not promotion of the
rejected intermediate Stage279 implementation. Its exact source remains preserved
in prior kit manifests and archives.

The commit includes source, regression/MPI/campaign helpers and the reviewed
55507/55508 evidence, scorecard and minimum-run campaign policy. Full local
.stage279-validation and .stage280-validation build/return archives are preserved
outside Git. The small job55506 audit archive inside the helper tests is a
required regression fixture. Stage266 remains accepted; main and historical
branches are unchanged. No merge, branch deletion or remote push is required.
