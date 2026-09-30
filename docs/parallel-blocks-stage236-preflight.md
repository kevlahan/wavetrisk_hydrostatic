# Stage236 preflight handoff — RK closure remains unimplemented

## Status

This continuation verified and traced the Stage235 implementation. **No solver
migration was implemented or qualified.** Stage235 remains the latest accepted
numerical result. This handoff is an additional boundary audit, not a replacement
for the final Stage235 handoff and not evidence of completed Stage236 work.

Authoritative directory: `/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.
Branch: `codex/performance-investigation`, HEAD
`4ea640128a71a857f9b8c6b9411dc34fc0aeb641`. Protected main remains
`6a245db7f1632b0205c2bebe95e5f0cd535c2582`. Existing uncommitted implementation
and historical untracked files were preserved. No commits or staging occurred.
The automatically created 26d8 checkout was not used as implementation source.
Write access to the authoritative directory was successfully requested through
normal sandbox escalation; there is no outstanding permission blocker.

## Verified evidence

`parallel-blocks-stage236-preflight.json` records the source-bound audit.
All 69 frozen build inputs match their Stage235 manifest and the live worktree.
All 115 frozen source identities and their separately recorded resolved contents
match; all four frozen executables match. The four live Stage235 helper hashes
also matched the final report. The manifest's digest hashes symlink text for a
symlink and separately hashes resolved contents. A naive contents-only comparison
incorrectly reports differences for test selectors and generated physics links.
No completed Stage235 solver checks were rerun. No new numerical acceptance,
performance measurement, or remote job was performed.

## Joined RK boundary that still needs implementation

The JSON supplies source locations and hashes for the following reviewed chain:

1. `capture_block_domain_multistage_candidate_tendency` retains the current
   tendency before rolling blocks back to the existing timestep base on stages
   after the first. `begin_local_block_accumulated_tendency_trial` forms the next
   candidate and retains its base in `block_tendency_trial`. Reuse this storage.
2. `RK_sub_step_compatibility` independently reads the Domain timestep base and
   fills mass/velocity `q1` or final `q`; temperature still copies scaffold and
   uses `prepare_temperature_boundary_stage` for sparse owner-produced updates.
3. `retain_block_native_multistage_candidate` publishes retained temperature.
   `WT_after_step` then reads Domain scaling through boundary exchange and, on
   the final stage, coarse velocity restriction. Its native scalar callback
   imports those boundary values before native scalar/vector wavelet production.
4. `prepare_native_inverse` seeds noncanonical operands from provisional
   scaling/wavelet fields, but **fixed scaffold comes separately from global
   sol**. On the final stage sol aliases the updated destination; on provisional
   stages it remains the step base. A single cached alias image cannot replace
   these distinct versions correctly.
5. The native inverse still publishes Domain scaling, completes provisional
   all-level aliases, and checks publication counts. Native tendency consumers
   subsequently use the explicit provisional operand epoch. Final consumer
   aliases have a separate deferred publication lifetime.

Consequently, deleting mass/velocity returns or skipping compatibility RK alone
would leave live consumers without their inputs. Moving only numerical interior
updates would also leave temperature closure, forward boundary import, and inverse
scaffold inputs tied to q1/sol. Do not replace this chain with another complete
Domain-layout base or stage snapshot.

## Next implementation and acceptance

Implement a joined retained closure using the existing block base plus explicitly
registered, versioned alias/scaffold inputs. Connect it to forward boundary and
inverse seeding consumers, preserve ordered alias copies and final/provisional
version distinctions, then remove mass/velocity RK returns after their last
consumers move. Update transaction/publication guards to reflect actual ownership,
without weakening independent reference checks. Keep the new selector off by
default and retain the feature-off route.

Qualify fresh frozen sources with exact intermediate oracles, old-operand poison
or absence across the entire RK/forward/inverse cycle, rejected missing and stale
producer reads, and full same-rank 4/41-rank adaptive/remap/restart/checkpoint/VTK
comparisons. Rank40 must execute owned blocks with no Domains. Stage235's poison
test restores primary compatibility fields for RK, so it does not prove this new
boundary absent. Preserve that distinction in any acceptance report.

After RK closure, upstream physics/soil, adaptation, restart, reference aliases
and output still require migration before their numerical storage can be removed.
No speedup or scaling claim is supported. Do not submit remote jobs now. Later
small bbserv correctness work must exclude bb04; nk has two 128-core sockets.
Later retrieval instructions must separately include archive creation/verification
and scp from bbserv to `kevlahan@murier.math.mcmaster.ca:~/Desktop/`.
