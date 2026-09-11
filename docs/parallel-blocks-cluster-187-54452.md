# Stage 187 cluster result: job 54452

The user supplied the completed bbserv campaign log and full results.json on
2026-09-11. All ten runs report PASS and the campaign status is passed. This supports
repeatable Stage 187 gains on the short cluster restart fixture.

| Metric | Pre-187 baseline | Stage 187 | Reduction |
| --- | ---: | ---: | ---: |
| Median elapsed time | 32.0893 s | 30.5763 s | 4.71% |
| Median ordinary-step time | 1.65 s | 1.52 s | 7.88% |

Three measured runs per version used alternating order A B B A A B. Every
candidate elapsed time (30.4810–30.6625 s) is below every baseline elapsed time
(32.0762–32.1993 s). Warmups and checked validation runs are excluded.

The session configuration is RK4, PARAM=param_J5, max_level=7, 83 ranks,
resume=3, dt_write=0.01 and time_end=0.315. Both versions received the common
checkpoint I/O workaround. The expected new checkpoint is 4, followed by reload
and continued timesteps. The [full campaign record](parallel-blocks-cluster-187-54452-campaign.json)
archives frozen input and executable hashes, placement and detailed gates.
All ten runs completed 12 steps and report a post-restart remap after checkpoint
4. Placement was bb[01-02], with 42 and 41 tasks and one thread per rank.
All 16 field categories, topology, headers and thresholds have zero reported
differences and no nonfinite values. There are eight independent comparisons
and two reference self-checks. Current local production sources match the
candidate and common I/O patch source manifests; cluster executable hashes are
consistent across each build kind in the campaign.

In the supplied runner, PASS follows log/restart coverage validation and the
exact checkpoint gate. Checked-on compares against checked-off; optimized runs
compare against baseline-warmup. The two initial references self-check layout
and finiteness. These gates passed without relaxing numerical tolerances.
They do not establish checked-versus-optimized or legacy equivalence.

The checkpoint I/O workaround successfully completed this solver campaign,
including independent checkpoint payload comparisons. This supersedes the
previous pending cluster validation status for this campaign. It does not
prove that every filesystem or workload is free of I/O failures. Node load and
memory/swap evidence were not included in the pasted log.

Retain Stage 187 and the common I/O workaround on this evidence. These results
measure the short restart case; longer production-run gains remain unmeasured.
Machine-readable transcription and derived percentages are in
[the result record](parallel-blocks-cluster-187-54452-results.json).

The retained pre-187 baseline is commit `6d5c8a3a14cbba71eb53324e9a6f2c7526ec087f`;
Stage 187 is commit `073713bcdd972cc926b1b86247d2f8c9491d1db6`. The common
checkpoint I/O fix, probes, cluster protocol and this evidence are recorded
in the following separate commit. Apply the common I/O fix to both comparison
commits to reproduce the cluster sources. The rejected zero-import experiment
is excluded. Known legacy and checked-versus-optimized discrepancies remain
explicit follow-up work, not closed correctness claims.

Before committing, all 61 helper tests passed on the combined working tree.
The isolated baseline and candidate commit snapshots also passed all 54 and
56 available helper tests respectively, and matched their supplied source
manifests. The snapshots used the original repository only for reading the
pinned historical Git diff required by one legacy-protocol test. Whitespace in frozen input fixtures and
archived patch context is preserved to retain their recorded hashes.
