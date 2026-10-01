# Stage250 logging calibration — job 55240

Quieting the per-rank storage diagnostic did **not** recover the large J5
regression. At 83 ranks on bb13, four independent full-run medians per mode
averaged **1.661949 s/step with logging** and **1.679132 s/step without**.
The A-B-B-A paired logging savings (loud minus quiet) were −42.657, +17.010,
−29.848 and −13.241 ms/step, averaging −17.184 ms. Their mixed sign and
spread give no evidence of a positive logging saving at this resolution.
Both modes used the same patched candidate binary; quiet mode retained the
complete Domain-storage assertion and suppressed only its write.

The separate main-quiet-quiet-main confirmation measured **0.541002 s/step**
for read-only main and **1.664006 s/step** for the quiet candidate: a
**1.123003 s/step residual gap**, or **3.076 times** main's elapsed time.
The two adjacent-pair regressions were +207.24% and +207.92%. This agrees
with the Stage249 current-binary gap of 1.136932 s within 14 ms. No
promotion or additional Domain removal follows.

All fifteen full J5 runs completed twelve steps, two remaps and 62 VTK
members with the same adaptive trajectory. All candidate runs reproduced the
previously qualified checkpoint hash; main's checkpoint hash was stable.
Twelve cross-version numerical checkpoint gates passed. The loud logs contain
the expected 1,162 storage-absence lines per run, while quiet and main logs
contain none; the quiet candidate still has the expected bounded output
diagnostics. The 83-rank placement used one physical core per rank, 42/41
across the two bb13 sockets. Both candidate modes have the same binary hash.
The independent [reviewer](../test/parallel_block_profile/review_cluster_250.py)
recomputed all fifteen run medians, both paired results and the diagnostic
counts and verified all 2,531 review members. Archive SHA-256:
`7779c3bbe599bad44b423efc18e1c68a7c32f6bbdc0eadb22a26192c9888da33`.
Its [machine-readable result](parallel-blocks-stage250-review.json) records
the individual runs. As with Stage249, the small archive omits checkpoint,
VTK and binary payloads; their recorded hashes and numerical reports were
checked rather than regenerated from omitted payloads.

This calibration isolates the diagnostic **write**, not the storage assertion
scan or other candidate work. The next performance step is a same-source,
same-node J5 feature ladder through the six early owner features that are
enabled in both Stage248 variants. The Stage248 feature-off baseline was
already about 1.629 s/step, versus 0.754893 s for the Stage230 block
candidate; that cross-campaign contrast is a search lead, not a feature-level
attribution. First establish a Stage230-like all-six-off control from the
current source, then enable columns, primitives, diffusion, mass owner,
velocity owner and initial tendency in dependency order. Validate every
variant numerically before timing; confirm any large adjacent jump with
independent A-B-B-A runs. Keep main unchanged and Domain removal paused.
