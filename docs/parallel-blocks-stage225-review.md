# Stage225 cluster review: job 54729 passed

All eleven scheduled correctness runs passed on bb13. The complete planned
matrix is now established for the restarted J4/J5 resident block pathway:
41 tasks over 40 Domains and 167 tasks over 160 Domains, including checked
builds, optimized builds and the existing numerical oracles. This establishes
the central capability sought by the block refactor in these tested cases.

Verified all 3,052 inventoried archive files. Six cluster build source manifests
match the outgoing Stage225 kit, their build-log hashes match, and the recorded
binaries match the identities used by each run. Build/run host identity and
physical-core placement were checked at 40, 41, 160 and 167 ranks across both
128-core sockets. Git main is unchanged.

## Raw evidence review

Every run completed 12 steps, two checkpoint reloads and two vertical remaps.
The adaptive trajectories match within each fixture. Raw per-rank lifecycle
records were reparsed and agree with the campaign report: every candidate rank
completed block dynamics and steps, all ranks agree on global metadata, and
Domain ownership covers the mesh at every recorded event.

| Case | Zero-Domain ranks | Blocks owned by those ranks | Result |
| --- | --- | --- | --- |
| J4, 41 ranks, native | 40 | 1–4 | Checked, optimized and oracle runs pass |
| J4, 41 ranks, empty root | 0 | 3 | Checked run passes; checkpoint source is rank 1 |
| J5, 167 ranks, native | 160–166 | 1–6 | Checked, optimized and oracle runs pass |

Both oracle runs explicitly enabled dynamics, adaptation and vertical-remap
validation. Logs confirm completion of both remap oracle checks. The Stage224
empty-rank extent failure does not recur.

All 17 numerical comparison reports pass exact comparison, including
cross-task-count comparisons. Reports show zero differing field values and no
nonfinite values. Each fixture's recorded final checkpoint hash is identical
across all its runs:

- J4: `270a3d6aec3ce3bc58c19fcad8f1c932adeeeb69d0e5674fa89dd234c93f25ab`
- J5: `1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53`

The review archive deliberately excludes raw checkpoints and binaries. Their
recorded hashes and detailed numerical reports were checked; checkpoint bytes
were compared on bbserv, not compared again locally. No new solver runs were
needed for this review.

## Recommended next step

Consolidate the validated cleanup and the two correctness fixes on the
investigation branch. Replace the campaign-only task-limit bypass with a
maintainable, explicitly supported resident-block configuration, retaining the
restriction for execution paths not yet validated. Remove campaign-only
per-rank file output from the production implementation while preserving these
regressions and evidence. Validate that integration at the existing task-count
boundary before expanding scope.

The ordinary production source still has its original startup guard; this
review has not changed it or committed/merged anything. Success does not yet
cover cold starts, arbitrary physics/configurations, multi-node execution,
arbitrary higher task counts, or ranks with zero blocks. These cases all
required useful block ownership on every rank.

There is no speedup claim from this campaign: checkpoint work, checked builds,
oracles and lifecycle file output were intentional correctness work. After
integration, use an uninstrumented production comparison to assess scaling;
exclude first restart steps and saving, use ordinary-step medians, and follow
the agreed screening/confirmation policy. Further broad instrumentation is
not required to establish the capability already demonstrated here.
