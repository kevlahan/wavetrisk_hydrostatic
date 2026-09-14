# Stage 191 cluster campaign 54561

The campaign passes all 74 recorded gates. Stage 191 improves the ordinary-step
median by 5.96% and elapsed time by 3.22% versus Stage 190. Installation-only
reproduces the gain; mass overlap alone shows no measurable benefit in this case.
The ordinary-step ratio to legacy drops from 2.558x to 2.406x. This removes about
9.8% of the excess time above legacy, so most of the architectural gap remains.

## Correctness and provenance

All 33 runs complete 12 steps, reload checkpoints 3 and 4, and remap after
checkpoint 4. The 74 gates comprise 35 exact checkpoint reports (three reference
self-checks), 31 workload comparisons, and eight fixed legacy numerical screens.
Checked oracle-on/off, checked/optimized, Stage 190/191, all ablations and
instrumented/uninstrumented legacy pass their required exact comparisons.
The review additionally re-parses all solver logs and verifies the same printed
workload across every run, including the two checked runs.

Legacy remains numerically screened rather than bitwise equal to block code.
The maximum absolute differences are velocity 3.72324e-12, mass 2.50111e-12,
and temperature 7.85803e-10, within the unchanged 1e-10/1e-10/1e-7 limits.
Soil fields, thresholds, headers and topology remain exact. Checkpoint payloads
were deliberately excluded from the archive: this review validates the supplied
comparison reports rather than re-running comparisons on checkpoint bytes.

The separately supplied build record equals the complete embedded record.
Candidate source hashes match the current local implementation. Every run's
recorded binary hash matches its build role. All five compiler-log hashes match
and the expanded Fortran compilation commands use contraction off: O2/native
for optimized builds and O0/runtime checks for checked. The build record confirms
23 common external physics object identities, J5 and RK4, compiled on bb13.
No compiler settings or solver source were changed during this review.

## Repeatability and attribution

Six timing triples rotate the three implementations within one exclusive nk
allocation; warmups, checks, controls, ablations and profiles are excluded from
headline medians.

| Implementation | Ordinary-step median | Elapsed median |
| --- | ---: | ---: |
| Stage 190 | 1.425 s | 32.88993 s |
| Stage 191: both changes | 1.340 s | 31.83032 s |
| Legacy Domain + common I/O fix | 0.557 s | 22.50012 s |
| Installation only (three ablations) | 1.340 s | 31.79020 s |
| Mass overlap only (three ablations) | 1.430 s | 32.98929 s |

All six candidate elapsed times (31.756--31.860 s) are below all six baseline
times (32.823--33.106 s). All candidate ordinary-step medians are 1.34 s versus
baseline 1.42--1.43 s. The both-off candidate control returns 1.43 s and 33.12230 s.
Ablations follow the main timing triples, so they are supporting attribution,
not equally balanced estimates of tiny differences. There is no evidence that
overlap adds to the installation gain, or a basis to claim a universal regression.

Recommendation: retain the installation optimization. Keep overlap explicitly
experimental or disable it by default in a subsequent reviewed change; do not
advertise it as a validated performance gain. This review leaves the tested
candidate and its selectors unchanged. No commit was made in this review.

## Profiles and memory

Each profile supplies all 83 ranks and 12 steps. Excluding the checkpoint step
leaves 11 ordinary steps. Exclusive timing conservation errors are at most
1.64 microseconds in the block profiles. Re-parsing all raw legacy files exactly
reproduces the embedded legacy profile summary.

Construction/migration rank-mean exclusive time falls from 0.17798 to 0.10104 s
per ordinary step. This 0.07694 s reduction accounts for nearly all the 0.07605 s
net change in the profiled rank-mean step, consistent with the timing ablations.
Production eliminates all three proof-byte counters and installed copy bytes;
13,862,170,508 packed-payload-equivalent bytes transfer allocation ownership.
These are ownership counters, not measured DRAM bytes or permanent memory saved.

Other candidate groups remain substantial: inverse and its transfers 0.25254 s,
scalar input/restriction 0.19997 s, ghost handling 0.08374 s and publication
0.04899 s. Native mass wait changes only from 0.07830 to 0.07596 s, while the
whole mass/velocity group is effectively unchanged (0.09993 versus 0.10040 s).

Legacy's rank-mean profiled step is 0.60910 s versus candidate 1.33376 s. Do not
substitute these instrumented means for the unprofiled medians. Legacy places
substantial inverse/scalar communication inside Domain boundary accounting:
legacy boundary handling is 0.44850 s versus candidate 0.28862 s. Therefore
subtracting similarly named inverse/scalar rows would exaggerate their excess.
The disjoint region table is preserved for investigation, with this limitation.

All runs use the same 83 distinct cores on bb13, with zero node swap-in/out
increments and abundant available memory. Major faults are recorded; this is
not a zero-page-fault claim. Timing-run maximum rank RSS is about 398--401 MiB
for Stage 190, 379--382 MiB for Stage 191 and 164--169 MiB for legacy. The median
sum of rank peaks decreases from 23.010 to 21.681 GiB (5.78%); legacy is 11.021 GiB.
These sums are not simultaneous node peaks and exclude waiting Python parents.

## Next work and evidence collection

The remaining 0.783 s ordinary-step excess calls for larger ownership/lifetime
changes. Prioritize eliminating a complete Domain/block producer-consumer round
trip, using the paired profiles to select the chain. Reuse surviving block
storage/plans with full dependency invalidation as a separate experiment. The
remaining construction/migration group is only about 0.101 s; it cannot alone
close the remaining gap. Keep previous-stage and legacy comparisons throughout.

The single-archive collection process succeeded: 5,702 files, including all
profiles, resource records and five compiler logs, arrived without name
collisions. Per the user's request, use it for future campaigns. The reusable
`collect_cluster_logs.py` adds a SHA-256 inventory, identifies results by Slurm
job ID and excludes checkpoints/executables. Transfer to
`kevlahan@murier.math.mcmaster.ca:~/Desktop/`.

Evidence is preserved in `parallel-blocks-cluster-191-54561-campaign.json` and
`parallel-blocks-cluster-191-54561-analysis.json`. Reproduce the review using
`analyze_stage191_campaign.py` against the extracted archive. Raw files remain
in `/private/tmp/stage191-review-54561/`; archive SHA-256 is in the analysis JSON.
