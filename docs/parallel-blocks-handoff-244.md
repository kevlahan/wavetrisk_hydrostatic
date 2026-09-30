# Native adaptation qualification

Stage 244 completes native numerical adaptation and topology reconstruction. The frozen candidate is `.stage244-validation/candidate3`; the compact source-bound manifest is `parallel-blocks-stage244-results.json`.

All ten cases passed: full 12-step production trajectories on 4 and 41 ranks, exact adaptation and independent RK/physics oracles, legacy control, epoch and generation rejection, full legacy SOL/WAV poisoning, and native owner readback on 4 and 41 ranks. Each full trajectory includes two remaps, 62 exact VTK members and checkpoint numerical validation. Fifteen arithmetic kernels passed the source audit.

Retired boundary storage must remain represented in native alias metadata because later topology closures can reuse those addresses. No missing-old-node fallback was introduced. Failed attempts remain preserved.

Remaining work: native between-step aliases/diagnostics/remap lifecycle, direct RK coefficient handoff, bounded output/restart adapters, and removal of persistent production mirrors. This is not yet the final bbserv performance candidate. Main is unchanged; no remote performance run was started. Local builds and solver runs were sequential.
