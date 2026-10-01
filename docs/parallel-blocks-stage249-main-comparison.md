# Stage249 direct comparison with legacy main

The last fresh direct J5 comparison was Stage229 job 54737. Its separate
unobserved 83-rank runs averaged 0.758268 s per ordinary step for the then
Stage203 block candidate and 0.542883 s for legacy: a 215.386 ms gap. The
legacy observer exceeded its overhead-calibration limit, so the phase
attribution in that campaign is provisional. The uninstrumented gap is valid
historical context, not a measurement of the Stage247/248 RK migration.

This campaign compares the current qualified mirror-free migration, sealed as
`bbserv-candidate5.tar.gz`, directly with a read-only export of current `main`
(`6a245db7f1632b0205c2bebe95e5f0cd535c2582`). Numerical solver sources
under `src`, `Makefile` and `data` are unchanged between that main commit and
Stage229's pinned legacy source. As in Stage229, the main export applies only
`legacy_checkpoint_io_190.patch` so the established restarted fixture can be
used reliably on bbserv. This modifies checkpoint I/O, not ordinary-step
numerical kernels; ordinary timing excludes saving/restart steps. The build
shares the candidate's 23 external Simple-physics objects between both solver
variants and uses the same compiler and FP contraction setting on bb13.

The kit validates full 12-step J5 runs at 16 and 83 ranks, two remaps, both VTK
events, same-version exact repeated outputs and the unchanged cross-version
checkpoint numerical policy. The candidate must reproduce the prior v7 J5
state and checkpoint hashes. At 83 ranks, it runs an independent A-B-B-A
confirmation, adding at most one batch if the 3% guard is ambiguous. Timing is
the median of the nine ordinary steps per run; the two or four independent
runs per variant supply the comparison. This protocol measures the current
main-relative gap without claiming a phase-level breakdown or automatic
promotion. It is scoped to restarted compressible J5 RK4 with Simple physics,
not every solver mode or rank count.

The local machine has four usable cores. Packaging and build validation run
serially (`make -j1`); cluster timing must run on bb13 through partition `nk`.
The [kit README](../test/parallel_block_profile/main_comparison_249/README.md)
has the exact transfer, launch and collection commands.
