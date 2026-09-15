# Stage 196: paired legacy/block diagnostic campaign

This stage changes observers only. The baseline is committed Stage 195
`5dfe0db1`, with Stage 193 and Stage 194 selectors explicitly enabled, Stage 192
sparse installation disabled, and contraction off. Legacy retains the fixed
source revision and common checkpoint I/O fix. No merge or optimization patch
must be reconstructed manually.

## Transfer and build

From the local machine:

```sh
scp /private/tmp/wavetrisk-stage196-cluster-kit.tar.gz kevlahan@bbserv:~/
```

On bbserv, use the same compiler/MPI/netCDF/LAPACK/Python environment as Stage 194:

```sh
cd "$HOME"
tar -xzf wavetrisk-stage196-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage196-cluster-kit"
tree="$HOME/wt-stage196"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage196-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage195-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage195-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 5dfe0db1
build_job=$(sbatch --parsable "$kit/build_196.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-196-build-${build_job}.log"
```

Build once and wait for `PASS: five fresh J5 builds ...`. Missing logs while a
job is pending are normal; inspect `squeue -j "$build_job"` before resubmitting.
Both build and test use exclusive nk/bb13. Keep `-march=native` and contraction
off. The five executables are untouched Stage 195, instrumented Stage 195,
checked instrumented Stage 195, untouched legacy with I/O fix, and instrumented
legacy. External physics object hashes must agree across all five builds.

## Run the diagnostic campaign

After the build succeeds:

```sh
cd "$HOME"
source "$HOME/stage196-state.sh"
out="$HOME/stage196-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_196.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$HOME/stage196-state.sh"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-196-${job}.log"
```

The job uses 83 tasks, the same J5/J7 input, checkpoint 3 -> 4/reload and subsequent
remapping. It has 58 runs: two checked runs, eight warmups, and six rounds of
four paired measurement modes. Pair order alternates and mode order rotates.
Expect roughly half an hour or more; the allocation limit is two hours.

- Bare: untouched block and legacy executables; authoritative performance references.
- Minimal: identical diagnostic binaries with common outer phases/WT transactions,
  but operation timers and counters disabled.
- Coarse: boundary/transfer operation scopes only; measure generic Domain handling,
  block ghosts, inverse/scalar/mass transport, writeback, and their posting/waits.
- Detailed: all 77 named operation scopes, with exclusive self wall time attributed
  to dynamics WT, tendency preparation, or the rest of dynamics. This is the new
  operation observer, not the old detailed profiler. Both old profilers remain off.

All diagnostic modes use the same common Stage 195 outer observer in mode 1.
`WAVETRISK_OPERATIONS196=0/1/2` selects minimal/coarse/detailed operation observation.
Uninstrumented references also calibrate disabled-hook overhead. Never interpret
coarse-mode unclassified time as measured arithmetic. Detailed numerical scopes
can include traversal/setup; the report labels these and retains residuals.

The observer never adds MPI barriers. It uses the identical source in both
implementations, buffers all records, and writes after the simulation loop.
Every rank must supply all steps and conserve its phase totals. Overflow,
missing records, wrong modes, inconsistent event flags, or failed numerical gates
stop the campaign. All observer modes must match their untouched reference
checkpoints exactly. The existing legacy numerical screen is unchanged.

## Read and return results

`results.json` contains per-run workflow summaries, repeated paired differences,
profiling overhead relative to minimal mode, phase excess with conservation,
correctness gates, build identity, affinity and memory evidence. Positive signed
phase differences mean blocks took longer. They do not by themselves prove that
the work is avoidable. The per-rank `workflow195-rank-*.txt` files retain the common phase and transaction
envelopes. New `operations196-rank-*.txt` files contain dense operation self times,
call/work counts and sparse logical block-message/byte records for each step and
caller context. Both are included automatically by the collection script. The same physical outer
phases are setup, advance, physics, adaptation/reconstruction, remap and output.

Transaction 1 is the complete WT_after_step envelope. Transaction 2 deliberately
labels implementation-specific tendency preparation: its calls/scopes must be
audited before subtracting its cost across implementations. It is not included
again in the outer phase sum. Operation self times must sum to the common dynamics and WT envelopes for each
rank and step. WT and outside-WT are disjoint; the tendency tag is a non-additive
diagnostic. No globally synchronized arrival times are collected.

After the job stops, including on failure, create the archive before transferring:

```sh
cd "$HOME"
source "$HOME/stage196-state.sh"
python3 "$kit/collect_cluster_logs.py" --stage 196 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage196-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that one archive. Preserve checkpoints and executables until review.
Do not select an optimization based on the largest block-only timer. The first
review must assess observer calibration, pair variation, matched phase excess,
rank timelines and residual uncertainty together. Hardware cache/bandwidth
counters are not collected by this portable campaign; memory causes remain
hypotheses unless supported by additional measurements.
