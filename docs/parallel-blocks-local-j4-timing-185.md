# Local J4/J6 development and timing follow-up

Development can continue locally while bbserv is unavailable. The evolved J4/J6
fixture provides useful correctness, checkpoint/restart and memory evidence.
Local optimized timing can screen changes, but the current measurements cannot
support a precise speedup claim or predict performance with 83 cluster tasks.

## Sequential optimized repeat experiment

Four fresh runs used the frozen Stage 185 seed and input, four MPI ranks,
one thread per library, and all oracles/profiling disabled. Execution order was
Stage 183 baseline, Stage 184 candidate, candidate, baseline. No solver builds
or other simulation runs were started concurrently. All four runs completed
the single checkpoint/reload and following remap; checkpoint 3 matches exactly
across all four, including fields, headers, topology and thresholds.

| Run order | Elapsed (s) | Median ordinary step (s) | System decompression (GiB) | Swapout (MiB) |
| --- | ---: | ---: | ---: | ---: |
| Baseline 1 | 137.55 | 10.80 | 11.64 | 0 |
| Candidate 1 | 127.18 | 9.63 | 13.54 | 0 |
| Candidate 2 | 122.08 | 9.21 | 4.61 | 0 |
| Baseline 2 | 125.92 | 9.60 | 6.34 | 0 |

The baseline changes by 11.63 seconds between repeats and the candidate by
5.10 seconds. Their elapsed-time ranges overlap. Two samples per executable
and changing memory pressure do not establish a precise performance effect.
All four fail the unchanged conservative VM screen, despite zero swapout:
compression/decompression can distort timing without disk swapping. These
counters are system-wide and cannot all be attributed to the solver.

The candidate is promising, and the exact optimized numerical comparisons
support continued development. The apparent timing gain remains provisional.
Do not convert the old 124-versus-120 result or these repeats into an accepted
percentage speedup.

## How to use local timing

- Keep the seed, endpoint, adaptive settings, ranks, build flags and oracle
  settings identical within each comparison.
- Run optimized baseline and candidate sequentially in alternating or balanced
  order on an otherwise idle machine. Repeat enough to establish the spread;
  the four-run experiment here is an initial variability check.
- Compare ordinary step timings and total elapsed time separately. A saved
  checkpoint and reload are part of this fixture's elapsed cost.
- Retain memory monitoring and exact checkpoint gates. Treat changes smaller
  than the observed variation, or results with substantial memory pressure,
  as inconclusive; a passed VM screen alone would not establish repeatability.
- Use clear, repeatable local improvements to guide local kernel and memory
  work. Confirm communication-heavy changes and overall scaling using the
  standard 83-task bbserv configuration when available.

## Remaining checked mismatch

Read-only localization of the existing checked RK4 oracle-on/off checkpoint
failure places all 208 differences at level 4, at 205 distinct horizontal
velocity locations across 31 Domains. There are no differing level-5/6
coefficients. The maximum remains 7.105427357601002e-15. This narrows the next
diagnostic target to the coarse-level velocity path; it does not identify the
first divergent operation or establish causality. No numerical code or
tolerance was changed, and full checked acceptance remains open.

## Evidence

- Raw timing runs: `/private/tmp/wavetrisk-stage185-local-timing-01/`
- Timing driver: `/private/tmp/wavetrisk-stage185-local-timing.py`
- Durable timing identities, gates and summaries:
  `parallel-blocks-local-j4-timing-185-results.json`
- Oracle location summary:
  `parallel-blocks-local-j4-oracle-locations-185-results.json`
- Full diagnostic locations: `/private/tmp/wavetrisk-stage185-oracle-locations.json`

These follow-up files augment the Stage 185 protocol and preserve its failed
checked gate. Production sources and both optimized binaries are unchanged.
