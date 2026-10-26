# Stage325: exact 256-owner locality qualification

This is a correctness and changed-work gate. It schedules **two serial full solver builds and three 256-rank solver launches**, each initialization plus four steps. It schedules **zero performance runs**. The optimized control and candidate share the same binary. The checked candidate enables the existing numerical/reference checks.

The exact GCC 13.2.0 C/Fortran and MVAPICH 4.0 identities, native runtime/linkage and physics-object checks pass before the first solver. Two small gate executables are also compiled; 17 callback cases, one MPI identity launch and one 256-rank placement probe precede the three solver launches. The chosen placement must equal the source-bound prior bb13 placement.

The kit reuses completed job 55567 qualification of the unaffected job 55561 source, the private seed originally created using resume=-1, and the independent legacy four-step endpoint/trajectory. It verifies those bytes through the retained sealed review and existing binding helper, without rereading old raw captures. Every new run receives a separate read-only byte copy of that seed; no user checkpoint is searched or overwritten. Missing or changed baseline artifacts stop the campaign. There is no fallback rebuild or new bootstrap.

The three runs are checked-L256, qual-C0 (optimized override off) and qual-L256 (optimized exact locality map). The new control must match its qualified predecessor exactly. Both optimized arms must match each other exactly at initialization/final checkpoint. All three must satisfy the independent legacy trajectory and science gate. Complete compact semantic coverage, ordered physical program/provider closure and the changed-work ledger must pass; a checksum-only endpoint agreement is insufficient.

The new compact mode preserves integer plan records and encodes real records with source-bound SHA256 rather than retaining full floating payloads. It is qualification instrumentation and its elapsed times cannot be interpreted as performance. The source enforces a 128 MiB physical stream cap per rank. The runner monitors compressed output at 64 MiB per rank, 2 GiB per run and 3 GiB for the campaign, checking every half-second and before/after each run. These aggregate thresholds are not atomic filesystem quotas: partial output may overshoot during polling/termination, while source per-rank sink limits independently bound writes. Captures, checkpoints, raw work TSV files and full peer-vector files remain native. Bound summaries and compact per-rank review counters are returned; this does not independently reparse the omitted raw files on the review machine. Collection directly streams a bounded metadata archive, with no staging tree and no recursive raw-data copy; selected files are limited to 64 MiB in total. Originals are preserved if the cap rejects collection.

Planning estimate 15–30 minutes; hard allocation cap 45 minutes. Native compact encoding/checking cost is not measured, so the estimate is uncertain. Timeout or any failed gate produces partial evidence for review; there is no automatic retry, expanded capture or timing continuation.

On bbserv, from the extracted kit:

```bash
python3 check_kit.py &&
bash launch.sh
```

After the job terminates:

```bash
bash "$HOME/rk325-locality-qualification-kit-v1/status.sh"
bash "$HOME/rk325-locality-qualification-kit-v1/collect.sh"
```

The collector retains `rk325-review-JOBID.tar.gz` on bbserv and attempts its normal Desktop transfer. If transfer fails, rerun only collect.sh. Do not resubmit. A successful native result is still `HOLD_PERFORMANCE_PENDING_INDEPENDENT_QUALIFICATION_REVIEW` and does not establish the primary cause or a speedup.
