# Stage325 verification-only recovery v2 of job 55569

This package verifies the completed captures from failed job 55569. It corrects the source-bound RHS component and unowned transform-dependency contracts in the offline checker only. Job 55571 stopped at rank 1 because the former checker rejected the valid unowned sentinel tuple. The replacement requires the complete tuple: owner=-1, statekey=0, level=-1, mask=0; it does not accept arbitrary negative owners. The original sealed kit, source, binaries, captures, checkpoints, raw work tables and failed job ledger stay unchanged. There are zero compilation and solver-launch paths.

One task on bb13/nk requests one CPU for at most 35 minutes. The site reports MEMORY=1 as a scheduler placeholder, so this package uses the site-compatible --mem=0 and enforces a 16-GiB soft and hard process virtual-memory limit before Python starts. The runner verifies both limits; this is not a 16-GiB Slurm memory reservation. The runner stops at 34 minutes to preserve its report before allocation expiry. Original elapsed 392 seconds plus predecessor verification job 55571's 13 seconds plus the maximum 2100-second recovery allocation is 2505 seconds, below the original 2700-second budget. Queue and collection time are separate. No automatic retry, allocation expansion or performance continuation exists.

Before verification, the runner checks the original kit manifest, source files, 3889 bound returned metadata files, both native executables, checkpoints, private seeds and original work-vector certificates. It requires all three original solver runs' source-bound science, placement and per-run work checks to have passed. Those science checks are reused; this recovery does not independently decode the checkpoints again.

The 768 captures total 2,554,376,142 compressed bytes. Their recorded sizes/timestamps and exact inventory must match the 55569 return. One pre-verification read must match every capture and raw-work SHA256 recorded by job 55571. Those 1,536 historical hashes are included in the sealed package, together with the predecessor failure, elapsed-time accounting, and returned-archive checksum. Changed bytes are rejected even when their size and modification time match. The semantic replay hashes the compressed bytes again while reading them, checks against those historical anchors and emits uncompressed stream digests. Original job 55569 did not save capture hashes. The historical binding begins with the read before job 55571's failed verification; it does not claim hashes existed during original execution. Full original metadata/source/binary/checkpoint identity is checked again afterward.

Semantic coverage, closed RHS phase contracts, provider reconciliation, ordered physical programs/values and all 48 storage-order changes are checked. Progress and failure reports are written atomically in the new recovery output. SQLite scratch is constrained to the new output tree and the checker's 4 GiB limit. A second verifier recomputes work accounting through a new shadow directory: only read-only input links point at old work TSV files; regenerated reports and full vectors stay in the new tree.

Collection streams at most 64 MiB of recovery reports and helper/binding files directly into a compact archive. It excludes original captures/checkpoints, raw TSV and regenerated full peer vectors. It does not create a staging copy or transfer the original evidence. Successful verification remains subject to independent review and establishes no performance benefit.

```bash
bash "$HOME/rk325-verification-recovery-v2/launch.sh"
```

```bash
bash "$HOME/rk325-verification-recovery-v2/status.sh"
bash "$HOME/rk325-verification-recovery-v2/collect.sh"
```

Controller state, intent, logs, and compact archive names use v2 and preserve all v1 records. A new submission requires both source job 55569 and predecessor 55571 to have ended.
