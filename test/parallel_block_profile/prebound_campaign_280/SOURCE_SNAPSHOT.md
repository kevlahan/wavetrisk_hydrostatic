# Stage280 campaign source and qualification

The exact five source-family manifests are in source-files.json. Stage266, main, v7 (previous) and Stage279 are byte-for-byte retained from the sealed Stage279 v3/job55507 kit. Candidate is the exact255-file locally qualified Stage280 source; only parallel_block_rhs_exchange.f90 and parallel_block_transform_backend.f90 differ from Stage279. IDENTITY.json records both the two-file repair and29-file total accepted-source overlay. The sixth checked build uses candidate source.

This helper directory alone is not a runnable kit. Sealed runnable source archives, binaries, local logs and compiler evidence are retained in .stage280-validation. The runnable kit has six serial builds and60 solver runs. ENTRY_GATE.json qualifies local campaign entry, not cluster performance or promotion. TEST_VALIDATION.json binds42 campaign helper tests,48 MPI exchange cases, actual local publication rechecks and cross-Python summary verification.

Run helper tests with:

```bash
python3 -m unittest discover -p 'test_*.py'
```

Stage266 remains accepted. Legacy main is untouched. Final lifecycle, granularity/load assignment, total-memory parity, legacy parity and>=640 scaling remain open. No campaign upload/submission or branch merge follows from preparation.
