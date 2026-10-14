# Qualified structural repair v7

The candidate in source-files.json is the exact251-path source measured in
job55504. Its predecessor is the previous entry. The six-file solver delta and
baseline identities are in IDENTITY.json; no solver content was changed for Git.

The complete sealed kit (including source/fixture archives) is retained locally
outside Git with SHA256
78d70c48b8f89b2329108e0427c4d54c7240f5ffce9fb7c9f7c91b9af1d804c5.
README.md describes its launch workflow. This tracked helper directory alone is
not a complete launch kit. Unit tests are self-contained in this directory:

```bash
python3 -m unittest test_diagnostic test_attribution test_contract
```

Compact qualification evidence lives in qualification/. Full raw logs, binaries,
checkpoints and immutable archives remain in .stage278-validation outside Git.
The entry gate and source checks in the sealed kit remain unchanged. V7 surpassed
266 absolute time on the tested fixture but did not attain main parity or repair
160-to256 scaling; it is a candidate branch, not a production promotion.
