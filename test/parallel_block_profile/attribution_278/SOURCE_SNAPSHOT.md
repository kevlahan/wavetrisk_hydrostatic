# Job 55503 predecessor source snapshot

This commit records the attribution candidate measured by job55503, built on
accepted Stage266 ff740744608731a05d53771c32a3055309f3a586. It is the exact
predecessor of the structural v7 repair, not a production promotion.

The candidate entry in source-files.json binds all251 source/fixture paths.
The original sealed kit SHA256 is
85e4960e80a82bb3373621e1803cfea4f7f0ad308725fec5b219decd1011fd23.
Helpers are tracked for reproducibility; launching needs the sealed kit and J5
fixture described in README.md. Large source/fixture archives and raw outputs
remain outside Git. The four identity manifests also support local helper tests:

```bash
python3 -m unittest test_diagnostic test_attribution test_contract
```

Run that command from this helper directory. The following v7 commit must modify
exactly six solver files; Stage273/277 partition/cache work is excluded.
