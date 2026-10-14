# Resident provisional transform candidate

Base is job55504/v7, commit de558ed03798f721991f38eeecc00f00eafc8370.
The exact candidate and unchanged Stage266/main/v7 control manifests are in
source-files.json. The overlay is source-bound to local qualification and the
sealed campaign entry gate. This helper directory alone is not a runnable kit.

Run helper qualification with:

```bash
python3 -m unittest discover -p 'test_*.py'
```

Full local logs, binaries, checkpoints, fixture and sealed kit remain outside
Git in .stage279-validation. The source scope is the provisional transform only;
final lifecycle closure, finer ownership, legacy parity and >=640 scaling remain
open. No branch merge or production promotion follows from preparing this kit.
