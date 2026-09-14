# Cluster review archive: default collection process

User preference established after Stage 191 job 54561: collect all diagnostic
results and logs in one `stageSTAGE-review-JOBID.tar.gz`, preserving directory
names, then transfer it from bbserv to `kevlahan@murier.math.mcmaster.ca:~/Desktop/`.
Use this process for future cluster campaigns. Do not request individual files
unless evidence is missing. Exclude checkpoint payloads and executables; retain
those on the cluster until review completes.

Include results.json, per-run logs and analyses, all rank resource/affinity files,
all raw legacy rank profiles, build identity records, every relevant compiler
log, and the campaign Slurm log. Future kits include `collect_cluster_logs.py`,
which also embeds a file inventory with SHA-256 hashes.

After the job has stopped, run on bbserv (example for Stage 191/job 54561):

```sh
python3 "$kit/collect_cluster_logs.py" --stage 191 --job 54561
scp "$HOME/stage191-review-54561.tar.gz" \
  kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

The helper locates results by their recorded Slurm job ID, not the newest
folder's timestamp. For a custom location or a setup failure before results.json
exists, pass `--results "$out" --tree "$tree"`. It rejects an existing archive
instead of overwriting it. The originally transferred Stage 191 kit predates this
helper; that campaign used the equivalent inline Python collection command.
