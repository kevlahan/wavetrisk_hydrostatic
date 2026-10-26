# Stage325 full native verification, v5 and v5p1

**Job 55576 passed full native qualification and independent terminal review.**
The supported case is the saved job 55569 source at fixed level 5, `tol=0`,
256 ranks, 30 layers, Simple physics and four steps. This closes the correctness,
provider, physical inverse and work gates for L256. It supports the authorized
verification commit. Performance and the cause of negative scaling remain open.

Job 55576 completed on bb13 with exit `0:0` in 441 scheduler seconds, from
17:08:26 to 17:15:47 on 2026-10-26 (scheduler timestamps). Its mandatory Linux
preflight passed all 292 tests without skips. The six closed phases establish
512 optimized capture reads, 254 provider families per arm, 508,840 matching
features per arm, all 1,024 physical inverse bridges and 768 work reads. The
bridge checks include 1,380 physical node instances whose storage labels moved.
The earlier checked closure is reused narrowly from job 55573.

The 1,864,493-byte review contains 227 manifest-listed payloads plus its manifest,
with 10,080,450 logical bytes. Its SHA256 is
`67e5823b502e5c40d36b37144f71fc7752b2b5132f3d82ef59ae1b5a8864ecd0`.
All archive entries, sealed package bytes, phase dependencies and returned
artifacts were verified. Both retained native database hashes match their cache
certificates. A separate 12.217-second read-only audit freshly checked all 1,536
historical capture/work hashes, 258 sources, 79 original kit files, 3,889 metadata
files, both binaries and all six checkpoints. Proven device/inode aliases bind
the submission's `/home` paths to the runner's canonical `/1/home` paths.

The actual step inherited a 16 GiB soft limit and an unlimited hard limit. The
corrected bootstrap installed and read back exactly 16 GiB for both limits;
the worker inherited and verified them. This observation describes 55576, not
the unrecorded limits of 55575. The supervisor observed no signals or surviving
owned child processes. Sampled peaks were 2,558,805,416 derived bytes and
3,402,747,904 aggregate virtual bytes, below the respective 4 GiB and 16 GiB
guards. Aggregate checks remain sampled limits, not hard aggregate quotas.

Work accounting confirms retained blocks rising from 6 to 400 and remote blocks
falling from 634 to 240 per reconstruction. Serialized remote send bytes fall
from 408,995,792 to 155,278,016, while 17 direct collectives per rank remain.
Replacement work, storage moves, local publication and boundary imports remain
accounted. Counter reductions establish neither elapsed benefit nor the primary
cause of the 160-to-256 slowdown.

Independent review passed with no qualification blockers. The original science,
placement, exact GCC 13.2/MVAPICH 4 builds and checkpoint checks are reused through
their preserved hashes; no solver was rerun. The compact review verifies the
source-reviewed checker's completed evidence chain and native cache hashes; it
does not constitute a second independent execution of every semantic comparison.
Sealed reports retain their original pending qualification fields. The external
terminal review supplies the qualification decision.

The original allocation total remains 2,656 seconds. This separately authorized
attempt used 441 of 1,800 seconds, for a cumulative 3,097 seconds. The unused
1,359 seconds do not authorize another submission. Monitoring and collection are
complete. No retry, resumed run, solver or performance campaign is automatic.
See `test/parallel_block_profile/locality_325/NATIVE_V5P1_TERMINAL_REVIEW.json`,
`validation/independent-review55576.json` and the bounded `validation/native55576`
evidence capsule. All failures, packages, captures and checkpoints are preserved.

## Historical preparation and failed startup

The user authorized full native qualification on2026-10-26. This authorizes one
verification-only bbserv allocation against the preserved job55569captures.
There are zero builds, solver runs and performance runs. Original source,
binaries, captures, checkpoints, packages and earlier failures remain protected.
The previously used2656seconds remain recorded; the new allocation has its own
1800secondcap. No automatic retry, resume or budget extension is permitted.

The exact sealed archive is
`4cac649fa13cf02731dd4cf4fad939e2e40719a75cc79bb93006beb36165aa7e`;
its195filemanifest is
`e62059e17e54e4454396f5a21e8f6b789b2e286262888a6ff38dfc4da97403fb`.
The payload occupies9,203,858logicalbytes. The independent design/package review
is `test/parallel_block_profile/locality_325/NATIVE_V5_DESIGN_REVIEW.json`.
It approves this exact package conditional on288native Linux tests passing with
zero skips, followed by the unchanged controller's mandatory preflight.

Local fresh-package tests passed288cases with six Linux process tests skipped.
The native preflight must exercise those six tests, including real child and
grandchild cleanup, resource stops, repeated signals and final-report retention.
The initial package preflight caught a missing root source fixture before any
submission; its log and preliminary packages are retained separately.

An additional saved-capture audit found that170checked-buildfamily identities
differ from optimized family identities through transform epoch labels. The
adapter now derives and compares the two exact optimized rank0inventories using
historical compressed-byte anchors and the actual source-bound reader. Both
rank0captures and the two retained optimized rank101captures agree on the exact
254familyinventory. No epoch/key normalization or weaker family comparison was
introduced. The checked256stream/254providerclosure remains narrowly reused from
55573. See `validation/native-v5-family-preflight.json`.

## Allocation and evidence gates

The allocation is one task with16CPUs onbb13/nk, exclusive with--mem=0. A Python
verification step uses srun with--cpu-bind=cores and validates its actual16CPU
mask before reading captures. Both address-space limits are hard-set to16GiB.
An external supervisor samples the conservative aggregate process VM and all
derived work/log/report bytes. These sampled aggregate checks are not hard
aggregate quotas. SQLite page limits and bounded writes apply within a4GiB
logical-byte allowance with reserved supervisor/log headroom.

The verifier runs in its own child process group. Its alarm is1740seconds;
the supervisor deadline is1750seconds, inside the1800secondallocation. There is
no early Slurm TERM warning. Cleanup targets only the owned group, uses bounded
TERM/KILL escalation and retains certificates/databases even on interruption.
An unreapable process records a failure without discarding the terminal report.

A successful software gate requires512new optimized stream/read proofs, exact
native cache census for each arm,254provider proofs per arm,508840features per
arm, all1024strict physical inverse bridges,768work read proofs and a final
fresh1536input/source/checker preservation check. All six named certificate
hashes must agree with actual certificate files and both completion receipts.
Independent terminal review additionally verifies the bounded archive, sealed
package identity, dependency/artifact manifests and retained database hashes on
bbserv. No database, capture, checkpoint, raw TSV or large work vector is copied
into the bounded32MiBreview.

Full qualification and the previously authorized conditional commit remain held
until that complete terminal evidence supports them. This gate establishes no
performance gain or cause of the160-to256negative scaling by itself.

## Job55575 terminal review and corrected candidate

All288Linux preflight tests passed with zero skips before submission. Job55575
then failed at the allocation guard on2026-10-26T16:19:05 before any verifier
state, database, capture read or work replay began. Start and end timestamps
coincide and Slurm records elapsed0seconds; the observed Python guard took about
0.000646seconds. This is rounded accounting, not a claim of literally zero
wall time. The failure archive contains203verified entries and has SHA256
`a8a1147a12fed8096e7611ea51f6f0eb5bb4cd40b418f12966ba0c0f306681b2`.

The guard found that the step's inherited soft/hard address-space limits did not
equal the required16GiBpair. Their actual values were not recorded, so the
evidence does not identify which limit differed or how it changed. Slurm's
PropagateResourceLimits setting alone does not establish the mechanism. The
launch incorrectly relied on batch-shell settings being sufficient inside the
new srun step. This failure supplies no new semantic/scaling result.

The corrected v5p1 candidate explicitly installs both limits at the Python CLI
boundary inside the step and reads them back before allocation/capture checks.
It logs actual before/target/after values and installation failures before any
output directory is created. The library supervisor does not irreversibly lower
its caller's hard limit. Native tests include actual isolated Linux installation,
grandchild inheritance and CLI bootstrap before rejection of missing allocation.
Existing exact limit, placement, resource, source, coverage and comparison gates
remain mandatory. Original v5source/package/state/failure evidence are preserved.

The final candidate archive is
`99fe77316934f2991f65ecf2e0c22ee15051f490ccb8f38a27c9fb78bc9f86b7`,
manifest`1a47b001cad28d82de231623f3c85e7b271e744369f9549d8a7b6b30d56dd196`.
Its197files contain9,218,067logicalbytes. All292local package tests passed with
eight Linux-only skips. At this preparation point, a standalone292case Linux
preflight and approval for one corrected submission were still required.

A separate5.179secondread-only postfailure audit reverified all1536historical
capture/work hashes,79originalkitfiles,258sources,3889metadatafiles, original
binaries and six private checkpoint hashes. This audit is recorded separately
from Slurm elapsed time. Original2656seconds and the remaining1800seconds of
the new allocation allowance are unchanged. Full qualification and conditional
commit remain held.

The corrected standalone Linux preflight subsequently passed all292cases with
zero skips (47.693seconds including its wrapper). The sealed manifest remained
unchanged. Independent exact-package/bootstrap review found no repair blockers;
see `NATIVE_V5P1_DESIGN_REVIEW.json` and `NATIVE_V5P1_READY.json`. The corrected
package was uploaded without submission. The user subsequently approved one
corrected30minuteverification-only job; its terminal result is recorded above.
The historical sealed README's causal wording about55575does not establish its
unrecorded limit tuple or the exact Slurm mechanism. Preserve the narrower
finding from the actual failure log.
