# Stage 191: reduce installation copying and expose mass overlap

Validated source snapshot: `37ac536405afbb8b2f8d59037db126950af574cd`. Cluster job 54561 passes all 74 gates, including fresh
full-oracle restart validation. Ordinary-step time improves by 5.96%; ablations
attribute the gain to installation. See `parallel-blocks-cluster-191-54561.md`.
The local evidence and pre-cluster limitations below are retained as history.

## Changes and reasons

The Stage 190 ordinary step was 1.43 s versus legacy's 0.558 s on nk (2.56x).
Its profile attributes about 0.178 s to block construction/migration, including
0.044 s of serialization/copy proofs, and 0.078 s to exposed native mass wait.
These are candidate costs, not guaranteed savings or measured legacy excess.
Even eliminating both whole groups would leave most of the 0.872 s gap. This
candidate tests two concrete causes while adding the missing paired legacy data.

- Installation transfers all 36 owned allocatable components from temporary
  source/receive blocks into final local blocks. It retains all 13 identity/layout
  metadata fields and structural checks. Production source, receive and copy
  serialization proofs are disabled together. Dynamics/adaptation oracle modes
  retain copies and proofs, including the source arrays needed by full validation.
- Mass exchange now has begin/finish operations. The fine-level tendency is posted
  before independent velocity-source work, then completed before its coarse-level
  consumer. MPI requests and buffers persist in the plan; pending buffer reuse,
  mismatched layer completion and plan retirement with requests pending fail.
  Full dynamics oracle mode retains synchronous exchange. Payloads and arithmetic
  are unchanged.
- Each change has an independent runtime selector, checked for rank agreement once
  and logged. Both default on. The cluster campaign measures Stage 190, the full
  candidate and legacy directly, plus installation-only/overlap-only variants and
  a both-off control. Production profiles verify zero proof/copied bytes and
  positive moved bytes.
- A separate frozen legacy timer patch adds phase accounting without replacing
  the authoritative legacy timing executable. Its checkpoint must match legacy
  exactly. All rank records are required, with exclusive time conservation checks.

The installation counter reports packed payload-equivalent bytes whose ownership
was transferred, not measured DRAM traffic or permanent memory saved. MPI packing
and Domain-to-block capture still exist. Post-install source-patch reads found
while investigating lifetime belong to full ghost-source validation, not the
writeback plan as the earlier proposal said; retaining copies in validation
preserves those reads.

## Local evidence

Fresh GNU Fortran 15.2 J4/J6 RK4 builds use contraction off. Optimized and checked
runs each complete 12 steps, checkpoint 2 -> 3/reload and two remaps. Their saved
checkpoint matches Stage 190 and the saved independent Stage 186 RK4 oracle
exactly. The fixed legacy numerical screen passes without changed tolerances.
The optimized detailed run transfers 8,611,644,132 payload-equivalent bytes with
zero source/receive/copy proof bytes and zero installed copied bytes.

The instrumented legacy J4 run completes the same 12 steps/restart, matches the
unprofiled legacy checkpoint exactly, and supplies complete 65-region profiles
for all four ranks. Actual frozen legacy timer patch application and all source
hashes are verified before its fresh local build.

The helper suite has 96 tests (one MPI test skipped inside the sandbox); separate
MPI-enabled component execution passes both tests. Component tests use the actual
Fortran move/exchange routines at O0 and O2. They cover all owned allocations and
nondefault array bounds, delayed rank arrival, repeated alternating layers,
pending reuse and incorrect completion. A fresh additional full-oracle run was
stopped after initial steps because of its local cost; it is not a passing full
oracle result. Full oracle/restart remains required on the cluster.

Local VM measurements invalidate speed comparisons: the optimized run incurred
about 44 GB of system-wide decompression. Do not infer reduced elapsed time or
RSS from this campaign. Detailed identities and comparisons are preserved in
`parallel-blocks-stage191-local-results.json`.

## Larger changes needed to close the gap

1. Retain storage for surviving subtrees during adaptation; then reuse stencil and
   routing plans only when their complete dependency signatures remain valid.
   Refresh every field even when geometry is reused. Use geometric subtree keys,
   not compacted local indices. This component of the original proposal is not
   implemented in Stage 191.
2. Move a complete shared producer/consumer chain onto final block owners and
   remove its Domain-to-block-to-Domain transfers. Choose the first chain from
   the paired profiles. Keeping both complete producers would defeat the change.
3. Test field-oriented scalar storage to remove repeated strided record gathers,
   then revisit work-weighted ownership after separating fixed Domain work from
   movable block work. Avoid another arithmetic-only restriction stage while
   reconstruction and representation transfers dominate.

Use the paired profiles and ablations to order these changes. Preserve direct
previous-stage/legacy comparisons and unchanged numerical gates throughout.
See `parallel-blocks-cluster-191.md` for the five-build, 83-rank nk protocol.
