# Stage325 independent provider and accounting review

Final verdict: no unresolved concrete blocker found in the reviewed paths.
This is an offline source and validator review, not a native 256-rank
qualification result or evidence of a performance improvement.

The review followed actual source paths for fixed-Domain providers, block
storage and patch addressing, numerical-group permutation, ordered ghost
installs, transform publications, inverse routes, and complete writeback
coverage. It checked the added compact record kinds and shape extents against
the preserved coverage validator, and checked the work ledger against the
changed live routes. Root reports the final source reproduces exactly and the
checked compilation passes; those checks were not rerun by this subtask.

Concrete gaps identified during review are resolved:

- Local transform publication now requires a positive canonical row, the
  matching physical key and an authoritative owned row. A matching key alone
  could have admitted a nonauthoritative local provider.
- Writeback now requires the exact complete remote catalogue lists, retained
  count and installed count derived from the selected map. Checking only
  records present on the wire could have missed an omitted block.
- The live fixed-Domain-to-block boundary importer now has a physical provider
  certificate covering all 640 installed blocks, including retained direct
  imports. It checks both directions of boundary metadata, source and level,
  exact node/value extents, ordering, displacements, peer capacities and local
  storage. Boundary scalar/vector and retained local work also have dedicated
  live ledger hooks with the original component and verification guards.
- Extra compact real records explicitly require real kind 3 and the full shape
  product; physical-key records explicitly require integer kind 1. Otherwise
  a common malformed or wrongly typed addition could evade the original
  validator after additions were stripped for compatibility.

The integrated boundary helper is byte-identical to the reviewed helper:
SHA256 `0d675429b55d7f0c382d8c2eff787903e707ba2fafba9504ae57775ddf2e1705`.
Its focused fixtures cover both exact maps and eleven corruptions of physical
metadata, capacities, widths, coverage and ordering. The semantics agent
reports all 30 integrated semantic tests pass, including wrongly typed real
and key additions. The five accounting test methods passed separately,
including parser corruptions, cross-rank byte conservation, both exact maps,
and a coordinated retained-source extent mismatch across qualification runs.

Source checks found no further concrete mismatch in one-based block/patch
addresses versus zero-based node offsets, dynamic patch widths, seven-field
boundary storage rows, actual installed block order, fixed-Domain source rank
identity, or the current numerical-group permutation. Added inverse node-map
and component-order scopes are included in fixed-program comparisons.

The remaining gate is the bounded target-hardware campaign: checked L256 and
optimized C0/L256 initialization plus four steps, with complete native
semantic, work, placement and science qualification. A local initialization
smoke test cannot replace that gate. Work counters and locality changes do not
identify elapsed-time attribution or establish useful scaling recovery.
