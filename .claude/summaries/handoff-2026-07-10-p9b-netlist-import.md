# Handoff — Phase 9B liionpack CSV import (2026-07-10)

## Outcome

The missing third parser now exists: bounded liionpack DataFrame CSV input
compiles atomically into `CompiledPackTopology`. This is an implementation
milestone, not P9-G2 closure; sanitizer-instrumented fuzz drivers, corpora, and
CI campaigns remain next.

## Contract and architecture

- Required columns are `desc,node1,node2,value`; reordered columns, RFC-style
  quoted fields, CRLF, BOM, and extra coordinate columns are accepted.
- `V*` rows preserve input cell order and node1-positive orientation. Positive
  `R*` rows are literal branches. Zero-ohm rows union endpoints before sorted
  dense remapping. Exactly one `I*` row is excluded from the graph and defines
  positive/negative terminals.
- Input is bounded to 4 MiB, 100,000 rows, 32 columns, 65,536 bytes per field,
  and 127 bytes per descriptor. Sparse labels up to `UINT32_MAX` never size an
  allocation.
- `PackTopologyInternal.hpp` now provides the single electrical validator used
  by both import finalisation and `SolverWorkspace`; the validator's solver
  semantics were moved without alteration.
- Parsing and exact file reads build only locals. Output publication is a final
  move, and persistent allocation failure preserves the prior topology.

## Evidence

- Stub/red baseline: 3 failed valid/file assertions, 121 hostile assertions
  already green.
- Final Debug and Release: NetlistCsv 728/728, parser allocation 913/913,
  PackTopology 66/66, PackSolver 667/667.
- A resistor-free 2s2p CSV receives the combinator's ladder offsets and cell
  order. The representative resistor graph is independently accepted by
  `SolverWorkspace`.
- No long simulation ran. PackSolver's bounded unit cases were preceded by a
  verbatim validator-move review.

## Explicit limitations

- Literal `Ri*` is added beside the SLIDE cell's model-owned Thevenin
  resistance; no liionpack behavioural parity is claimed.
- V/I magnitudes and source descriptors/node labels are discarded by the
  compiled representation, so lossless CSV round trip is not yet possible.
- The fixed trust budget intentionally excludes ordinary exports of the
  largest 100,000-cell packs.

## Next

Implement raw-byte Experiment, BPX, and NetlistCsv libFuzzer targets. Each must
check deterministic parsing, structural validity on success, and deep sentinel
identity on every rejection. Commit regression corpora/dictionaries and run a
bounded core-only Clang ASan+UBSan campaign in CI.
