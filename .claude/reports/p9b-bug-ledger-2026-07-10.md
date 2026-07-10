# Phase 9B adversarial bug ledger — 2026-07-10

This is the live evidence ledger required by P9-G4. Every confirmed defect must
have a short regression written before its production fix. Refuted candidates
remain recorded so a later pass does not revive them without new evidence.

## Confirmed defects

| ID | Severity | Subsystem | Pre-fix red evidence | Fix | Post-fix evidence | State |
|----|----------|-----------|----------------------|-----|-------------------|-------|
| P9-B01 | High | `CyclerV2` transaction | The registered throwing post-advance event callback returned `Invalid_parameters` after mutating the arena; the initial `[P9]` run contributed state/output mismatches to 13 failed assertions. | Restore the full arena on integrator, observable, control, indicator, bisection-trial, and final-event failures before any event-aware step is committed. | Debug and fast-math Release `[P9]`: 42/42 assertions, including byte-exact rollback; full Experiment binary 183/183 in each build. | FIXED |
| P9-B02 | Medium | Experiment validation/parser | Direct segments accepted invalid enum and NaN/Inf metadata; scaled durations could overflow after the only finiteness check. Tests were registered before the validator was added. | Exhaustively validate public descriptors before output/state initialization; reject invalid integrator values; recheck the duration after unit scaling. | Debug and Release direct-validation/event suite 42/42; documented grammar/atomicity suite 56/56; full subsystem 183/183 in each build. | FIXED |
| P9-B03 | High | `Recorder` thin ordering | With capacity one, steps 0, 2, duplicate 2, and backward 1 all returned success; the thinning count reached 3 instead of remaining 1. Three registered assertions failed. | Track the last accepted cadence point independently of stored slots, update it before both store and thin outcomes, and reset it on `clear()`. | Focused Debug/Release 20/20; full Recorder 98/98 in each build. | FIXED |
| P9-B04 | Medium | `Recorder` derived metadata | Bit-injected NaN elapsed time and `1e300 A / 1e-300 m²` both returned success and consumed a slot; four registered assertions failed. | Validate elapsed time and each derived density with the integer-barrier `is_finite` guard before advancing the watermark or writing storage. | Focused fast-math Release plus Debug 20/20; full Recorder 98/98 each. | FIXED |
| P9-B05 | Low | `RecorderConfig` enum | `static_cast<BackpressurePolicy>(255)` returned success and marked the recorder configured; two registered assertions failed after B03/B04 were already green. | Exhaustively accept only `stop` or `thin` before allocations/commit. | Focused Debug/Release 20/20; configuration remains atomic. | FIXED |
| P9-B06 | High | `CompressedRecording` allocation bound | A valid 64-byte header claiming `UINT64_MAX` one-value snapshots threw `length_error`; `CHECK_NOTHROW` and status assertions failed. | Bound snapshot count by `(file_size - file_header)/block_header` before any vector allocation and map `length_error` to `Invalid_parameters`. | Focused Debug/Release 67/67; no exception and failed open leaves the reader invalid. | FIXED |
| P9-B07 | Medium | compressed codec wire value | CRC-valid file and block codec fields of `256` narrowed to `CompressionCodec::none`; open returned success and two assertions failed. | Accept only the exact 32-bit wire values for `none`/`zstd` before casting; availability remains a second check. | Focused Debug/Release 67/67; full async 233/233 each. | FIXED |
| P9-B08 | Medium | `AsyncRecorder` derived density | Finite current/area whose quotient overflowed returned success and produced one written snapshot; two assertions failed. | Validate the quotient with `is_finite` before locking or publishing the slot. | Invalid enqueue, successful empty finish, zero snapshots in both builds. | FIXED |
| P9-B09 | Low | `AsyncRecorderConfig` enum | Policy value 255 configured successfully and left the object configured; two assertions failed. | Exhaustively accept `block` or `thin` before allocations/file creation. | Focused Debug/Release 67/67; configuration is unchanged. | FIXED |
| P9-B10 | High | `AsyncRecorder::configure` rollback | Faulting the measured long-path copy returned `Numerical_failure` but left `configured()==true`; one registered assertion failed. The member publish also preceded context ownership transfer. | Copy the path into a local before all resource creation; swap only during the no-throw member-commit phase. | Allocation executable Debug/Release 15/15; failed object remains unconfigured and then configures/finishes successfully. | FIXED |
| P9-B11 | High | Mode C convergence | Tiny `α` made current updates fall below tolerance while KCL remained violated; the old solver returned `Success`. Terminal-only KCL also misses the registered two-layer internal-node construction. | Require both max current delta and max KCL residual over every non-reference node to meet caller tolerance before publication. | Tiny-gain and terminal-satisfied/internal-violated regressions reject atomically; ModeC 52/52 each build. | FIXED |
| P9-B12 | Medium | Mode C constraint diagnostic | The former `initial_drift·(1-α)^k` field was presented for general coupled graphs although that contraction is topology-dependent; the replacement test was initially tautological with acceptance. | D-28 limits the analytic contraction to P4-G2's one-unknown network and derives it independently; arbitrary graphs gate on measured all-node KCL. Roundoff is labelled an estimate. | P4-G2 clears the independent oracle in Debug/Release. | FIXED |
| P9-B13 | High | solver finite publication | Finite layer values in a two-series ladder sum to infinity while both current deltas remain exactly zero; disabling the new ladder/publication guards makes the registered test publish `Success`. The initial one-cell candidate was rejected as non-discriminating. | Bit-safe finite checks cover residuals, stamps, Eigen outputs, node/current updates, layer accumulation, and final candidate publication. | Mutation-red ladder test; final PackSolver 667/667 in each build. | FIXED |
| P9-B14 | Low | `PackSolveMode` discriminant | Value 255 fell through to relaxation and could return `Success`, changing the public solution. | Exhaustively validate the mode before diagnostics or candidate mutation. | Invalid-mode regression returns `Invalid_parameters`; prior solution unchanged. | FIXED |
| P9-B15 | High | compiled topology/netlist trust | Invalid `PackNodeKind`, branch kind/node/cell indices, false connectivity, sparsity, batch/lane bijections, and ladder tables could reach indexed solver paths. | Exhaustive cold validation, graph traversal, exact sparsity/cell bijections, dense lane proof, and ladder orientation proof before commit. | Topology 50/50 and solver malformed-netlist regressions in both builds. | FIXED |
| P9-B16 | Low | cross-mode diagnostics | A ladder solve retained sparse residual `3.0` and stale constraint fields. | Reset the complete diagnostic value on every valid solve and repopulate only mode-owned fields. | Mode-switch regression observes zero non-owned fields. | FIXED |
| P9-B17 | High | mutable workspace ownership | Non-const `workspace()` allowed move-out/reconfigure after solver configuration, leaving a null or topology-mismatched internal workspace and possible UB. | Expose only `const SolverWorkspace&`; `PackSolver::invalidate()` remains the mutation seam. | Compile-time contract was red before removal and now enforces const ownership. | FIXED |
| P9-B18 | Medium | hostile solver dimensions | `adjacency(node_count)` occurred before proving a connected graph could exist for the branch count; signed/Eigen index bounds were not enforced. | Reject node counts beyond signed index range or `branches+1` before size-proportional allocation. | Hostile `INT_MAX+1` reconfigure rejects atomically and preserves the prior solution. | FIXED |
| P9-B19 | Low | successful reconfigure state | Reconfigure zeroed solution vectors but retained the old terminal-voltage scalar, exposing a mixed stale solution while `has_solution=false`. | Reset terminal and candidate scalar diagnostics during successful configuration. | Reconfigure regression observes all-zero unpublished solution state. | FIXED |
| P9-B20 | Medium | public byte shuffle | Width zero terminated the Release test process; undersized or overlapping spans relied on Debug-only assertions and admitted divide-by-zero/OOB UB. | Return and propagate `Status`; validate width, equal/divisible sizes, and non-overlap before indexing. | Full AsyncRecorder Debug/Release 238/238; bitwise involution retained. | FIXED |
| P9-B21 | High | compiled-curve indexing/fast-math | The test-first Debug `[P9]` run failed 6/21 assertions: a denormal domain reached an infinite reciprocal, extreme finite ordinates produced non-finite slopes, and NaN derivative lookup reached the float-to-index path. The first guarded implementation then access-violated in Release because ThinLTO replaced a constexpr-NaN return branch with an `llvm.assume`. | Prove spacing, slope, range, bin ratio, reciprocal, sample, error, and scaled index representable before conversion/commit. Classify query bits through a reference boundary and construct the NaN sentinel opaquely so `-ffinite-math-only` cannot delete the guard. | Full CompiledCurve Debug and fast-math Release: 163/163 assertions each, including raw IEEE sentinel bits, invalid tolerance, BPX resolution, and bit-exact legacy OCV interpolation. | FIXED |
| P9-B22 | High | `PackStepper` transaction publication | A later-batch Euler or exponential failure restored both arenas but leaked the rejected solver current/node/terminal solution, diagnostics, and cell/boundary heat vectors; 12 postcondition assertions failed across the two modes. | Preallocate and checkpoint all public solver/heat publication alongside arena bytes; restore it on every `Status` failure while invalidating only the numerical workspace cache. | Focused Debug/Release rollback gate 33/33; full PackStepper 172/172 in each build. | FIXED |
| P9-B23 | High | thermal derived arithmetic | Finite temperatures/conductances overflowed either an edge product or endpoint sum and returned `Success`; both registered sections were red. | Assemble finite differences, products, and transactional endpoint totals into preallocated trial buffers; publish heat/edge flux only after the whole graph validates. | Focused thermal/transaction gate plus full PackTopology 66/66 in each build. | FIXED |
| P9-B24 | Medium | PackStepper configure side effect | A topology rejected later for thermal-composition mismatch had already changed the caller batch's trusted lane period from 6 to 2. | Derive the period locally, allocate/configure every candidate first, then apply it immediately before statically no-throw member publication; explicitly handle the optimization's `Status`. | Failed-config regression preserves period 6; P2-G1 remains zero-allocation, 6/6 each build. | FIXED |
| P9-B25 | High | mutable/hostile thermal ownership | A corrupted incidence sign remained within the accepted enum range and returned `Success` with non-conservative heat. Separately, `PackStepper::solver()` exposed a mutable solver that callers could reconfigure away from the stepper topology; its compile-time contract was red. | Validate each incidence against its owning endpoint with preallocated low/high membership bits, and expose only `const PackSolver&` from the stepper. | Mutation-style incidence test rejects atomically; compile-time solver ownership contract and full Debug/Release topology/stepper suites pass. | FIXED |

| P9-B26 | High | Experiment parser resource amplification | Two individually legal `* 6000` steps committed 12,000 segments, while a padded 500-byte step repeated 10,000 times retained about 5 MB of duplicated source text; the pre-fix resource suite failed 6 assertions and preserved neither documented bound. | Cap expanded segments at 10,000, each step at 65,536 bytes, aggregate retained text at 4 MiB, and drive-cycle names at 1,024 bytes using overflow-safe pre-allocation arithmetic. Replace allocation-heavy regexes with the same explicit grammar and reject the empty `Run (A)` overlap before subtraction. | Debug/Release grammar 72/72 and full Experiment 199/199. Removing the minimum-name guard makes 3 registered assertions fail; late expansion allocation failure returns `Numerical_failure` without publication under persistent OOM. | FIXED |
| P9-B27 | High | BPX JSON/file parser trust | Leading-zero numbers, invalid raw UTF-8, vertical-tab/form-feed whitespace, a 4 MiB-plus source, and an ignored 65,537-value array were accepted or allowed unbounded tree/file construction; present-but-invalid optional fields were silently omitted. The registered pre-fix slices produced 11 failures. | Enforce RFC JSON number/whitespace/UTF-8 grammar, a 4 MiB wire cap, a 65,536-value tree cap, bounded exact file reads with EOF verification, and exact `Status` propagation for required, optional, derived, curve, activation, state, and default insertions. | Debug/Release ParameterSet 1,495/1,495, including valid raw UTF-8, finite-input derived overflow, sparse oversized file rejection, and unchanged output on every hostile input. | FIXED |
| P9-B28 | High | parser allocation/status atomicity | Optional BPX storage failures were collapsed to success; `ParameterSet::set()` allocated during canonicalisation outside its catch; parser catch diagnostics could allocate again during OOM; and BPX stream construction/open sat outside the file catch boundary. | Make `set()` and the complete file reader exception-translating transactions, use allocation-safe best-effort diagnostics, and target repeatable parser-owned allocations rather than arbitrary CRT/library ordinals. | Debug/Release allocation binary 906/906. Persistent failure covers a late Experiment reserve, direct map-node and canonical-name allocation, and all 176 BPX map-node-size occurrences. Restoring the old optional-status swallow fails at occurrence 149 with published partial output. | FIXED |
| P9-B29 | Medium | Experiment parser/runner semantic mismatch | The committed `Hold at 4.2 V until 3.8 V` fuzz seed exited 77 because parsing returned `Success`, while `CyclerV2::run()` rejects a voltage-controlled segment terminated by a voltage event. The ordinary pre-fix grammar test failed 7 assertions as the invalid success also replaced the prior output. | Run the shared internal `validSegment` semantic gate after syntactic parsing and before expansion reserve/publication, returning a step-local diagnostic on mismatch. | Full Debug and fast-math Release Experiment binaries pass 203/203. The new seed passes both optimized and unoptimized fuzz replays; post-fix Linux Experiment campaigns complete 149,399 Debug and 291,223 Release executions in 60 seconds with no sanitizer finding. | FIXED |

## Refuted candidates

| Candidate | Adversarial evidence | Result |
|-----------|----------------------|--------|
| Cycler chooses a later event when two roots fall in one accepted step. | Independent indicators at 0.25 s and 0.75 s with a 1 s trial select the named 0.25 s event. | REFUTED |
| Cycler misses an event exactly on a step boundary. | Indicator root exactly at 1 s with a 1 s step terminates at exactly 1 s. | REFUTED |
| Pack topology accepts empty/zero groups or duplicate thermal endpoints. | Compile-time validation paths reject the malformed descriptors before workspace mutation. | REFUTED (read-only audit; mechanised tests pending ledger expansion) |
| Ordinary finite thermal spans and one-cell/zero-thermal-edge packs are intrinsically invalid. | Compile validation and the scalar topology algebra support these cases; no failing invariant was found. | REFUTED (read-only audit) |
| Sparse external CSV node labels require allocation through the maximum label. | Labels up to `UINT32_MAX` are sorted, zero-wire representatives are contracted, and only the two dense representatives are allocated in the registered case. | REFUTED (mechanised Debug/Release) |
| A zero-ohm CSV resistor can simply be dropped. | Dropping the row disconnects or changes the graph; union-find contraction preserves the ideal-wire topology and admits the sparse-label regression. | REFUTED (mutation-sensitive construction) |

## Import compatibility limits

- `NetlistCsv` imports liionpack's DataFrame CSV schema, not its LTSpice
  `.cir`/`.txt` reader format.
- Every positive `R*` row, including `Ri*`, remains a literal resistor in
  addition to the SLIDE cell model's Thevenin resistance. This is not a
  behavioural-parity claim.
- V/I magnitudes, original resistor descriptors, and external node labels are
  not retained in `CompiledPackTopology`, so only electrical equivalence—not a
  lossless source round trip—is currently representable.
- The fixed 4 MiB/100,000-row trust budget intentionally excludes ordinary
  exports of the largest 100,000-cell packs. A future validated limits policy
  is required before claiming that import scale.

## P9-G2 fuzz-gate evidence

- Three raw-byte libFuzzer drivers cover Experiment, BPX, and NetlistCsv. Each
  parses twice for deterministic status/diagnostics/output, starts from a
  fully populated poison value, requires exact poison preservation on every
  rejection, and requires poison removal plus structural validity on success.
- Netlist success is not delegated to the production validator. The harness
  independently recomputes endpoint ranges, cell-ID bijection, path/location
  uniqueness, resistor finiteness, BFS connectivity, exact sorted sparsity,
  ladder classification/orientation, and the empty imported thermal graph.
- Oracle mutation checks are red on the first valid seed: append instead of
  replace in Experiment exits 77; merging BPX into the old output exits 77;
  publishing stray netlist thermal scratch exits 77. All mutations were
  reverted and the clean targets rebuilt. A second independent review found
  that `std::isfinite` could be folded away by Release fast-math and that
  callback poison identity was under-specified. The oracle now classifies IEEE
  exponent bits through an opaque reference/volatile load and uses named
  callback targets with exact markers. Publishing a bit-constructed quiet NaN
  makes the final Release target exit 77; the reverted clean seed replay passes.
  Review also found and fixed the
  Experiment byte adapter treating a normal trailing newline as an empty step.
  The added voltage-control/voltage-event seed then exposed P9-B29 before gate
  closure; all Experiment campaign counts in this report are post-fix.
- The instrumented core is a distinct `slide_core_fuzz` target. A combined
  Windows configuration built both the ordinary Debug core smoke and the
  ASan+UBSan fuzzer, proving sanitizer/static-CRT settings do not cross the ABI
  boundary. Ordinary full-tree Debug and Release `slide_core` builds and the
  core-only 1/1 smokes remain green.
- Linux/Clang 18 ASan+UBSan+LSan, exact 60-second campaigns with committed
  dictionaries and 65,536-byte mutation bounds:

  | Build | Experiment | BPX | NetlistCsv | Slowest input | Peak RSS |
  |-------|-----------:|----:|-----------:|--------------:|---------:|
  | Debug | 149,399 | 73,981 | 94,933 | 0 s | 526 MiB max |
  | fast-math Release | 291,223 | 546,881 | 556,946 | 0 s | 576 MiB max |

- Separate Linux Debug/Release seed replays cover generated 65,537-byte
  Experiment steps and 4,194,305-byte BPX/netlist documents. Windows Clang 21
  Debug/Release ASan+UBSan campaigns and the same size replays also pass; only
  Windows disables `detect_container_overflow` because the prebuilt libFuzzer
  runtime and MSVC STL container annotations use incompatible ABIs. Heap
  redzones, UBSan, and strict string checks remain enabled there.
- `.github/workflows/core-fuzz.yml` reproduces the Linux Debug/Release matrix,
  full Linux annotations/leak detection, limit replays, three 60-second
  campaigns, fixed seeds, writable corpus copies, and failure artifacts. It
  is committed and syntax-checked, but the hosted jobs cannot run before a
  push; no hosted-green claim is made.

## Validation protocol

- All scenarios are parser-only or at most one accepted model step; no long
  trajectory was used.
- Debug and Release are separate binaries. Release uses the repository's
  fast-math flags, and NaN/Inf tests construct IEEE bit patterns directly.
- P9-G2 is closed by the evidence above. Full-suite sanitizer/TSan evidence,
  measured Status-branch coverage, and remaining subsystem audits are still
  required before P9-G1, P9-G3, and P9-G4 can close.
