# SLIDE Development Discussions

This document tracks design decisions, architecture evolution, and session notes for long-term pair programming collaboration on SLIDE.

---

## Active Design Decisions

### [2026-01-29] Initial Codebase Review & Collaboration Setup

**Status**: Decided

**Context**: First comprehensive review of SLIDE codebase to establish pair programming infrastructure for multi-year development.

**Key Findings**:
1. Architecture is well-structured: StorageUnit → Cell → Module → Battery hierarchy
2. C++20 codebase with modern practices (Deep_ptr, std::span, constexpr)
3. Good test infrastructure (Catch2) but some tests commented out
4. Python/MATLAB bindings planned but not implemented
5. TODO.md has 365 lines of unstructured items

**Decision**: Create comprehensive documentation and tracking infrastructure:
- Updated CLAUDE.md with SLIDE-specific domain knowledge
- Updated style guides (cpp-style.md, python-style.md) for SLIDE
- Create skills for common operations
- Reorganize TODO.md with priorities
- Use this discussions.md for ongoing decision tracking

---

### [Pending] Python Bindings Strategy

**Status**: Open

**Context**: Python bindings are a high-priority feature request for broader adoption.

**Options**:
1. **pybind11** - Mature, well-documented, widely used
   - Pros: Large ecosystem, good numpy integration
   - Cons: Compile times, header-only can increase binary size

2. **nanobind** - Modern successor to pybind11 by same author
   - Pros: Faster compile, smaller binaries, better C++17/20 support
   - Cons: Smaller community, fewer examples

3. **SWIG** - Language-agnostic wrapper generator
   - Pros: Can generate Python and MATLAB simultaneously
   - Cons: Less Pythonic output, steeper learning curve

**Decision**: TBD - Need to evaluate compile time impact and numpy integration quality

---

### [Pending] `redistributeCurrent_new` Performance

**Status**: Open

**Context**: The current redistribution algorithm for parallel modules can require up to 2500 iterations (TODO.md line 79), causing significant slowdowns.

**Options**:
1. **Newton-Raphson solver** - Quadratic convergence
2. **Direct linear algebra** - Solve coupled equations directly
3. **Adaptive tolerance** - Accept less precision when converging slowly
4. **Better initial guess** - Use previous solution as starting point

**Decision**: TBD - Need profiling data to identify root cause

---

## Architecture Evolution Log

### Version 3.x Roadmap

**Current**: v3.0.0 (slide-pack merged, C++20)

#### Planned Milestones
- [ ] **v3.1.0**: Python bindings (pybind11/nanobind)
- [ ] **v3.2.0**: MATLAB MEX interface
- [ ] **v3.3.0**: PyBaMM-compatible Experiment interface
- [ ] **v3.4.0**: SUNDIALS solver integration (optional)
- [ ] **v4.0.0**: GPU acceleration (CUDA, optional)

#### Completed Milestones
- [x] v3.0.0 (Current)
  - Merged SLIDE and slide-pack
  - Upgraded to C++20
  - Added Cell_ECM, Cell_Bucket
  - Converted shared_ptr → unique_ptr/Deep_ptr
  - Added std::span for state access
  - Status class for error handling
  - Model_SPM optimization: 34s → 12s (2.8x speedup)
  - File I/O elimination: 30s improvement

---

## Performance Benchmarks History

| Date | Version | Test | Time | Hardware | Notes |
|------|---------|------|------|----------|-------|
| 2026-01-29 | v3.0.0 | 5000 1C CC cycles | <1 min | - | Baseline from README |
| 2026-01-29 | v3.0.0 | CC+CV cycles | <2 min | - | Baseline from README |
| 2026-01-29 | v3.0.0 | EPFL battery 1hr CC | ~2 sec | - | From TODO.md |
| 2026-01-29 | v3.0.0 | EPFL battery + aging | ~3.5 sec | - | From TODO.md |

---

## Known Technical Debt

| Priority | Issue | File | Notes |
|----------|-------|------|-------|
| High | `redistributeCurrent_new` 2500+ iterations | Module_p_impl.cpp | Performance bottleneck |
| High | MSVC 3x slower than Clang | CMakeLists.txt | Vectorization flags |
| Medium | Raw `assert()` in tests | Cycler_test.cpp | Should use Catch2 REQUIRE() |
| Medium | Mixed `#define` and `constexpr` | settings.hpp:101 | DATASTORE_BATT macro |
| Medium | `StorageUnit::copy()` returns raw ptr | StorageUnit.hpp | Should return unique_ptr |
| Medium | C-array `double Tneighb[]` params | Module.hpp:284 | Use std::span |
| Low | Inconsistent Status codes | Various | Some int, some Status enum |
| Low | Global `settings::isParallel` | settings.hpp | Encapsulate thread control |
| Low | `therm.Qcontact` 6x overestimated | Module thermal | Calculation bug |

---

## Session Notes

### [2026-01-29] Initial Codebase Exploration

- **Focus**: Complete exploration of SLIDE library, documentation setup, collaboration infrastructure
- **Decisions**:
  - Style guides updated in-place (not new files) for SLIDE
  - Completed TODO items archived to separate COMPLETED.md
  - Session notes use summary bullets + rationale format
- **Next Steps**:
  - Create skill files for common operations
  - Restructure TODO.md into prioritized format
  - Run test suite to verify build still works
- **Blockers**: None

---

## 2026-07-07 — v4 architecture decided (Fable session)

- **Decision**: Full v4 refactor architecture fixed and recorded in [/PLAN.md](../PLAN.md) §3–§4 (14-entry
  decision log with rejected alternatives). Headlines: SoA state arena with archetype batches (D-01),
  compile-time model composition with per-batch dispatch (D-02), netlist compiled flat to sparse MNA with
  three solver modes — sparse Newton / Thomas ladder / waveform-relaxation+Baumgarte (D-03..05), exponential
  modal propagator replacing forward Euler on diffusion (D-07), Strang multirate (D-08), strangler migration
  with digit-diff parity gates (D-14).
- **Evidence**: three-scout audit 2026-07-07; findings condensed in PLAN.md §2 with file:line citations.
  Chebyshev nch≠5 root-caused (centre-node `-0.5` sign + Eigen inverse aliasing; already fixed on `Claude`).
- **Next Steps**: Phase 0 bug fixes (2 Opus agents, in flight); Phase 1 core data model after Volkan reviews
  PLAN.md §3/§4 and §7 open questions.
- **Blockers**: PLAN.md §7 Q1–Q7 await Volkan's answers (assumptions recorded, work proceeds on them).

---

## 2026-07-08 — Phase 1: Q8 Release re-confirm + production diffusion kernel

- **Q8 fully closed.** P1-G0 pilot re-run in Release/-O3 (`build-release`): decisive band rel ≤ 1e-12
  HOLDS (max_rel 5.63e-15); exact bit-identity (Debug: 0) FALSIFIED under -O3 cross-TU FMA contraction
  (max_abs 2.26e-17). `-ffp-contract=off` rejected — can't fix cross-TU on the test target alone, and
  forcing it globally recompiles legacy (§5.1 forbids). Bonus hypothesis falsified & recorded; decisive
  band untouched. Pilot's exact-zero CHECK scoped to Debug. See `handoff-2026-07-08-phase1-diffusion.md`.
- **Production `SpectralDiffusion<NCH>` landed** (`src/core/SpectralDiffusion.hpp`): vectorised SoA
  forward-Euler diffusion. Validated vs the legacy-shaped oracle — Debug bit-identical (proves math
  identity), Release rel 3.8e-15 ≤ 1e-12. Two-statement inner form; once-allocated scratch (PC-1).
- **Design decision surfaced (OPEN for Volkan/Fable):** ThermalLumped + ageing kernels are NOT
  self-contained — they need the observable-reconstruction layer (c_surf = C·z + D·flux, overpotential,
  dOCV, Rdc; §3.7/D-10) and the §3.11 BatchView/StepCtx interface, both DEFERRED this session (diffusion
  took plain spans as a stopgap). Recommend designing that layer + doing **P1-G3 (Chebyshev/Carslaw &
  Jaeger oracle)** — which exercises the same C/D output path — before porting thermal/ageing.
- **Review channels down:** advisor unavailable, Fable agent out of usage credits — no external review
  this session. Decisions rest on the Debug bit-identity arbiter + derived accumulation bounds.

## 2026-07-10 — P8-G0 optionality and portability

- **Finding:** the existing `SLIDE_CORE_ONLY` configuration built successfully but registered zero tests, so its apparent CTest success could not support P8-G0. An independent review also found that `${CMAKE_SOURCE_DIR}/src` broke public includes when SLIDE was embedded and that CUDA's implementation macro leaked publicly.
- **Decision:** add an opt-in dependency-light smoke plus a nested-superproject external consumer, use `PROJECT_SOURCE_DIR`, keep CUDA metadata private, prefer installed Eigen before the pinned fallback, and add a Linux/macOS/Windows core-only CI matrix. D-25 defines the registered “dependency-free” phrase as optional-toolchain-free; Eigen remains the required cold-path linear-algebra dependency.
- **Evidence:** root core-only smoke 1/1, nested consumer 1/1, optional-off Debug 49/49, optional-off Release 49/49, and rebuilt installed CPython 3.13 wheel 9 passed/2 expected optional skips on Windows. Existing installed-wheel CI remains Linux/macOS/Windows × Python 3.10/3.13; `cmake/**` now triggers it. Cross-platform jobs require a push before GitHub can execute them.
- **Next:** P8-G5 tested v4 documentation.

## 2026-07-10 — P8-G5 executable v4 documentation

- **Finding:** the published tree was v3-only, the old docs workflow did not execute examples, and a nominal local Jekyll success could still emit unthemed pages. Adversarial validation also found a missing Liquid-tag plugin, a blank Doxygen main page, a broken custom Doxygen header, shell-dependent wheel globs, and several over-broad compatibility claims.
- **Decision:** make `docs/v4` dominant while retaining explicit v3 banners; treat fenced quickstarts as source code through one standard-library extractor; check local links, heading fragments, and required Jekyll layout metadata; pin the remote theme and scope Pages/OIDC permissions to deployment. Keep MATLAB on the licensed machine rather than claiming an unverified hosted MEX lane.
- **Evidence:** the fresh C++ external consumer passes 1/1; installed-wheel Python and MATLAB R2025b each return 7 finite samples and final voltage 3.879196 V. Doxygen 1.14 returns zero generator errors and renders the v4 main page (123 retained v3 source-comment warnings); production Jekyll emits 8 themed v4 pages with `/SLIDE` URLs and `master/docs` edit links. Ruff, yamllint, actionlint, link/front-matter validation, and `git diff --check` pass.
- **Review:** an independent read-only adversarial pass checked documentation claims against CMake, C++, Python, and MATLAB sources and reported no remaining P8-G5 blocker after fixes.
- **Next:** Phase 9A audit-debt closure (AUD-1, AUD-2, AUD-4; AUD-3 remains Q11).

## 2026-07-10 — Phase 9A audit-debt closure

- **Finding:** P2-G1 carried unexplained `1e-10` V/I assertions despite Phase 9A's `1e-12` requirement. The tightened Debug/Release gate passes, but worst current drift uses 36.95% of the band, proving the correction is material. An independent read-only audit also confirmed that `09e8ec5` introduced the Phase-5 implementation and every numeric band together.
- **Decision:** retain the Phase-5 bands as useful empirical regression sentinels while explicitly labelling them post-hoc. The `0.2 µAh` scale is consistent with `20 µA × 30 s = 0.1667 µAh`; the 45-iteration event bisection leaves at most a `2.84e-14 s` bracket for the tested one-second step. Neither derivation turns the values into universal accuracy guarantees.
- **Decision:** D-26 waives mandatory CVODE only for the exact frozen diagonal modal subflow. P3-G1 now evaluates the closed form through long-double `exp` and `expm1(rate·h)/rate`, with an exact zero-rate limit and no copied production Taylor branch. The unrun full 1C/current-step voltage comparison remains explicitly unclaimed.
- **Decision:** D-27 resolves Q11 by waiving the historical quiet-host/operator condition for v4.0, while retaining PAY-1/2/4 as qualified development-host evidence. The claim that quiet conditions could only improve a ratio was removed. PAY-4 raw JSON is now committed; no weak rerun on the same busy host was performed.
- **Evidence:** targeted Debug and Release CTest each pass 3/3 (`P2G1_pack`, `core_Experiment`, `core_ExponentialModal`). Release Phase-5 reproduction observes 0.101674 µV, 12.479451 µA, and 0.116939 µAh against the unchanged empirical envelopes. Full derivation: `.claude/reports/p9a-audit-debt-2026-07-10.md`.
- **Next:** Phase 9B systematic adversarial bug-hunt.

## 2026-07-10 — Phase 9B Experiment adversarial pass

- **Confirmed:** direct `ExperimentSegment` values bypassed parser invariants; invalid enums and NaN/Inf metadata could enter execution. Duration parsing checked finiteness before, but not after, unit scaling.
- **Confirmed:** `CyclerV2` checkpointed before an event-aware trial but did not restore on several post-advance failures. A throwing custom event callback therefore returned an error while leaking the advanced arena.
- **Fix:** one exhaustive cold-path validator runs before output/state initialization; invalid integrator values are rejected at configure time; every fallible operation between trial advance and commit restores the byte-exact arena. The event-corner tests also refuted incorrect-first-root and exact-breakpoint concerns.
- **Evidence:** the registered tests were written first and the pre-fix `[P9]` run produced 13 assertion failures. The corrected Debug and `-Ofast` Release binaries each pass the complete Experiment subsystem: 183 assertions in 9 cases. The focused slice is 42/42 P9 assertions plus 56/56 parser-atomicity assertions.
- **Next:** continue the PackSolver/topology, PackStepper, recorder, curve, parser-fuzz, sanitizer, and Status-coverage passes; do not close 9B until P9-G1..G4 all pass.

## 2026-07-10 — Phase 9B synchronous Recorder pass

- **Confirmed:** once the fixed-capacity recorder entered `thin`, ordering was checked only against the last stored slot. Duplicate and backward omitted steps therefore succeeded and inflated the thinned count.
- **Confirmed:** `record()` validated total current but copied a non-finite elapsed time and could derive infinite current density from individually finite area/current. An invalid `BackpressurePolicy` also configured successfully by falling through as `stop`.
- **Fix:** mirror the async monotonic watermark, reset it on `clear()`, validate time and current-density quotients before advancing it or touching storage, and exhaustively validate the policy enum during atomic configuration.
- **Evidence:** the initial registered P9 run failed 7 assertions; after those fixes, the added invalid-enum test failed 2 assertions while the earlier cases passed. Final focused Debug/Release gates pass 20/20, and sequential full Recorder binaries pass 98/98 in each build. Parallel full-binary execution was deliberately rejected as evidence because both configurations share fixed temporary filenames and collided at the filesystem layer.

## 2026-07-10 — Phase 9B AsyncRecorder format/input pass

- **Confirmed:** a valid-CRC 64-byte header with `snapshots=UINT64_MAX` reached `vector` construction and let `length_error` escape the `Status` API. A CRC-valid raw codec value `256` narrowed to byte-sized `none` and decoded successfully.
- **Confirmed:** invalid async backpressure values configured a worker, and finite `1e300 A / 1e-300 m²` was published as infinite density and written as a snapshot.
- **Fix:** prove the file can hold at least one 64-byte block header per claimed snapshot before allocation, catch `length_error`, reject raw codec values before narrowing, exhaustively validate the policy, and validate the density quotient before taking the ring lock/publishing a slot.
- **Evidence:** the test-first Debug run failed 8 assertions, including the escaped exception. Focused Debug/Release passes 67/67; sequential full async binaries pass 233/233 each. The suite uses at most eight model steps and a bounded 50-enqueue pressure case.
- **Deferred deliberately:** configuration path-copy fault atomicity and public shuffle preconditions remain open until their own red fault/UB tests are registered; this commit does not claim them.

## 2026-07-10 — Phase 9B async configure fault atomicity

- **Confirmed:** deterministic failure of the late `path_ = path` allocation returned `Numerical_failure` after `batch_` had already been published. `configured()` therefore remained true without a worker/output member, and a zstd build would also strand the still-local context.
- **Test design:** warm one-time library state, measure the allocation size of copying the deliberately long path, and fail only that matching allocation. This avoids injecting into unrelated standard-library internals and produced one precise pre-fix assertion failure without timeouts.
- **Fix/evidence:** copy the path locally before slot/file/context creation and swap it during no-throw member commit. The allocation executable passes 15/15 in Debug and Release, including successful reuse after failure and the original zero-allocation enqueue check.
- **Still open:** public byte-shuffle preconditions, close-time errors, and concurrent finish/hook contracts remain separate ledger candidates.

## 2026-07-10 — Phase 9B PackSolver/Mode-C adversarial pass

- **Confirmed:** Mode C accepted tiny-gain stagnation from current delta alone while KCL remained large. It measured only terminal KCL, so an internal residual could be hidden. The historical general `(1-α)^k` diagnostic was not valid for coupled networks.
- **Confirmed:** invalid solve/node/branch metadata reached plausible fallthroughs; malformed compiled indices, connectivity, sparsity, lane mappings, and ladder tables were trusted. Sparse/ladder/relaxation paths could derive non-finite values without a bit-safe publication guard, and diagnostics leaked across modes.
- **Independent review blockers fixed before commit:** remove mutable `workspace()` ownership (compile-time red contract), reject graph dimensions before proportional allocation/signed indexing, clear stale terminal voltage on successful reconfigure, replace the non-discriminating one-cell overflow test, and add an internal-node KCL construction whose terminal equation is satisfied.
- **Numerical correction:** strict-FP compensated accumulation reduces PAY-3 100k KCL drift from `4.986991736e-6` to `2.54658516e-9` Debug and `2.61934474e-9` Release. The first Kahan attempt failed Release because `-Ofast` reassociated it; compiler-scoped strict-FP controls fixed that rather than loosening the `1e-8` gate.
- **Evidence:** pre-fix invalid-mode/kind, tiny-gain, diagnostic, and publication checks failed; the two-series overflow regression is mutation-red when the new ladder/publication guards are disabled. Final Debug and Release binaries pass PackTopology 50/50, PackSolver 667/667, and ModeC 52/52; the 100k gates take 7.9/8.7 s Debug and 2.1/2.7 s Release.

## 2026-07-10 — Phase 9B byte-shuffle public contract

- **Confirmed red:** Release terminated on `byteShuffle({}, {}, 0)`; the only guards were Debug assertions, and undersized/overlapping output could write out of bounds.
- **Fix:** both transforms return `Status`, validate width/shape/divisibility/non-overlap, and have writer/reader callers propagate unexpected failures. The direct involution test now checks successful statuses.
- **Evidence:** full AsyncRecorder Debug and Release binaries pass 238/238; the degenerate test is eight bytes and performs no model step.

## 2026-07-10 — Phase 9B compiled-curve index safety

- **Confirmed red:** the no-simulation Debug regression failed 6/21 assertions. A denormal domain published an infinite accelerator reciprocal, extreme finite ordinates admitted a non-finite segment slope, and NaN derivative lookup reached a float-to-integer conversion.
- **Optimized-build finding:** the first guarded implementation passed Debug but access-violated in Release. ThinLTO had turned the finite guard into `llvm.assume` because its rejected branch visibly returned a constexpr NaN under `-ffinite-math-only`.
- **Fix:** validate every derived build/sample/index quantity before conversion or commit, use a reference-based IEEE classifier for public floating inputs, and construct the invalid sentinel through an opaque no-inline integer path. Keep the ordinary interpolation expression unchanged.
- **Evidence:** full CompiledCurve Debug and fast-math Release binaries each pass 163/163 assertions, including raw sentinel bits, invalid tolerance, extreme/tiny inputs, BPX endpoint resolution, and bit-exact legacy OCV interpolation. No model step or trajectory was run.

## 2026-07-10 — Phase 9B PackStepper and thermal transaction pass

- **Confirmed red:** both thermal overflow sections returned `Success`; a rejected configure changed the caller's trusted period 6→2; later-batch Euler/exponential failures restored arenas but leaked solver publication, diagnostics, and heat (12 failed postconditions). Review then found a corrupted incidence sign returning non-conservative `Success` and a compile-time mutable-solver ownership hole.
- **Fix:** checkpoint the full public publication surface in preallocated storage, restore it while invalidating only workspace cache validity, assemble heat into trial buffers with endpoint/incidence proof, and defer caller metadata until candidate construction is complete. `PackStepper::solver()` is const-only.
- **Evidence:** Debug and fast-math Release each pass PackTopology 66/66, PackStepper 172/172, and the P2-G1 zero-allocation gate 6/6. The rejection test takes one attempted step in each integrator; no trajectory or performance simulation was run.

## 2026-07-10 — Phase 9B Experiment/BPX parser hardening

- **Confirmed red:** aggregate Experiment repetition committed 12,000 segments and padded repeats amplified retained source text; BPX accepted non-JSON number/whitespace/UTF-8 forms, unbounded ignored values, and invalid present optionals. Pre-fix focused slices produced 6 Experiment and 11 BPX failures.
- **Architecture:** replace the two Experiment regex matches with explicit bounded grammar, cap semantic and retained expansion independently, and parse BPX behind both wire-byte and tree-value budgets. File absorption now proves exact bounded reads and EOF before parsing.
- **Allocation contract:** `ParameterSet::set`, both parser transactions, and the complete BPX file-open/read path translate allocation and length failures without publication. Error diagnostics are best-effort and cannot throw a second allocation failure through the `Status` API.
- **Adversarial evidence:** persistent size-targeted faults cover a late Experiment reserve, direct map-node and canonical-name allocation, and every BPX allocation matching the measured `ParameterSet` node size. Reintroducing the optional-status swallow fails at occurrence 149/176; removing the empty drive-name guard fails 3 assertions.
- **Validation:** Debug and fast-math Release each pass Experiment 199, ParameterSet 1,495, parser-allocation 906, ThermalLumped 25, CompiledCurve 163, PackSolver 667, and PackStepper 172 assertions. No long trajectory was run; Experiment simulations are bounded unit scenarios.
- **Next:** implement and fuzz the missing liionpack-compatible netlist CSV parser, then complete the Experiment/BPX/netlist P9-G2 corpus and CI driver.

## 2026-07-10 — Phase 9B liionpack CSV importer

- **Schema decision:** follow upstream's `desc,node1,node2,value` DataFrame CSV and accept extra generated coordinate columns. `V*` is a cell from node1 (+) to node2 (−); positive `R*` remains a resistor; zero `R*` contracts endpoints; one excluded `I*` row defines the positive/negative terminals.
- **Architecture:** arbitrary graphs compile directly into a local `CompiledPackTopology`; a shared internal finalizer now owns metadata plus the exact validator used by `SolverWorkspace`. Sparse uint32 labels are sorted and densely remapped only after deterministic union-find contraction.
- **Evidence:** the initial stub made 3 valid/file assertions red. Final Debug and fast-math Release each pass 728 importer assertions, PackTopology 66, PackSolver 667, and the expanded persistent parser-allocation gate 913. The pure-cell 2s2p CSV receives the same ladder offsets/cell order as the combinator fast path.
- **Claim limits:** literal `Ri*` can add resistance beyond the model-owned cell resistance; V/I magnitudes are ignored after validation; descriptors/node labels are not retained for lossless export; and the 4 MiB budget does not admit ordinary 100,000-cell tables. This is graph-schema compatibility, not liionpack waveform parity.
- **Next:** add three libFuzzer drivers, committed corpora/dictionaries, and a bounded Clang ASan+UBSan campaign; P9-G2 remains open until that evidence is green.

## Quick Links

- [CLAUDE.md](CLAUDE.md) - Main runbook
- [cpp-style.md](cpp-style.md) - C++ conventions
- [python-style.md](python-style.md) - Python conventions
- [develop/TODO.md](../develop/TODO.md) - Active development items
- [develop/COMPLETED.md](../develop/COMPLETED.md) - Archived completed items
- [CHANGELOG.md](../CHANGELOG.md) - Release history
