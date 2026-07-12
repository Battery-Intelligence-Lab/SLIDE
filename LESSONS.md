# SLIDE Engineering Lessons

This file records reusable findings from the v4 implementation and validation work. Detailed design decisions remain in `PLAN.md`; session chronology remains in `.claude/discussions.md` and `.claude/summaries/`.

## Single-source numerical kernels still need code-generation provenance

- An algebraically identical inline helper can change fast-math loop IR, vector-loop/code-generation shape, and final bits. Freeze whole traces before refactoring, inspect optimized code, and use a narrowly documented shared statement/expression kernel when an ordinary function boundary breaks a required operation-order contract.
- Recorded bits are compiler, CRT, ISA, and build-configuration evidence—not portable mathematical truth. Make them explicit opt-in, retain the historical producer commit, and never bless new post-refactor hashes merely because another toolchain differs; portable analytic, finite-difference, parity, and structural gates remain mandatory.
- A single-source gate must enumerate every production consumer and keep independent oracles outside that source. The first CPU/CUDA/Dual map missed `SpectralDiffusion`; expanding the structural scan found it, while forbidding the raw legacy oracle from including the shared header preserves its ability to catch common-mode wiring errors.
- Preflight multiplication without evaluating the dangerous product in the guard. Once finite nonnegative factors are established, a factor at most one cannot cause overflow; for larger factors compare safely against the maximum, treat a rounded equality conservatively, and validate every derived parameter scale before either compute or RHS publication.
- Differentiated observations need the full primal algebra, including temperature corrections that vanish at a reference fixture. A reference-temperature exact trace cannot expose an omitted entropic tangent, so pair it with an off-reference difference-of-differences arbiter.

## Portability gates must exercise a consumer

- A successful `SLIDE_CORE_ONLY` build followed by a zero-test CTest run is not validation. Keep a dependency-light executable smoke that checks the state arena, owned thread pool, and disabled optional capabilities.
- An in-tree target cannot prove that public CMake usage works. Configure SLIDE as a nested subproject and compile headers through `target_link_libraries(slide_core)`; this exposed the incorrect use of `CMAKE_SOURCE_DIR`, which changes meaning in a superproject.
- “Dependency-free” must be scoped precisely. The v4 core is free of optional CUDA, MATLAB, zstd, Arrow, TBB, and legacy dependencies, but Eigen remains required for cold-path spectral and sparse solves. Prefer an installed package and keep a pinned source fallback.
- CI path filters are part of correctness: changes under `cmake/**` must trigger both core portability and installed-wheel jobs.

## Documentation is an executable interface

- Keep runnable examples in one place: mark the fenced Markdown block, extract that exact text, and compile or run it away from the source tree. A copied example file can pass while the published page has already drifted.
- Validate rendered structure as well as tool exit codes. Jekyll returned success for pages without `layout: default` but emitted unthemed HTML; Doxygen returned success while its configured main page was outside `INPUT` and therefore blank. Gates now require the Jekyll front matter and inspect both rendered landing pages.
- A documentation theme is a build dependency. Pin its immutable revision and load every Liquid-tag provider explicitly; the prior remote theme referenced `github_edit_link` without enabling `jekyll-github-metadata`.
- Cross-platform commands must avoid shell-expanded wheel globs. Select the local artifact with pip's `--no-index --find-links=dist` interface instead.
- Do not turn a licensed local MATLAB result into an unverified hosted-runner claim. The same extractor supports the gate, but MEX execution remains on an explicitly provisioned licensed machine.
- Doxygen's 123 current warnings all originate in retained v3 `src/` comments. The v4 docs gate removes generator errors and new-page/configuration warnings without pretending that this legacy debt is already zero.

## Numerical gates need provenance as well as tight numbers

- Never silently loosen a registered numerical band in the assertion. P2-G1's `1e-10` voltage/current checks looked harmless beside tiny observed drift, but the corrected `1e-12` current gate has only 2.71× Debug headroom and can catch a materially smaller regression.
- A threshold introduced in the same commit as its implementation is a post-hoc empirical regression sentinel, not a preregistered accuracy oracle. Keep useful sentinels, but label their provenance and derive only what the contracts support.
- Dimensional consistency is not logical implication. The Phase-5 `0.2 µAh` envelope is consistent with `20 µA × 30 s = 0.1667 µAh`, yet a final-current check cannot bound the whole integrated trajectory; the charge assertion remains independent.
- Prefer an exact analytic oracle over a tolerance-converged numerical reference when the isolated subproblem is genuinely linear and closed form. Independence still matters: evaluate through a different branch and operation order, and state what shared coefficients define the ODE rather than pretending the oracle validates them.
- Scope every waiver narrowly. D-26 replaces CVODE only for the frozen diagonal modal subflow; it does not validate a generic RHS adapter, nonlinear integration, or splitting. D-27 waives a historical timing protocol, not the need to label development-host timings as qualified evidence.
- Busy-host noise has no guaranteed direction in a ratio. Preserve raw timing ranges, repeat both tools, and require a named stable host before making an unqualified performance claim.

## A transaction boundary includes every fallible post-step observation

- Saving state before an integrator call is insufficient if voltage reconstruction, algebraic control, or an event callback can fail after the state changes. Restore on every path until the step is both observed and committed.
- Validate public value-semantic descriptors before constructing output or touching an arena. Exhaustive enum checks are required because an invalid scoped-enum value otherwise reaches a plausible default branch in Release.
- Revalidate derived quantities after unit conversion and arithmetic. Finite input is not evidence that `value * scale` is finite, especially under the project's fast-math build; adversarial tests should inject IEEE values by bit pattern rather than relying on optimizer-sensitive `infinity()` construction.
- A dropped item is still part of an ordered stream. Thin/backpressure paths must advance the same monotonic watermark as stored items, or duplicate and out-of-order retries become silently accepted.
- Validate serialized discriminants at their wire width before narrowing to an in-memory enum. A CRC proves bytes were preserved; it does not prevent `256` from aliasing byte-sized value `0` after a premature cast.
- Bound allocation counts by the bytes that must physically exist in the input before constructing vectors. Arithmetic-overflow checks alone still permit `vector(max_size)` and an escaping `length_error` when a tiny file claims an impossible count.
- Exception atomicity is easiest when every throwing copy/allocation happens in locals before the first public member is committed. This includes already-constructed movable runtimes such as executors: publishing one early still creates a hybrid object if later scratch allocation fails. Move/swap only through a statically checked no-throw commit phase; deterministic allocation faults should prove both continued use of the prior configuration and successful later reconfiguration.

## Iteration convergence requires the equation residual

- A small update can mean convergence or merely a tiny relaxation gain. Mode C must gate on both current-update size and the KCL residual at every independent node; terminal KCL alone can be exactly zero while an internal node violates conservation.
- Do not promote a one-unknown contraction identity into a theorem for a coupled graph. Keep the exact `(1-α)^k` oracle on the topology where it is derivable, and measure the governing residual directly elsewhere.
- Fast-math may reassociate even a textbook compensated sum. A compensation kernel needs a scoped strict-FP implementation plus Debug/Release evidence; its roundoff diagnostic should be labelled an estimate, not a portable theorem.
- A read-only diagnostics accessor must not also expose movable/reconfigurable solver ownership. Const public surfaces are a correctness boundary, not only an API-style preference.
- A regression must distinguish the old bug. The first extreme-value test already failed on the old solver and was rejected during review; the replacement keeps current deltas finite while terminal accumulation alone overflows, and mutation testing makes it red.
- An assert is not a public error contract. Low-level transforms reachable by callers need validated sizes, widths, and aliasing plus a propagated `Status`; otherwise Release turns a simple bad argument into divide-by-zero or out-of-bounds UB.

## Fast-math safety includes control-flow provenance

- Validate every floating intermediate before converting it to an integer. Finite knots do not imply a finite range, slope, reciprocal bin width, or scaled lookup index; extreme and denormal inputs independently break those implications.
- A bit-safe finite predicate is necessary but not sufficient under `-ffinite-math-only`. If the guarded branch visibly returns a constexpr NaN, ThinLTO can declare that branch unreachable and replace the check with an assumption. Keep invalid sentinels opaque (or compile the boundary with strict FP) and inspect raw IEEE bits in the optimized regression.
- Preserve the normal finite arithmetic path while hardening exceptional values. The CompiledCurve fix retains the existing operation order, so the legacy non-knot OCV oracle remains bit-exact rather than being weakened to a tolerance after the change.

## Transaction rollback includes publication, not only physics rows

- Restoring every arena byte is insufficient if a failed multi-batch step leaves the rejected solver solution, diagnostics, or heat outputs visible. Checkpoint the complete public observation surface; invalidate derived factorization caches rather than pretending their execution history rolled back.
- Caller-owned optimization metadata belongs outside the fallible configure phase. Allocate and validate candidates first, apply the lane-period optimization last, and make the following member publication mechanically no-throw.
- A graph's sizes and enum ranges do not prove its incidence semantics. Validate that each thermal edge appears exactly once at its low endpoint with `+1` and once at its high endpoint with `-1`; preallocated membership bits preserve the zero-allocation hot path.
- Do not return a mutable subsolver from a composed transaction owner. Reconfiguration through that reference can invalidate topology and checkpoint invariants even if each class is locally valid.

## Untrusted parsers need grammar, resource, and failure budgets

- A per-token repetition limit is not an aggregate allocation bound. Prove expanded element count and retained byte count before `reserve()` or copying; include duplicated source and metadata, not only the semantic payload.
- JSON library-like behaviour must still enforce the wire grammar: only four whitespace bytes are legal, numbers cannot have leading zeroes, and raw UTF-8 must reject overlong encodings, surrogates, invalid continuations, and code points beyond U+10FFFF.
- Optional means absent is allowed, not present-but-invalid is ignored. Preserve `Invalid_parameters` versus `Numerical_failure` through every required, optional, derived, and default insertion.
- Allocation-fault tests should target repeatable allocations owned by the transaction. Arbitrary global ordinals can enter CRT or standard-library internals and abort instead of testing the intended boundary. Persistent failure after the selected allocation also tests whether the error diagnostic itself can fail safely.
- Under fast-math, keep finite classifiers behind a reference boundary as well as an integer bit barrier. Individually finite parsed operands can still overflow in derived arithmetic and must be revalidated before storage.
- External graph node labels are identifiers, not allocation sizes. Sort unique labels, contract ideal wires with union-find, then map representatives densely before building adjacency or solver workspaces.
- Compatibility claims need a semantic boundary. Literal liionpack `Ri*` rows coexist with the SLIDE cell's own Thevenin resistance, while V/I magnitudes are runtime-owned; importing the graph schema does not establish waveform parity or lossless round-trip fidelity.

## Fuzz gates need independent publication oracles

- A sanitizer lane is not real until production and test compile commands carry the requested flags, linked executables contain the runtime symbols, and verbose CTest discovery points only into that lane's build directory. Shared source-root runtime outputs can make a perfectly green lane execute another configuration's binary.
- Embedded-NUL and hostile-byte fixtures must derive their view length from the backing array (`sizeof(array)-1`), never a hand-counted constant. The parser test itself otherwise becomes the out-of-bounds read that ASan reports.
- TSan owns global allocation interceptors. An allocation-counter executable that deliberately replaces global new/delete is structurally incompatible with that runtime; keep it in ASan/UBSan and run TSan on the actual concurrency suites, recording the exclusion explicitly.
- Instrument a separate library target. Adding sanitizer and Windows runtime/iterator flags privately to a shared static library still changes its object ABI and can break every non-fuzzer consumer in the same build graph.
- A deterministic double parse is not an atomicity proof. Seed the output with validly copyable poison across every public field, require byte/value identity on failure, and require every poison marker to disappear on success; otherwise append or merge publication can pass both parses.
- Do not call the production validator as the only success oracle. Independently recompute graph connectivity, cell bijection, sparsity, ladder classification, and cold thermal metadata so a correlated validator defect cannot bless its own output.
- Mutation-test the oracle, not only production guards. Deliberate Experiment append, BPX merge, and netlist scratch publication must crash on the first valid seed; a green mutation means the fuzzer is observing too little.
- A committed corpus must actually enter success paths. The Experiment byte adapter originally turned a normal final newline into an extra empty step, so its documented seed only exercised rejection until the adapter was corrected.
- Parser success must imply runner admissibility. A one-line voltage-control/voltage-event seed exposed syntax that parsed successfully but was rejected by `CyclerV2`; reuse the semantic descriptor validator before parser publication.
- Small mutation campaigns and resource boundaries are different tests. Keep fast 65,536-byte campaigns, then replay generated 65,537-byte and 4 MiB-plus inputs so documented limits are reached without committing multi-megabyte blobs.
- Disable C++ module dependency scanning at project scope when the project contains no modules. A target-only setting fixed fuzz compilation but not `FindThreads`: CMake 3.31 also scans its `try_compile` probes and otherwise requires the separately packaged `clang-scan-deps` executable.
- Scope platform caveats exactly. Windows LLVM required disabling incompatible MSVC STL container annotations, while full ASan+UBSan+LSan ran on Linux; cross-platform evidence should state that difference rather than blending the lanes.

## A thread pool is not integrated until production work overlaps

- A standalone pool and speed diagnostic do not satisfy an architecture that promises parallel batches. Use a bounded overlap oracle in the real Thevenin path; the original implementation was green in isolation while production stayed completely serial.
- Parallelize only ownership-disjoint archetype kernels. Gather inputs, scatter outputs, assemble thermal edges, reduce residuals, and cross substep boundaries in canonical order so worker scheduling cannot alter pack numerics.
- Prove ownership disjointness at configuration. Two distinct archetype slots are not evidence of two distinct objects; reject duplicate type-erased identities and duplicate concrete batch pointers before any worker can touch scratch or arena state.
- Make the execution context movable without moving live worker synchronization objects. A heap-owned persistent pool inside a movable `BatchExecutor` lets `PackSolver` retain exception-atomic candidate publication.
- “First failure” must name an order. Scheduling-first CAS is nondeterministic; lowest task index matches serial semantics and can be selected after all tasks finish.
- Source-order reduction syntax is not an order guarantee under fast-math. A long cancellation oracle made the unprotected loop return 5.0 instead of +0.0; strict-FP scope and volatile scalar stages are part of the implementation contract.
- Hardware discovery may legally return zero. Normalize that explicitly, cap workers to available tasks, and never let an exception cross a worker entry point: join first, then rethrow through the caller's channel.
- One worker and many workers need the same failure contract. The first legacy fix completed all tasks only in the threaded branch; a 1-vs-2 regression exposed the serial branch stopping at its first throw.
- Const ownership views should stay const. When diagnostics or benchmarks need one stateful operation, expose that operation narrowly on the owner instead of returning a mutable child that can invalidate topology and workspace invariants.

## Exact failure coverage needs source-level identity, not line coverage

- A failure arm can share a source line and region with success, especially for a conditional Status return. Ordinary line totals do not prove the failure token executed; scan the actual arm and intersect its token interval with LLVM segment counts.
- Do not combine unrelated standalone executables in one llvm-cov export: same-named inline or template functions can carry incompatible mapping hashes. Export each test against only its own profile, union counts by source-site identity, and use a whole-archive zero-count anchor to prove every production mapping exists.
- Bind a report to the code that produced it. Record source, compiled-source, manifest, compile-command, and binary identities; require one fresh profile per registered test; reject pending build work and mismatched Clang, llvm-cov, and llvm-profdata majors.
- An exception is a proof obligation, not an ignored line. Pin the exact site plus source/context hashes, allow only a named structural class and reason, and fail if the site becomes covered, inactive, unmapped, stale, or pushes the cap over budget.
- Finite positive resistance does not guarantee a usable conductance: `1/R` can underflow or flush to zero. Validate the derived coefficient used by the algebra, in every solver mode.
- Representable sizes are not allocation budgets. Bounded readers must account for simultaneously retained input, decoded payload, output, and scratch storage before allocation, and translate allocation or I/O failures without partially publishing state.

## Fast-math validation must classify the coefficient that actually exists

- A source-level `value > 0` is not a reliable guard under finite-math optimization. With FTZ enabled, `1 / DBL_MAX` materialises as +0, and Clang can delete the floating comparison; classify the derived coefficient through an opaque integer-bit barrier.
- Floating-environment tests must be scoped transactions. Save MXCSR, enable FTZ only for the discriminating call, and restore it through RAII even when an assertion aborts the section.
- Use exact power-of-two witnesses to isolate arithmetic order. `2^1023 * (step/8)` is representable for all eight source stages, while `(2^1023 * 2)/8` overflows at the intermediate; fast-math may reassociate the mutant safely, so a conforming Debug red result remains independent evidence.
- Resource-amplification tests should measure the largest requested block and throw before `malloc`, not rely on an OOM or final `Status`. A 128-byte wire header was proven to request 96,000,047 bytes without allocating it.
- A mutation that stays green is evidence against the claimed defect. Removing strict-FP from the already-separated PackSolver validation TU did not break the current Release test, so that policy remains defensive architecture rather than an extra ledger bug.
