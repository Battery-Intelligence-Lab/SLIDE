# SLIDE Engineering Lessons

This file records reusable findings from the v4 implementation and validation work. Detailed design decisions remain in `PLAN.md`; session chronology remains in `.claude/discussions.md` and `.claude/summaries/`.

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
- Exception atomicity is easiest when every throwing copy/allocation happens in locals before the first public member is committed. Move/swap only through the no-throw commit phase; deterministic allocation faults should also prove the object can be configured successfully afterward.

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
