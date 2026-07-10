# Handoff — Phase 9B production ThreadPool audit (2026-07-10)

## Outcome

The pool is now a production batch runtime rather than a diagnostic-only class.
`PackTheveninSystem` linearizes independent archetypes concurrently and
`PackStepper` advances independent batches concurrently behind a barrier for
every substep. All cross-batch gathers, scatters, thermal assembly, reductions,
and publication remain canonical and serial.

## Architecture

- D-30 introduces movable `BatchExecutor`: it owns at most one heap-stable
  persistent `ThreadPool`, caps requested workers to batch count, reports the
  selected count, and executes inline when exactly one worker is selected.
- Thevenin views expose an opaque object identity and PackStepper checks concrete
  pointers; duplicate batch handles reject during cold configuration, preventing
  concurrent access to one scratch/cache or arena.
- `PackSolver` owns the executor, so standalone solves and composed pack steps
  share the same production path. Thread creation occurs during exception-atomic
  configuration and maps failure to `Status`; accepted steps allocate nothing.
- The pool returns the failure from the lowest task index after all tasks finish.
  A typed stack context preserves const callbacks without casting away const.
- `fixedOrderSum` disables reassociation/contraction locally and uses volatile
  scalar stages. This is required by the repository's fast-math Release build.
- `Threads::Threads` is a public core link requirement. The retained v3
  `slide::run` is still a per-call strangler facade, but zero-worker discovery,
  partial construction, and escaping worker exceptions are now safe.
- The repository declares its no-C++-modules policy at project scope. This is
  required for CMake 3.31's Clang `FindThreads` `try_compile` probes as well as
  ordinary targets; a target-only property still required `clang-scan-deps`.

## Adversarial evidence

- The deterministic-failure test was red before the fix: index 1 published
  first, but index 0 was the required serial-order failure.
- A standalone const-lvalue callback compile probe failed before and succeeds
  after the typed context change.
- The production overlap test was red with `max_active=1`; it now observes at
  least two simultaneous Thevenin archetype calls.
- Removing strict-FP protection makes the optimized cancellation oracle return
  bit pattern `0x4014000000000000` (5.0) rather than +0.0.
- Pack solutions are bit-identical with 1, 2, and 7 workers.
- Four real two-archetype thermal pack steps are byte-identical with 1 and 2
  workers across both arenas, solution, diagnostics, cell heat, and boundary heat.
- A parallel two-archetype accepted step changes the global allocation counter
  by exactly zero after warm-up.

Debug and fast-math Release each pass ThreadPool 600, PackSolver 704,
PackStepper 208, and P2-G1 allocation 14 assertions.
Fresh WSL Ubuntu/Clang 18 core-only configure/build/smoke passes 1/1 with
`Threads::Threads` discovered and propagated to the public consumer.

## Next

First fault-inject standalone `PackSolver::configure`: independent review found
that it publishes executor/topology members before later scratch-vector
allocations, so a thrown allocation may violate reconfiguration atomicity.
After that focused transaction gate, P9-G1 remains open: build isolated
instrumented core copies and run the full Debug suite under ASan+UBSan plus
serial focused ThreadPool/AsyncRecorder/PackStepper tests under TSan. Do not
claim hosted results until pushed.
