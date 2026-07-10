# Handoff — Phase 9B parser hardening (2026-07-10)

## Outcome

The Experiment and BPX parser adversarial pass is complete. P9-B26 through
P9-B28 are fixed and recorded. Phase 9B remains active: P9-G2 still needs the
netlist CSV parser plus all three bounded fuzz campaigns, and P9-G1/G3/G4 are
not yet closed.

## Changes

- Experiment parsing now caps 10,000 expanded segments, 4 MiB of aggregate
  retained text, 65,536 bytes per step, and 1,024 bytes per drive-cycle name.
  All arithmetic is checked before allocation; explicit string grammar replaces
  two regexes without changing documented inputs.
- BPX JSON now enforces strict number, whitespace, and UTF-8 grammar, caps the
  source at 4 MiB and the parsed tree at 65,536 values, and validates exact
  bounded file reads through EOF.
- Every required, optional, derived, activation, curve, state, and default
  `ParameterSet` insertion propagates its exact `Status`. Present-invalid
  optional values no longer disappear into defaults.
- Parser, direct `ParameterSet::set`, and BPX file allocations are
  exception-translating transactions. Best-effort diagnostic assignment remains
  safe when allocation continues failing.
- The shared IEEE finite classifier now takes a reference boundary (with a
  volatile overload) so Release finite-math cannot attach a by-value finite-only
  assumption before bit classification.

## Mutation and validation evidence

- Restoring the old optional-field `Numerical_failure` swallow fails the
  allocation test at matching occurrence 149 of 176 and publishes 29 entries.
- Removing the drive-cycle minimum-length guard makes `Run (A)` succeed and
  fails 3 registered grammar assertions.
- Debug and Release: Experiment 199/199, ParameterSet 1,495/1,495, parser
  allocation 906/906, ThermalLumped 25/25, CompiledCurve 163/163, PackSolver
  667/667, and PackStepper 172/172.
- Allocation faults are persistent after the selected owned allocation, so the
  same tests exercise diagnostic double faults. Arbitrary global allocation
  ordinals were rejected because one enters a Windows CRT/library abort dialog
  instead of testing the parser transaction.
- No long simulation ran. The full Experiment binary contains only bounded unit
  scenarios and was preceded by a static simulation-facing delta review.

## Next

Implement the missing liionpack-compatible netlist CSV parser with an atomic
compiled-topology output, a committed hostile corpus, and a bounded fuzz driver.
Then add matching Experiment/BPX fuzz drivers and CI wiring to close P9-G2
before moving to the ThreadPool production-integration defects.
