---
layout: default
title: Add a cell model
nav_order: 6
---

# Add a v4 cell-model family

There is no public runtime model-plugin ABI in v4.0. `AgingMechanismSpec` and the cold description types are metadata, not dynamic registration hooks. A new model family is a contributor change that must preserve the compiled-batch and pack contracts.

If only the chemistry or geometry changes, do **not** add a model family. Supply a BPX file or extend `ParameterSet` absorption, then compile the existing SPM family.

For genuinely different equations (for example SPMe or DFN), use this sequence:

1. **Write the equations and falsifiers first.** Define states, algebraic variables, conservation laws, units, admissible ranges, an independent oracle, and the shortest non-degenerate validation case.
2. **Add cold value types.** Extend or parallel `CellDesign`/`ElectrodeDesign` with named parameter packs. No file I/O, strings, maps, or expression evaluation may survive into the hot kernel.
3. **Declare the state layout.** Use `BatchBuilder` and `StateSpec` for every value read by a later step, including warm history and cumulative state. Choose the correct `StateRole` for ODE, algebraic, cumulative, and input rows.
4. **Implement scalar-generic kernels.** Operate on `BatchView`, `RhsViews`, and `StepCtx`; add into a caller-zeroed derivative. Keep observables in one shared path used by physics, recording, and sensitivities.
5. **Compose a concrete batch.** Instantiate the kernel family at compile time and register one cold selection per supported discretisation/composition. Inner lane loops must see concrete types, not virtual calls or `std::function`.
6. **Implement the pack boundary.** Provide the frozen-state terminal voltage and Thevenin linearisation required by `PackSolver`; keep topology and solver code independent of cell type.
7. **Expose all languages from the same selection.** Map C++, Python, and MATLAB options to the same registry entry. Reject unsupported options before allocation or stepping.
8. **Validate in layers.** Unit-test operators and conservation, compare one batch with an independent oracle, add rollback/restart and allocation gates, test heterogeneous lanes and degenerate counts, then add language parity.

Before merging, re-check the `PLAN.md` performance contract: zero accepted-step allocations, contiguous SoA state, no per-cell virtual dispatch, no exceptions/locks/I/O on the hot path, and memcpy rollback.

The legacy v3 `Cell` subclass instructions remain in the older Advanced pages only for façade maintenance; editing `Cell_user` or `main.cpp` does not register a v4 model.
