# Phase 0 (agent B) changelog fragment

Keep-a-Changelog style bullets for the CHANGELOG.md `Unreleased` section.
The architect merges these into CHANGELOG.md (agent A owns that file).

### Fixed

- Time-series cell data storage (`CellDataStorage<storeTimeData>`) now **appends** each
  timestep instead of dropping history: replaced the ill-formed `data.assign(data.end(), {...})`
  (an iterator + initializer-list call that matched no `std::vector::assign` overload) with
  `data.insert(data.end(), {...})`. Also moved the `#include "CellDataWriter.hpp"` out of the
  `namespace slide { ... }` block to global scope; the previous mid-namespace include pulled
  `<variant>` and other standard headers into `namespace slide` (corrupting them) and nested
  `slide::slide::CellDataWriter`, so the header could not be compiled at all when included.
  (`src/recording/CellDataStorage.hpp`)
- `Histogram::add()` no longer performs an out-of-bounds write on a default-constructed
  (empty-bins, `Nbins == 0`) histogram: it now returns early when there are no bins. This
  affects the global `slide::EmptyHistogram` and any histogram used before `initialise()`.
  (`src/types/Histogram.hpp`)
- `Cycler::CC` energy throughput (`ThroughputData::Wh`) now integrates with the trapezoid
  rule `(v_before + v_after)/2` over each step, sampling the terminal voltage before and
  after `timeStep_CC`. Previously it used a single `vi` sample that `Cycler::setCurrent`
  never assigns, so the CC energy throughput was always 0 Wh. (`src/procedures/Cycler.cpp`)
- Corrected the SEI degradation-model documentation to match the code (the source of truth):
  the unreachable-model error message now says "0 to 4" (cases 0-4 exist, not 0-3); the
  inline comments for SEI ids 1 and 2 had their literature references swapped relative to
  the actual formulas (id 1 = kinetics-limited/Ning & Popov, id 2 = kinetics + SEI-layer
  diffusion/Pinson & Bazant); and `DEG_ID`'s SEI list now documents id 4. No behavioural
  change. (`src/cells/Cell_SPM/Cell_SPM_degradation.cpp`, `src/cells/Cell_SPM/param/DEG_ID.hpp`)
- Documented and debug-guarded the invariant behind the degradation forward-Euler loop in
  `Cell_SPM::timeStep_CC`: the loop integrates every state index (including the algebraic
  current/voltage slots) and is only correct because `dState_degradation` leaves `d_st[i_I]`
  and `d_st[i_V]` at 0. Added `assert(d_st.I() == 0.0 && d_st.V() == 0.0)`. Digit-identical
  behaviour (Release: assert removed; Debug: invariant holds). (`src/cells/Cell_SPM/Cell_SPM_dstate.cpp`)
