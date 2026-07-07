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
