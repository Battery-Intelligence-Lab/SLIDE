/**
 * @file CellDataStorage_test.cpp
 * @brief Unit tests for the time-series specialization of CellDataStorage.
 *
 * Regression test for bug B1 (PLAN.md 2.4): the storeTimeData specialization used
 * `data.assign(data.end(), {...})`, which (a) is ill-formed (iterator + init-list
 * cannot bind the assign(InputIt,InputIt) overload) and (b) would REPLACE rather
 * than APPEND the time-series history. This test instantiates the specialization
 * (nothing else in the codebase does — Cell::storeData is a no-op and the cellData
 * member is commented out at Cell.hpp:40) and verifies each timestep is appended.
 */

#include "../../src/recording/CellDataStorage.hpp"
#include "../../src/types/State.hpp" // ThroughputData = State<0>

#include <catch2/catch_test_macros.hpp>

namespace {

//!< Minimal duck-typed stand-in for a Cell, exposing exactly what
//!< CellDataStorage<storeTimeData>::storeData reads.
struct FakeCell
{
  double i{ 1.5 }, v{ 3.7 }, soc{ 0.42 }, t{ 298.15 };
  double I() const { return i; }
  double V() const { return v; }
  double SOC() const { return soc; }
  double T() const { return t; }
  slide::ThroughputData getThroughputs() const
  {
    slide::ThroughputData th{};
    th.time() = 11.0;
    th.Ah() = 2.0;
    th.Wh() = 7.0;
    return th;
  }
};

} // namespace

using slide::settings::CellDataStorageLevel;

TEST_CASE("CellDataStorage<storeTimeData> appends 7 values per timestep", "[recording][B1]")
{
  slide::CellDataStorage<CellDataStorageLevel::storeTimeData> store;
  FakeCell cell;

  REQUIRE(store.data.empty());

  constexpr int nsteps = 3;
  constexpr int nfields = 7; //!< I, V, SOC, T, time, Ah, Wh
  for (int k = 0; k < nsteps; ++k)
    store.storeData(cell);

  //!< All three timesteps must survive (the buggy `assign` would leave <21 entries).
  REQUIRE(store.data.size() == static_cast<size_t>(nsteps * nfields));

  //!< First record's contents must be intact (would be clobbered by `assign`).
  REQUIRE(store.data[0] == cell.I());
  REQUIRE(store.data[1] == cell.V());
  REQUIRE(store.data[2] == cell.SOC());
  REQUIRE(store.data[3] == cell.T());
  REQUIRE(store.data[4] == cell.getThroughputs().time());
  REQUIRE(store.data[5] == cell.getThroughputs().Ah());
  REQUIRE(store.data[6] == cell.getThroughputs().Wh());

  //!< Third record starts at offset 14 and is identical (append, not overwrite).
  REQUIRE(store.data[14] == cell.I());
  REQUIRE(store.data[20] == cell.getThroughputs().Wh());
}
