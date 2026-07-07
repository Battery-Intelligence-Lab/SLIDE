/**
 * @file Histogram_test.cpp
 * @brief Unit tests for slide::Histogram.
 *
 * Regression test for bug B2 (PLAN.md 2.4): a default-constructed Histogram has
 * Nbins == 0 and an EMPTY bins vector, yet add() computed a clamped index and wrote
 * bins[i] -> out-of-bounds write (heap corruption / UB). The static slide::EmptyHistogram
 * and the default-constructed cooling histograms (cool_data.hpp) are exactly this shape.
 * The fix makes add() a no-op when there are no bins; this test both proves the empty
 * case is safe and that a properly-constructed histogram still bins correctly.
 */

#include "../../src/types/Histogram.hpp"

#include <catch2/catch_test_macros.hpp>

using slide::Histogram;

TEST_CASE("Default-constructed Histogram: add() is a safe no-op", "[Histogram][B2]")
{
  Histogram<> h; //!< empty bins (Nbins == 0)
  REQUIRE(h.size() == 0);

  //!< Pre-fix these calls perform bins[0]++ on an empty vector (OOB write); with the
  //!< MSVC debug STL (_DEBUG => _ITERATOR_DEBUG_LEVEL 2) that aborts the process.
  for (double x : { -10.0, 0.0, 0.5, 1.0, 42.0 })
    h.add(x);

  REQUIRE(h.size() == 0); //!< still empty, nothing recorded, no crash
}

TEST_CASE("EmptyHistogram global is safe to add() into", "[Histogram][B2]")
{
  REQUIRE(slide::EmptyHistogram.size() == 0);
  slide::EmptyHistogram.add(3.14); //!< must not corrupt the heap
  REQUIRE(slide::EmptyHistogram.size() == 0);
}

TEST_CASE("Properly-constructed Histogram bins values correctly", "[Histogram][B2]")
{
  //!< 10 equidistant bins over [0,10] -> dx = 1; internal storage is Nbins+2 = 12
  //!< (index 0 = underflow, index end = overflow).
  Histogram<> h(0.0, 10.0, 10);
  REQUIRE(h.size() == 12);

  h.add(-5.0); //!< underflow -> bins[0]
  h.add(3.5);  //!< 1 + floor(3.5) = bin 4
  h.add(3.9);  //!< same bin 4
  h.add(100.0); //!< overflow -> last bin (index 11)

  auto bins = h.viewBinValues();
  REQUIRE(bins[0] == 1);
  REQUIRE(bins[4] == 2);
  REQUIRE(bins[11] == 1);

  //!< total count preserved
  size_t total = 0;
  for (auto b : bins) total += b;
  REQUIRE(total == 4);
}
