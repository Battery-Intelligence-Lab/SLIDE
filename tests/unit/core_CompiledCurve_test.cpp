/**
 * @file core_CompiledCurve_test.cpp
 * @brief Tests for build-time curve canonicalisation (PLAN.md §3.11, D-16).
 */

#include "../../src/core/CompiledCurve.hpp"
#include "../../src/types/OCVcurves.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>

using namespace slide;

TEST_CASE("IndexedPiecewiseLinear exactly reproduces legacy nonuniform OCV interpolation",
          "[core][parameters]")
{
  auto curves = OCVcurves::makeOCVcurves(cellType::KokamNMC);
  core::IndexedPiecewiseLinear compiled;
  REQUIRE(compiled.build(curves.OCV_neg.x, curves.OCV_neg.y) == Status::Success);

  const std::array queries{ 0.037, 0.123456789, 0.5, 0.87654321, 0.963 };
  for (const double x : queries) {
    const double legacy = curves.OCV_neg.interp(x, false, true);
    const double value = compiled.eval(x);
    CAPTURE(x, legacy, value);
    REQUIRE(value == legacy);
  }
}

TEST_CASE("UniformLut validates a smooth injected curve to the registered tolerance",
          "[core][parameters]")
{
  constexpr std::size_t count = 65;
  std::array<double, count> x{}, y{};
  for (std::size_t i = 0; i < count; ++i) {
    x[i] = static_cast<double>(i) / static_cast<double>(count - 1);
    y[i] = std::exp(0.5 * x[i]);
  }

  core::UniformLut lut;
  REQUIRE(lut.build(x, y, 1e-6, 4096) == Status::Success);
  REQUIRE(lut.valid());
  REQUIRE(lut.max_relative_error() <= 1e-6);

  for (int i = 0; i <= 100; ++i) {
    const double query = static_cast<double>(i) / 100.0;
    const double expected = std::exp(0.5 * query);
    REQUIRE(std::abs(lut.eval(query) - expected) / expected <= 4e-5);
  }
}

TEST_CASE("Compiled curves reject malformed input without becoming valid", "[core][parameters]")
{
  const std::array valid_x{ 0.0, 0.5, 1.0 };
  const std::array valid_y{ 1.0, 2.0, 3.0 };
  const std::array x{ 0.0, 0.5, 0.4, 1.0 };
  const std::array y{ 1.0, 2.0, 3.0, 4.0 };
  core::IndexedPiecewiseLinear curve;
  core::UniformLut lut;

  REQUIRE(curve.build(valid_x, valid_y) == Status::Success);
  REQUIRE(lut.build(valid_x, valid_y) == Status::Success);
  REQUIRE(curve.build(x, y) == Status::Invalid_parameters);
  REQUIRE(lut.build(x, y) == Status::Invalid_parameters);
  REQUIRE_FALSE(curve.valid());
  REQUIRE_FALSE(lut.valid());

  const double nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  const std::array nan_y{ 1.0, nan, 3.0 };
  REQUIRE(curve.build(valid_x, nan_y) == Status::Invalid_parameters);
  REQUIRE(lut.build(valid_x, nan_y) == Status::Invalid_parameters);
}
