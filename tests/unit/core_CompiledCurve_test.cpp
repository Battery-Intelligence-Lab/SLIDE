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
#include <limits>

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

TEST_CASE("Indexed curves admit BPX adaptive endpoint resolution",
          "[core][parameters][BPX]")
{
  const std::array x{ 0.0, 1.0 / 65'536.0, 0.5, 1.0 };
  const std::array y{ 1.0, 1.1, 2.0, 3.0 };
  core::IndexedPiecewiseLinear curve;
  REQUIRE(curve.build(x, y) == Status::Success);
  CHECK(curve.eval(x[1]) == y[1]);
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

TEST_CASE("Compiled curves reject unrepresentable accelerators before integer conversion",
          "[core][parameters][curve][P9]")
{
  const double maximum = std::numeric_limits<double>::max();
  const std::array extreme_x{ -maximum, 0.0, maximum };
  const std::array ordinary_y{ 1.0, 2.0, 3.0 };
  const std::array ordinary_x{ 0.0, 1.0 };
  const std::array extreme_y{ -maximum, maximum };
  const double denormal = std::bit_cast<double>(UINT64_C(1));
  const double twice_denormal = std::bit_cast<double>(UINT64_C(2));
  const std::array tiny_x{ 0.0, denormal, twice_denormal };

  for (const auto &x : { extreme_x, tiny_x }) {
    core::IndexedPiecewiseLinear curve;
    core::UniformLut lut;
    CHECK(curve.build(x, ordinary_y) == Status::Invalid_parameters);
    CHECK(lut.build(x, ordinary_y) == Status::Invalid_parameters);
    CHECK_FALSE(curve.valid());
    CHECK_FALSE(lut.valid());
  }

  core::IndexedPiecewiseLinear curve;
  core::UniformLut lut;
  CHECK(curve.build(ordinary_x, extreme_y) == Status::Invalid_parameters);
  CHECK(lut.build(ordinary_x, extreme_y) == Status::Invalid_parameters);
}

TEST_CASE("Compiled curve domain edges never cast NaN to an index",
          "[core][parameters][curve][P9]")
{
  const std::array x{ 0.0, 1.0 };
  const std::array y{ 2.0, 4.0 };
  core::IndexedPiecewiseLinear curve;
  core::UniformLut lut;
  REQUIRE(curve.build(x, y) == Status::Success);
  REQUIRE(lut.build(x, y) == Status::Success);

  CHECK(curve.eval(-1.0) == y.front());
  CHECK(curve.eval(0.0) == y.front());
  CHECK(curve.eval(1.0) == y.back());
  CHECK(curve.eval(2.0) == y.back());
  CHECK(curve.derivative(0.0) == 2.0);
  CHECK(curve.derivative(1.0) == 2.0);

  const double nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  core::UniformLut invalid_tolerance_lut;
  CHECK(invalid_tolerance_lut.build(x, y, nan) == Status::Invalid_parameters);
  CHECK_FALSE(invalid_tolerance_lut.valid());

  SECTION("indexed evaluation")
  {
    const double result = curve.eval(nan);
    const auto bits = std::bit_cast<std::uint64_t>(result);
    CAPTURE(bits);
    CHECK((bits & exponent_mask) == exponent_mask);
  }
  SECTION("indexed derivative")
  {
    const double result = curve.derivative(nan);
    const auto bits = std::bit_cast<std::uint64_t>(result);
    CAPTURE(bits);
    CHECK((bits & exponent_mask) == exponent_mask);
  }
  SECTION("uniform evaluation")
  {
    const double result = lut.eval(nan);
    const auto bits = std::bit_cast<std::uint64_t>(result);
    CAPTURE(bits);
    CHECK((bits & exponent_mask) == exponent_mask);
  }
}

TEST_CASE("Compiled curve builders classify every representable arithmetic boundary",
          "[core][parameters][curve][coverage]")
{
  const double nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  const double denormal = std::bit_cast<double>(UINT64_C(1));
  const double maximum = std::numeric_limits<double>::max();

  SECTION("the first knot is validated independently")
  {
    const std::array x{ nan, 1.0 };
    const std::array y{ 0.0, 1.0 };
    core::IndexedPiecewiseLinear curve;
    CHECK(curve.build(x, y) == Status::Invalid_parameters);
  }

  SECTION("a finite range divided by a denormal spacing is rejected")
  {
    const std::array x{ 0.0, denormal, 1.0 };
    const std::array y{ 1.0, 1.0, 1.0 };
    core::IndexedPiecewiseLinear curve;
    CHECK(curve.build(x, y) == Status::Invalid_parameters);
  }

  SECTION("an unrepresentable indexed reciprocal is rejected")
  {
    const std::array x{ 0.0, denormal };
    const std::array y{ 1.0, 1.0 };
    core::IndexedPiecewiseLinear curve;
    CHECK(curve.build(x, y) == Status::Invalid_parameters);
  }

  SECTION("the finer uniform reciprocal is checked separately")
  {
    const std::array x{ 0.0, 1e-305 };
    const std::array y{ 1.0, 1.0 };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 1e-6, 4096) == Status::Invalid_parameters);
  }

  SECTION("sampling detects overflow in the legacy interpolation order")
  {
    const std::array x{ 0.0, 1e308 };
    const std::array y{ 0.0, 1e308 };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 1e-6, 3) == Status::Numerical_failure);
  }

  SECTION("sampling rejects rounded addition beyond finite endpoints")
  {
    const double xmin = std::bit_cast<double>(UINT64_C(0x7fc665fc6bf52b9e));
    const std::array x{ xmin, maximum };
    const std::array y{ 1.0, 1.0 };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 1e-6, 2) == Status::Numerical_failure);
  }

  SECTION("knot validation rejects an unrepresentable relative error")
  {
    const std::array x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
    const std::array y{ maximum, 0.0, -maximum, 0.0, maximum };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 1e-6, 2) == Status::Numerical_failure);
  }

  SECTION("midpoint validation independently rejects interpolation overflow")
  {
    const std::array x{ 0.0, 1e308 };
    const std::array y{ 0.0, 1e308 };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 1e-6, 2) == Status::Numerical_failure);
  }

  SECTION("a finite but inaccurate coarse LUT fails its requested tolerance")
  {
    const std::array x{ 0.0, 0.5, 1.0 };
    const std::array y{ 0.0, 1.0, 0.0 };
    core::UniformLut lut;
    CHECK(lut.build(x, y, 0.1, 2) == Status::Numerical_failure);
  }
}
