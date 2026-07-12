/**
 * @file core_SpmScalarKernels_test.cpp
 * @brief Direct scalar and Dual gates for the PC-10 single-source SPM algebra.
 */

#include "../../src/core/Dual.hpp"
#include "../../src/core/SpmScalarKernels.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>

using namespace slide;

namespace {

double independentModal(double z, double x, double input, double flux)
{
  const long double xl = x;
  const long double phi = xl == 0.0L ? 1.0L : std::expm1(xl) / xl;
  return static_cast<double>(std::exp(xl) * static_cast<long double>(z)
                             + phi * static_cast<long double>(input)
                                 * static_cast<long double>(flux));
}

template <class Scalar>
auto observableChain(Scalar reaction_reference,
                     Scalar concentration,
                     Scalar current_density,
                     Scalar collector_resistance_area)
{
  using namespace core::spm_scalar;
  const double temperature = 310.0;
  const double faraday = 96487.0;
  const double arrhenius = arrheniusFactor(298.15, temperature, 8.314);
  const auto reaction_rate = activatedValue(reaction_reference, 35000.0, arrhenius);
  const auto exchange = exchangeCurrent(reaction_rate,
                                        1.0,
                                        faraday,
                                        1000.0,
                                        concentration,
                                        30'000.0);
  const auto argument = activationArgument(-1.0,
                                           current_density,
                                           180'000.0,
                                           70e-6,
                                           exchange);
  const auto eta_negative = activationOverpotential(temperature,
                                                    8.314,
                                                    1.0,
                                                    faraday,
                                                    argument);
  const double negative_area = activeArea(180'000.0, 0.1027, 70e-6);
  const double positive_area = activeArea(150'000.0, 0.1027, 75e-6);
  const auto resistance = seriesResistance(2e-9,
                                           5e5,
                                           2e-4,
                                           3e-4,
                                           collector_resistance_area,
                                           negative_area,
                                           positive_area,
                                           0.1027);
  const double ocv = cellOpenCircuitVoltage(0.12,
                                            4.08,
                                            temperature,
                                            298.15,
                                            1.4e-4);
  return terminalVoltage(ocv,
                         eta_negative,
                         0.035,
                         resistance,
                         current_density * 0.1027);
}

} // namespace

TEST_CASE("PC-10 scalar modal kernel covers zero, Taylor, and expm1 branches",
          "[core][PC-10][scalar][modal]")
{
  using core::spm_scalar::advanceModal;
  constexpr double z = 0.725;
  constexpr double input = 1.125;
  constexpr double flux = -0.0375;
  const double below = std::nextafter(1e-7, 0.0);
  const double above = std::nextafter(1e-7, 2e-7);
  for (const double x : { 0.0,
                          -0.0,
                          below,
                          -below,
                          1e-7,
                          -1e-7,
                          above,
                          -above,
                          -0.35 }) {
    const double actual = advanceModal(z, 1.0, x, 1.0, input, flux);
    const double expected = independentModal(z, x, input, flux);
    CAPTURE(x, actual, expected);
    REQUIRE(actual == Catch::Approx(expected).epsilon(8e-15).margin(2e-16));
  }
  REQUIRE(advanceModal(z, 1.0, 0.0, 1.0, input, flux)
          == z + input * flux);
}

TEST_CASE("PC-10 scalar kernels propagate Dual tangents through both modal branches",
          "[core][PC-10][scalar][dual]")
{
  auto check = [](double diffusivity) {
    const core::Dual z{ 0.7, 1.2 };
    const core::Dual d{ diffusivity, -0.004 };
    const core::Dual flux{ -0.2, 0.05 };
    const core::Dual actual = core::spm_scalar::advanceModal(
      z, d, -2.0, 0.4, 1.1, flux);
    auto primal = [=](double t) {
      return core::spm_scalar::advanceModal(0.7 + 1.2 * t,
                                            diffusivity - 0.004 * t,
                                            -2.0,
                                            0.4,
                                            1.1,
                                            -0.2 + 0.05 * t);
    };
    constexpr double h = 1e-6;
    const double finite_difference = (primal(h) - primal(-h)) / (2.0 * h);
    CAPTURE(diffusivity, actual.value, actual.derivative, finite_difference);
    REQUIRE(actual.value == Catch::Approx(primal(0.0)).epsilon(2e-15));
    REQUIRE(actual.derivative
            == Catch::Approx(finite_difference).epsilon(2e-8).margin(2e-10));
  };
  check(1e-8); // Taylor branch
  check(0.03); // expm1 branch
}

TEST_CASE("PC-10 shared observable leaves preserve mixed scalar algebra",
          "[core][PC-10][scalar][observable]")
{
  constexpr double reaction_reference = 2.4e-12;
  constexpr double concentration = 14'000.0;
  constexpr double current_density = 22.0;
  constexpr double collector = 1.7e-4;
  const core::Dual actual = observableChain(core::Dual{ reaction_reference, 0.3e-12 },
                                            core::Dual{ concentration, 125.0 },
                                            core::Dual{ current_density, -0.4 },
                                            core::Dual{ collector, 2e-5 });
  auto primal = [](double t) {
    return observableChain(reaction_reference + 0.3e-12 * t,
                           concentration + 125.0 * t,
                           current_density - 0.4 * t,
                           collector + 2e-5 * t);
  };
  constexpr double h = 1e-5;
  const double finite_difference = (primal(h) - primal(-h)) / (2.0 * h);
  REQUIRE(actual.value == Catch::Approx(primal(0.0)).epsilon(3e-15));
  REQUIRE(actual.derivative
          == Catch::Approx(finite_difference).epsilon(2e-7).margin(2e-10));

  const std::array<core::Dual, 3> modes{
    core::Dual{ 1.0, 0.1 }, core::Dual{ -2.0, 0.2 }, core::Dual{ 0.5, -0.3 }
  };
  constexpr std::array coefficients{ 0.25, -0.5, 2.0 };
  const core::Dual surface = core::spm_scalar::concentrationOutput(
    3, modes.data(), 1, coefficients.data(), 0.4, core::Dual{ 0.2, 0.03 }, core::Dual{ 0.8, -0.02 });
  REQUIRE(surface.value == Catch::Approx(2.35).epsilon(2e-15));
  REQUIRE(surface.derivative == Catch::Approx(-0.6575).epsilon(2e-15));

  const core::Dual interpolated = core::spm_scalar::linearInterpolate(
    core::Dual{ 0.4, 1.0 }, 0.0, 1.0, 2.0, 5.0);
  REQUIRE(interpolated.value == Catch::Approx(3.2));
  REQUIRE(interpolated.derivative == Catch::Approx(3.0));
  REQUIRE(core::spm_scalar::diffusionRate(2.0, 0.3, -4.0, 1.5, -0.2)
          == Catch::Approx(-2.7));
}
