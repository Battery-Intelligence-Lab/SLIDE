/**
 * @file core_SpectralModel_test.cpp
 * @brief Cold spectral build gates and legacy-default coefficient parity.
 */

#include "../../src/slide.hpp"
#include "../../src/core/SpectralModel.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

using namespace slide;

namespace {

template <int NCH>
void require_registry_order_builds()
{
  core::CompiledSpectralModel<NCH> model;
  const core::PerDomain<double> radius{ 12.5e-6, 8.5e-6 };
  REQUIRE(core::compileSpectralModel<NCH>(radius, model) == Status::Success);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    REQUIRE(std::count(model.A[d].begin(), model.A[d].end(), 0.0) == 1);
    for (const double value : model.A[d])
      REQUIRE(value <= 0.0);
  }
}

} // namespace

TEST_CASE("Spectral compiler validates every registered NCH order", "[core][factory][spectral]")
{
  require_registry_order_builds<5>();
  require_registry_order_builds<8>();
  require_registry_order_builds<12>();
}

TEST_CASE("Spectral compiler preserves the legacy default model coefficients",
          "[core][factory][spectral]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  static_assert(NCH == 5);
  core::CompiledSpectralModel<NCH> actual;
  const core::PerDomain<double> radius{ 12.5e-6, 8.5e-6 };
  REQUIRE(core::compileSpectralModel<NCH>(radius, actual) == Status::Success);
  const auto &legacy = *Model_SPM<NCH>::makeModel();

  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    const auto legacy_domain = domain == core::Domain::neg ? neg : pos;
    REQUIRE(actual.zero_mode[d] == legacy.zero);
    for (int mode = 0; mode < NCH; ++mode) {
      CAPTURE(d, mode);
      REQUIRE(actual.A[d][static_cast<std::size_t>(mode)] == legacy.A[legacy_domain](mode));
      REQUIRE(actual.B[d][static_cast<std::size_t>(mode)] == legacy.B[legacy_domain](mode));
      REQUIRE(actual.x_inner[static_cast<std::size_t>(mode)] == legacy.xch(mode));
      for (int column = 0; column < NCH; ++column)
        REQUIRE(actual.state_transform[d][static_cast<std::size_t>(mode)]
                                      [static_cast<std::size_t>(column)]
                == legacy.V[legacy_domain](mode, column));
    }
    for (int row = 0; row < NCH + 1; ++row) {
      REQUIRE(actual.D[d][static_cast<std::size_t>(row)] == legacy.D[legacy_domain](row));
      for (int mode = 0; mode < NCH; ++mode)
        REQUIRE(actual.C[d][static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(mode)]
                == legacy.C[legacy_domain](row, mode));
    }
  }
  REQUIRE(actual.cc_coeff == legacy.cc_coeff);
  for (int row = 0; row < NCH + 1; ++row)
    REQUIRE(actual.Cc[static_cast<std::size_t>(row)] == legacy.Cc(row));
  for (int row = 0; row < 2 * NCH + 3; ++row)
    for (int column = 0; column < 2 * NCH + 3; ++column)
      REQUIRE(actual.integration[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)]
              == legacy.Q(row, column));
}

TEST_CASE("Spectral compiler scales geometry and fails atomically", "[core][factory][spectral]")
{
  core::CompiledSpectralModel<5> base;
  core::CompiledSpectralModel<5> scaled;
  const core::PerDomain<double> base_radius{ 12.5e-6, 8.5e-6 };
  const core::PerDomain<double> scaled_radius{ 25.0e-6, 17.0e-6 };
  REQUIRE(core::compileSpectralModel<5>(base_radius, base) == Status::Success);
  REQUIRE(core::compileSpectralModel<5>(scaled_radius, scaled) == Status::Success);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int mode = 0; mode < 5; ++mode) {
      if (base.A[d][static_cast<std::size_t>(mode)] != 0.0)
        REQUIRE(std::abs(scaled.A[d][static_cast<std::size_t>(mode)]
                           / base.A[d][static_cast<std::size_t>(mode)]
                         - 0.25)
                <= 1e-12);
    }
    REQUIRE(scaled.D[d][0] == 2.0 * base.D[d][0]);
  }

  scaled.A[0][0] = 42.0;
  const core::PerDomain<double> invalid_radius{ 0.0, 8.5e-6 };
  REQUIRE(core::compileSpectralModel<5>(invalid_radius, scaled)
          == Status::Invalid_parameters);
  REQUIRE(scaled.A[0][0] == 0.0);
}

TEST_CASE("Spectral validation rejects each malformed eigen-spectrum class",
          "[core][factory][spectral][validation]")
{
  int zero_mode = -1;
  Eigen::Vector<double, 5> spectrum;

  spectrum << 1e-7, -1.0, -4.0, -9.0, -16.0;
  REQUIRE(core::detail::validateSpectrum<5>(spectrum, 1.0, zero_mode)
          == Status::Numerical_failure);

  spectrum << 0.0, 1.0, -4.0, -9.0, -16.0;
  REQUIRE(core::detail::validateSpectrum<5>(spectrum, 1.0, zero_mode)
          == Status::Numerical_failure);

  spectrum << 0.0, -1.0, -4.0, -9.0, -16.0;
  REQUIRE(core::detail::validateSpectrum<5>(spectrum, 1.0, zero_mode)
          == Status::Numerical_failure);
}

TEST_CASE("Spectral compiler rejects unrepresentable radius scales atomically",
          "[core][factory][spectral][validation]")
{
  core::CompiledSpectralModel<5> output;
  output.A[0][0] = 42.0;
  const double maximum = std::numeric_limits<double>::max();
  REQUIRE(core::compileSpectralModel<5>({ maximum, maximum }, output)
          == Status::Numerical_failure);
  REQUIRE(output.A[0][0] == 0.0);

  output.A[0][0] = 42.0;
  const double minimum = std::numeric_limits<double>::min();
  REQUIRE(core::compileSpectralModel<5>({ minimum, minimum }, output)
          == Status::Numerical_failure);
  REQUIRE(output.A[0][0] == 0.0);
}

TEST_CASE("Compiled spectral finiteness validation covers every coefficient family",
          "[core][factory][spectral][validation]")
{
  core::CompiledSpectralModel<5> valid;
  REQUIRE(core::compileSpectralModel<5>({ 12.5e-6, 8.5e-6 }, valid)
          == Status::Success);
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(valid)
          == Status::Success);
  const double quiet_nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));

  auto corrupt = valid;
  corrupt.A[0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.B[0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.D[0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.C[0][0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.state_transform[0][0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.x_inner[0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.Cc[0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.cc_coeff = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);

  corrupt = valid;
  corrupt.integration[0][0] = quiet_nan;
  REQUIRE(core::detail::validateCompiledSpectralModelFiniteness(corrupt)
          == Status::Numerical_failure);
}
