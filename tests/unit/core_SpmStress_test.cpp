/**
 * @file core_SpmStress_test.cpp
 * @brief Direct legacy parity for Dai and Laresgoiti stress observables.
 *
 * REGISTERED BEFORE FIRST RUN:
 *  - Both electrode Dai maxima and negative Laresgoiti stress match the public
 *    Cell_SPM stress functions to relative error <=1e-12.
 *  - For four lanes with distinct constant concentrations, every node in both
 *    domains is uniform and all eight normalized Dai maxima satisfy
 *    |sigma_h| / (Omega*E/(1-nu)*c0) <= 1e-10.
 *  - Adding a lane-distinct constant to every node of a non-uniform profile
 *    changes each normalized maximum by <=1e-10, using
 *    Omega*E/(1-nu)*max_node(c) as the stress scale.
 *    Distinct lane values expose node-major stride/lane-index mixing.
 *
 * LIMITATIONS: both analytic invariants are blind to the material factor's
 * magnitude and sign and to the deliberate legacy node association. Those
 * remain covered only by the legacy replay; a closed-form non-uniform Dai
 * solution would be needed for an independent oracle.
 *
 * OBSERVED AFTER FIRST RUN (2026-07-23): maximum normalized residuals were
 * 1.4e-16 (uniform) and 5e-17 (offset) in Debug, and 1.3e-16 (uniform) and
 * 6e-17 (offset) in fast-math Release.
 */

#include "../../src/slide.hpp"
#include "../../src/core/SpmStress.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

using namespace slide;

namespace {

class StressReferenceCell : public Cell_SPM
{
public:
  template <int NCH>
  core::SpmStressParams<NCH> params(const Model_SPM<NCH> &model) const
  {
    core::SpmStressParams<NCH> p;
    for (int i = 0; i < NCH; ++i)
      p.x_inner[static_cast<std::size_t>(i)] = model.xch(i);
    for (int row = 0; row < 2 * NCH + 3; ++row)
      for (int column = 0; column < 2 * NCH + 3; ++column)
        p.integration[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] = model.Q(row, column);
    core::domain_value(p.partial_molar_volume, core::Domain::pos) = sparam.omegap;
    core::domain_value(p.partial_molar_volume, core::Domain::neg) = sparam.omegan;
    core::domain_value(p.youngs_modulus, core::Domain::pos) = sparam.Ep;
    core::domain_value(p.youngs_modulus, core::Domain::neg) = sparam.En;
    core::domain_value(p.poisson_ratio, core::Domain::pos) = sparam.nup;
    core::domain_value(p.poisson_ratio, core::Domain::neg) = sparam.nun;
    constexpr std::array x{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    constexpr std::array y{ 0.0, 5.5, 8.5, 9.5, 10.0, 10.0, 10.5, 13.0, 16.5, 21.0, 23.5 };
    REQUIRE(p.laresgoiti_negative.build(x, y) == Status::Success);
    return p;
  }
};

constexpr std::size_t core_index(slide::Domain domain)
{
  return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
}

} // namespace

TEST_CASE("SPM stress observables match legacy Cell_SPM", "[core][ageing][stress]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  StressReferenceCell cell;
  cell.setT(307.0);
  REQUIRE(cell.setCurrent(18.0, false, false) == Status::Success);
  auto *model = Model_SPM<>::makeModel();
  const auto params = cell.params(*model);
  REQUIRE(core::validateSpmStressParams(params) == Status::Success);

  double dai_positive{}, dai_negative{}, lares_negative{};
  Cell_SPM::sigma_type radial_positive{}, radial_negative{}, tangential_positive{},
    tangential_negative{}, hydrostatic_positive{}, hydrostatic_negative{};
  cell.getDaiStress(&dai_positive, &dai_negative, radial_positive, radial_negative, tangential_positive, tangential_negative, hydrostatic_positive, hydrostatic_negative);
  cell.getLaresgoitiStress(false, &lares_negative);

  core::SpmObservableScratch<NCH> observable_scratch{ 1 };
  auto observables = observable_scratch.view();
  for (const auto legacy_domain : { pos, neg }) {
    const auto d = core_index(legacy_domain);
    const auto concentration = cell.getC(legacy_domain);
    for (int node = 0; node < NCH + 2; ++node)
      observables.concentration[d][static_cast<std::size_t>(node)] = concentration[static_cast<std::size_t>(node)];
  }
  DPair surface{};
  REQUIRE(cell.getCSurf(surface, false));
  observables.surface_stoichiometry[core_index(pos)][0] = surface[pos] / 51385.0;
  observables.surface_stoichiometry[core_index(neg)][0] = surface[neg] / 30555.0;

  core::SpmStressScratch scratch{ 1 };
  const auto output = scratch.view();
  REQUIRE(core::computeSpmStress(params, observables, 1, output) == Status::Success);

  const std::array expected{ dai_positive, dai_negative };
  const std::array actual{ output.dai_maximum_hydrostatic[core_index(pos)][0],
                           output.dai_maximum_hydrostatic[core_index(neg)][0] };
  for (int domain = 0; domain < 2; ++domain) {
    const double scale = std::max(std::abs(expected[static_cast<std::size_t>(domain)]), 1e-30);
    CAPTURE(domain, expected[static_cast<std::size_t>(domain)], actual[static_cast<std::size_t>(domain)]);
    REQUIRE(std::abs(actual[static_cast<std::size_t>(domain)]
                     - expected[static_cast<std::size_t>(domain)])
              / scale
            <= 1e-12);
  }
  REQUIRE(std::abs(output.laresgoiti_negative[0] - lares_negative)
            / std::max(std::abs(lares_negative), 1e-30)
          <= 1e-12);

  auto explosive = params;
  explosive.partial_molar_volume[core::domain_index(core::Domain::neg)] =
    std::numeric_limits<double>::max();
  explosive.youngs_modulus[core::domain_index(core::Domain::neg)] =
    std::numeric_limits<double>::max();
  REQUIRE(core::computeSpmStress(explosive, observables, 1, output)
          == Status::Numerical_failure);
}

TEST_CASE("Dai hydrostatic stress obeys uniform and constant-offset invariants",
          "[core][ageing][stress][oracle]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr int lanes = 4;
  constexpr std::array concentrations{ 5.0e3, 1.5e4, 2.5e4, 3.0e4 };
  constexpr std::array offsets{ 1.0e3, 2.0e3, 3.0e3, 4.0e3 };

  StressReferenceCell cell;
  auto *model = Model_SPM<>::makeModel();
  const auto params = cell.params(*model);
  REQUIRE(core::validateSpmStressParams(params) == Status::Success);

  core::SpmObservableScratch<NCH> observable_scratch{ lanes };
  auto observables = observable_scratch.view();
  core::SpmStressScratch output_scratch{ lanes };
  const auto output = output_scratch.view();

  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int node = 0; node < NCH + 2; ++node)
      for (int lane = 0; lane < lanes; ++lane)
        observables.concentration[d][static_cast<std::size_t>(node * lanes + lane)] =
          concentrations[static_cast<std::size_t>(lane)];
  }

  REQUIRE(core::computeSpmStress(params, observables, lanes, output)
          == Status::Success);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    const double factor = params.partial_molar_volume[d] * params.youngs_modulus[d]
                          / (1.0 - params.poisson_ratio[d]);
    for (int lane = 0; lane < lanes; ++lane) {
      const double stress = output.dai_maximum_hydrostatic[d][static_cast<std::size_t>(lane)];
      const double scale = factor * concentrations[static_cast<std::size_t>(lane)];
      const double normalized = std::abs(stress) / scale;
      CAPTURE(d, lane, stress, scale, normalized);
      REQUIRE(normalized <= 1e-10);
    }
  }

  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int node = 0; node < NCH + 2; ++node)
      for (int lane = 0; lane < lanes; ++lane)
        observables.concentration[d][static_cast<std::size_t>(node * lanes + lane)] =
          concentrations[static_cast<std::size_t>(lane)] * (1.0 + 0.1 * node);
  }
  REQUIRE(core::computeSpmStress(params, observables, lanes, output)
          == Status::Success);
  core::PerDomain<std::array<double, lanes>> baseline{};
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane)
      baseline[d][static_cast<std::size_t>(lane)] =
        output.dai_maximum_hydrostatic[d][static_cast<std::size_t>(lane)];
  }

  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int node = 0; node < NCH + 2; ++node)
      for (int lane = 0; lane < lanes; ++lane)
        observables.concentration[d][static_cast<std::size_t>(node * lanes + lane)] =
          concentrations[static_cast<std::size_t>(lane)] * (1.0 + 0.1 * node)
          + offsets[static_cast<std::size_t>(lane)];
  }
  REQUIRE(core::computeSpmStress(params, observables, lanes, output)
          == Status::Success);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    const double factor = params.partial_molar_volume[d] * params.youngs_modulus[d]
                          / (1.0 - params.poisson_ratio[d]);
    for (int lane = 0; lane < lanes; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const double stress = output.dai_maximum_hydrostatic[d][i];
      const double difference = std::abs(stress - baseline[d][i]);
      const double max_node = concentrations[i] * (1.0 + 0.1 * (NCH + 1));
      const double scale = factor * max_node;
      const double normalized = difference / scale;
      CAPTURE(d, lane, baseline[d][i], stress, difference, scale, normalized);
      REQUIRE(normalized <= 1e-10);
    }
  }
}

TEST_CASE("SPM stress parameters reject missing curves and invalid mechanics",
          "[core][ageing][stress][validation]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  core::SpmStressParams<NCH> invalid;
  REQUIRE(core::validateSpmStressParams(invalid)
          == Status::Invalid_parameters);

  StressReferenceCell cell;
  auto *model = Model_SPM<>::makeModel();
  auto params = cell.params(*model);
  params.youngs_modulus[core::domain_index(core::Domain::neg)] = 0.0;
  REQUIRE(core::validateSpmStressParams(params)
          == Status::Invalid_parameters);
}

TEST_CASE("Stress history is explicit checkpointed state", "[core][ageing][stress]")
{
  core::BatchBuilder builder;
  const auto layout = core::declareStressHistory(builder);
  auto arena = builder.build(2);
  const auto neg = core::domain_index(core::Domain::neg);
  const auto pos = core::domain_index(core::Domain::pos);
  arena.at(layout.previous_dai[neg], 0, 0) = -12.0;
  arena.at(layout.previous_dai[pos], 0, 0) = 8.0;
  arena.at(layout.previous_laresgoiti_negative, 0, 0) = 4.0;
  arena.at(layout.interval, 0, 0) = 60.0;
  std::vector<double> snapshot(arena.raw().size());
  std::copy(arena.raw().begin(), arena.raw().end(), snapshot.begin());
  std::fill(arena.raw().begin(), arena.raw().end(), 0.0);
  std::copy(snapshot.begin(), snapshot.end(), arena.raw().begin());
  REQUIRE(arena.at(layout.previous_dai[neg], 0, 0) == -12.0);
  REQUIRE(arena.at(layout.previous_dai[pos], 0, 0) == 8.0);
  REQUIRE(arena.at(layout.previous_laresgoiti_negative, 0, 0) == 4.0);
  REQUIRE(arena.at(layout.interval, 0, 0) == 60.0);
  REQUIRE_FALSE(builder.is_ode_row(layout.previous_dai[neg].row_begin));
  REQUIRE_FALSE(builder.is_ode_row(layout.previous_laresgoiti_negative.row_begin));
  REQUIRE_FALSE(builder.is_ode_row(layout.interval.row_begin));
}
