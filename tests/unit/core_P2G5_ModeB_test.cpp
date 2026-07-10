/**
 * @file core_P2G5_ModeB_test.cpp
 * @brief Registered Mode-B analytical-current admission measurements.
 */

#include "../../src/core/EulerLegacy.hpp"
#include "../../src/core/PackSolver.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace slide;

namespace {

struct AffineBatch
{
  std::vector<double> ocv;
  std::vector<double> resistance;

  Status linearizeThevenin(std::span<const double> current,
                           std::span<double>
                             intercept,
                           std::span<double>
                             tangent_resistance)
  {
    if (current.size() != ocv.size())
      return Status::Invalid_parameters;
    std::copy(ocv.begin(), ocv.end(), intercept.begin());
    std::copy(resistance.begin(), resistance.end(), tangent_resistance.begin());
    return Status::Success;
  }
};

core::CompiledPackTopology parallelTopology(int lanes, std::string archetype)
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(lanes,
                                     core::cell({ .archetype = std::move(archetype) })) },
            topology)
          == Status::Success);
  return topology;
}

void applySpread(core::SpmBatch &batch)
{
  const auto &layout = batch.layout().spm;
  for (int lane = 0; lane < batch.n_lanes(); ++lane) {
    const double coordinate = batch.n_lanes() == 1
                                ? 0.0
                                : 2.0 * lane / static_cast<double>(batch.n_lanes() - 1) - 1.0;
    for (const auto domain : core::domains) {
      const auto d = core::domain_index(domain);
      batch.state().at(layout.active_fraction[d], 0, lane) *= 1.0 + 0.02 * coordinate;
      batch.state().at(layout.specific_surface_area[d], 0, lane) *= 1.0 + 0.02 * coordinate;
      batch.state().at(layout.specific_resistance[d], 0, lane) *= 1.0 + 0.05 * coordinate;
    }
    batch.state().at(layout.current_collector_resistance, 0, lane) *= 1.0 + 0.05 * coordinate;
  }
}

double measureSpmEnvelope(int lanes)
{
  auto input = test_support::make_legacy_kokam_input(0.55, 288.0, 298.0);
  core::SpmBatch sparse_batch, ladder_batch;
  REQUIRE(core::buildSpmBatch(input, {}, lanes, sparse_batch) == Status::Success);
  REQUIRE(core::buildSpmBatch(input, {}, lanes, ladder_batch) == Status::Success);
  applySpread(sparse_batch);
  applySpread(ladder_batch);

  const auto topology = parallelTopology(lanes, "spm");
  const std::array<core::TheveninBatchView, 1> sparse_view{
    core::TheveninBatchView::bind(sparse_batch, lanes)
  };
  const std::array<core::TheveninBatchView, 1> ladder_view{
    core::TheveninBatchView::bind(ladder_batch, lanes)
  };
  core::PackSolver sparse, ladder;
  REQUIRE(sparse.configure(topology, sparse_view) == Status::Success);
  REQUIRE(ladder.configure(topology, ladder_view) == Status::Success);
  core::EulerLegacy sparse_stepper{ sparse_batch };
  core::EulerLegacy ladder_stepper{ ladder_batch };
  std::vector<double> sparse_density(static_cast<std::size_t>(lanes));
  std::vector<double> ladder_density(static_cast<std::size_t>(lanes));

  double max_relative{};
  constexpr int steps = 120;
  for (int step = 0; step < steps; ++step) {
    const double branch_current = step < 60
                                    ? 16.0
                                    : 16.0 - 14.4 * (step - 59) / 60.0;
    const double applied = branch_current * lanes;
    REQUIRE(sparse.solve(applied, core::PackSolveMode::sparse_newton) == Status::Success);
    REQUIRE(ladder.solve(applied, core::PackSolveMode::ladder) == Status::Success);
    for (int lane = 0; lane < lanes; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const double reference = sparse.solution().cell_current[i];
      max_relative = std::max(max_relative,
                              std::abs(ladder.solution().cell_current[i] - reference)
                                / std::max(1e-12, std::abs(reference)));
      sparse_density[i] = reference / sparse_batch.electrode_area();
      ladder_density[i] = ladder.solution().cell_current[i]
                          / ladder_batch.electrode_area();
    }
    REQUIRE(sparse_stepper.step(sparse_batch, sparse_density, step, 1.0)
            == Status::Success);
    REQUIRE(ladder_stepper.step(ladder_batch, ladder_density, step, 1.0)
            == Status::Success);
  }
  return max_relative;
}

} // namespace

TEST_CASE("P2-G5a Mode B is exact for a linear 16p resistance spread",
          "[core][pack][mode-b][P2-G5]")
{
  constexpr int lanes = 16;
  AffineBatch sparse_batch, ladder_batch;
  for (int lane = 0; lane < lanes; ++lane) {
    sparse_batch.ocv.push_back(4.0);
    sparse_batch.resistance.push_back(0.1 * (0.995 + 0.01 * lane / 15.0));
  }
  ladder_batch = sparse_batch;
  const auto topology = parallelTopology(lanes, "linear");
  const std::array<core::TheveninBatchView, 1> sparse_view{
    core::TheveninBatchView::bind(sparse_batch, lanes)
  };
  const std::array<core::TheveninBatchView, 1> ladder_view{
    core::TheveninBatchView::bind(ladder_batch, lanes)
  };
  core::PackSolver sparse, ladder;
  REQUIRE(sparse.configure(topology, sparse_view) == Status::Success);
  REQUIRE(ladder.configure(topology, ladder_view) == Status::Success);
  REQUIRE(sparse.solve(160.0) == Status::Success);
  REQUIRE(ladder.solve(160.0, core::PackSolveMode::ladder) == Status::Success);
  double maximum{};
  for (int lane = 0; lane < lanes; ++lane)
    maximum = std::max(maximum,
                       std::abs(sparse.solution().cell_current[static_cast<std::size_t>(lane)]
                                - ladder.solution().cell_current[static_cast<std::size_t>(lane)]));
  REQUIRE(maximum <= 1e-8);
}

TEST_CASE("P2-G5 SPM Mode B envelope is admitted at 16p and 256p",
          "[core][pack][mode-b][spm][P2-G5]")
{
  const double error_16p = measureSpmEnvelope(16);
  const double error_256p = measureSpmEnvelope(256);
  std::printf("P2-G5 SPM envelope: 16p=%.17g 256p=%.17g\n",
              error_16p,
              error_256p);
  CAPTURE(error_16p, error_256p);
  REQUIRE(error_16p <= 1e-3);
  REQUIRE(error_256p <= 1e-3);
}
