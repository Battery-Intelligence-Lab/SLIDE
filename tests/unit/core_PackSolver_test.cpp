/**
 * @file core_PackSolver_test.cpp
 * @brief Phase-2 affine Thevenin, Mode A/B, and workspace-memory gates.
 */

#include "../../src/core/PackSolver.hpp"
#include "../../src/core/EulerLegacy.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

using namespace slide;

namespace {

struct AffineBatch
{
  std::vector<double> ocv;
  std::vector<double> resistance;
  int calls{};

  Status linearizeThevenin(std::span<const double> current,
                           std::span<double>
                             output_ocv,
                           std::span<double>
                             output_resistance)
  {
    if (current.size() != ocv.size())
      return Status::Invalid_parameters;
    ++calls;
    std::copy(ocv.begin(), ocv.end(), output_ocv.begin());
    std::copy(resistance.begin(), resistance.end(), output_resistance.begin());
    return Status::Success;
  }
};

core::CompiledPackTopology compile(const core::PackNode &root)
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription({ .root = root }, topology) == Status::Success);
  return topology;
}

} // namespace

TEST_CASE("Mode A solves heterogeneous affine parallel cells and reuses factorization",
          "[core][pack][mode-a][workspace]")
{
  const auto topology = compile(core::parallel(4, core::cell({ .archetype = "affine" })));
  AffineBatch batch{ .ocv = { 4.0, 4.1, 3.9, 4.2 },
                     .resistance = { 0.1, 0.2, 0.15, 0.3 } };
  const std::array<core::TheveninBatchView, 1> batches{ core::TheveninBatchView::bind(batch, 4) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(2.0) == Status::Success);
  REQUIRE(solver.diagnostics().iterations <= 2);
  REQUIRE(solver.workspace().symbolicFactorizations() == 1);
  REQUIRE(solver.workspace().numericFactorizations() == 1);

  double current_sum{};
  for (std::size_t cell = 0; cell < batch.ocv.size(); ++cell) {
    const double current = solver.solution().cell_current[cell];
    current_sum += current;
    REQUIRE(std::abs(batch.ocv[cell] - batch.resistance[cell] * current
                     - solver.solution().terminal_voltage)
            <= 1e-12);
  }
  REQUIRE(std::abs(current_sum - 2.0) <= 1e-12);

  const auto first_current = solver.solution().cell_current;
  REQUIRE(solver.solve(2.0) == Status::Success);
  REQUIRE(solver.diagnostics().iterations == 1);
  REQUIRE(solver.workspace().numericFactorizations() == 1);
  REQUIRE(solver.solution().cell_current == first_current);

  solver.invalidate();
  REQUIRE(solver.solve(2.0) == Status::Success);
  REQUIRE(solver.workspace().numericFactorizations() == 2);
  REQUIRE(solver.solution().cell_current == first_current);

  const int factorization_before_segment = solver.workspace().numericFactorizations();
  for (int step = 0; step < 100; ++step)
    REQUIRE(solver.solve(2.0) == Status::Success);
  REQUIRE(solver.workspace().numericFactorizations() - factorization_before_segment <= 10);
  REQUIRE(solver.workspace().numericFactorizations() == factorization_before_segment);
}

TEST_CASE("Mode B is equivalent to Mode A on a series-of-parallel ladder",
          "[core][pack][mode-a][mode-b][P2-G2]")
{
  const auto topology = compile(core::series(
    3, core::parallel(2, core::cell({ .archetype = "affine" }))));
  REQUIRE(topology.electrical.series_parallel_ladder);
  AffineBatch sparse_batch{ .ocv = { 4.0, 4.1, 3.9, 4.2, 4.05, 3.95 },
                            .resistance = { 0.10, 0.12, 0.09, 0.15, 0.11, 0.13 } };
  AffineBatch ladder_batch = sparse_batch;
  const std::array<core::TheveninBatchView, 1> sparse_view{
    core::TheveninBatchView::bind(sparse_batch, 6)
  };
  const std::array<core::TheveninBatchView, 1> ladder_view{
    core::TheveninBatchView::bind(ladder_batch, 6)
  };
  core::PackSolver sparse, ladder;
  REQUIRE(sparse.configure(topology, sparse_view) == Status::Success);
  REQUIRE(ladder.configure(topology, ladder_view) == Status::Success);
  REQUIRE(sparse.solve(3.0, core::PackSolveMode::sparse_newton) == Status::Success);
  REQUIRE(ladder.solve(3.0, core::PackSolveMode::ladder) == Status::Success);
  REQUIRE(ladder.workspace().numericFactorizations() == 0);
  REQUIRE(std::abs(sparse.solution().terminal_voltage - ladder.solution().terminal_voltage)
          <= 1e-12);
  for (std::size_t cell = 0; cell < sparse.solution().cell_current.size(); ++cell)
    REQUIRE(std::abs(sparse.solution().cell_current[cell]
                     - ladder.solution().cell_current[cell])
            <= 1e-12);
}

TEST_CASE("nested topology solves in one global loop and workspace amortizes factorization",
          "[core][pack][P2-G3][P2-G4]")
{
  const auto topology = compile(core::series(
    2, core::parallel(2, core::parallel(2, core::cell({ .archetype = "affine" })))));
  AffineBatch batch{ .ocv = { 4.0, 4.02, 3.98, 4.04, 3.95, 4.01, 4.03, 3.99 },
                     .resistance = { 0.08, 0.09, 0.11, 0.13, 0.07, 0.12, 0.10, 0.15 } };
  const std::array<core::TheveninBatchView, 1> batches{ core::TheveninBatchView::bind(batch, 8) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  for (int step = 0; step < 100; ++step) {
    REQUIRE(solver.solve(4.0) == Status::Success);
    REQUIRE(solver.diagnostics().iterations <= 8);
  }
  REQUIRE(solver.workspace().numericFactorizations() <= 10);
  REQUIRE(solver.workspace().numericFactorizations() == 1);
  REQUIRE(batch.calls <= 101); // first solve needs two affine confirmations; warm solves need one
}

TEST_CASE("Thevenin system dispatches once per archetype batch", "[core][pack][thevenin]")
{
  const auto root = core::parallel(std::vector{
    core::cell({ .archetype = "b" }), core::cell({ .archetype = "a" }), core::cell({ .archetype = "b" }), core::cell({ .archetype = "a" }) });
  const auto topology = compile(root);
  REQUIRE(topology.batch_archetypes == std::vector<std::string>{ "a", "b" });
  AffineBatch a{ .ocv = { 4.0, 4.0 }, .resistance = { 0.1, 0.1 } };
  AffineBatch b{ .ocv = { 4.0, 4.0 }, .resistance = { 0.1, 0.1 } };
  const std::array<core::TheveninBatchView, 2> batches{
    core::TheveninBatchView::bind(a, 2), core::TheveninBatchView::bind(b, 2)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(1.0) == Status::Success);
  REQUIRE(a.calls == solver.diagnostics().iterations);
  REQUIRE(b.calls == solver.diagnostics().iterations);
}

TEST_CASE("SPM batches expose a nonlinear Thevenin tangent to the pack solver",
          "[core][pack][spm][chord]")
{
  constexpr int lanes = 4;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, { .nch = 5 }, lanes, batch) == Status::Success);
  auto collector = batch.state().row(batch.layout().spm.current_collector_resistance.row_begin);
  collector[0] *= 0.8;
  collector[1] *= 1.0;
  collector[2] *= 1.2;
  collector[3] *= 1.4;

  const std::array<double, lanes> probe_current{ 8.0, 12.0, 16.0, 20.0 };
  std::array<double, lanes> intercept{}, resistance{}, voltage{};
  REQUIRE(batch.linearizeThevenin(probe_current, intercept, resistance)
          == Status::Success);
  std::array<double, lanes> density{};
  for (int lane = 0; lane < lanes; ++lane)
    density[static_cast<std::size_t>(lane)] = probe_current[static_cast<std::size_t>(lane)]
                                              / batch.electrode_area();
  REQUIRE(batch.terminalVoltage({ .time = 0.0, .dt = 0.0, .i_app = density }, voltage)
          == Status::Success);
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    REQUIRE(resistance[i] > 0.0);
    REQUIRE(std::abs(intercept[i] - resistance[i] * probe_current[i] - voltage[i])
            <= 1e-14);
  }

  const auto topology = compile(core::parallel(
    lanes, core::cell({ .archetype = "spm" })));
  const std::array<core::TheveninBatchView, 1> batches{
    core::TheveninBatchView::bind(batch, lanes)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(64.0, core::PackSolveMode::sparse_newton, 1e-10, 8)
          == Status::Success);
  REQUIRE(solver.diagnostics().iterations <= 8);

  double sum{};
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    density[i] = solver.solution().cell_current[i] / batch.electrode_area();
    sum += solver.solution().cell_current[i];
  }
  REQUIRE(batch.terminalVoltage({ .time = 0.0, .dt = 0.0, .i_app = density }, voltage)
          == Status::Success);
  REQUIRE(std::abs(sum - 64.0) <= 1e-10);
  for (const auto cell_voltage : voltage)
    REQUIRE(std::abs(cell_voltage - solver.solution().terminal_voltage) <= 1e-10);

  const int factorization_before = solver.workspace().numericFactorizations();
  core::EulerLegacy stepper{ batch };
  for (int step = 0; step < 100; ++step) {
    REQUIRE(solver.solve(64.0) == Status::Success);
    for (int lane = 0; lane < lanes; ++lane)
      density[static_cast<std::size_t>(lane)] = solver.solution().cell_current[static_cast<std::size_t>(lane)]
                                                / batch.electrode_area();
    REQUIRE(stepper.step(batch, density, static_cast<double>(step), 1.0)
            == Status::Success);
  }
  REQUIRE(solver.workspace().numericFactorizations() - factorization_before <= 10);
}
