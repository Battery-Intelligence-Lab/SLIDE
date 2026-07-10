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
#include <bit>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <limits>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace slide;

static_assert(std::is_same_v<
                decltype(std::declval<core::PackSolver &>().workspace()),
                const core::SolverWorkspace &>,
              "PackSolver must not expose mutable workspace ownership");

namespace {

struct BatchOverlapProbe
{
  void enter()
  {
    std::unique_lock lock{ mutex };
    ++active;
    max_active = std::max(max_active, active);
    ready.notify_all();
    if (!wait_claimed) {
      wait_claimed = true;
      ready.wait_for(lock,
                     std::chrono::milliseconds{ 100 },
                     [&] { return max_active >= 2; });
    }
  }

  void leave()
  {
    const std::lock_guard lock{ mutex };
    --active;
  }

  std::mutex mutex{};
  std::condition_variable ready{};
  int active{};
  int max_active{};
  bool wait_claimed{};
};

struct AffineBatch
{
  std::vector<double> ocv;
  std::vector<double> resistance;
  int calls{};
  BatchOverlapProbe *overlap{};

  Status linearizeThevenin(std::span<const double> current,
                           std::span<double>
                             output_ocv,
                           std::span<double>
                             output_resistance)
  {
    if (current.size() != ocv.size())
      return Status::Invalid_parameters;
    if (overlap != nullptr)
      overlap->enter();
    ++calls;
    std::copy(ocv.begin(), ocv.end(), output_ocv.begin());
    std::copy(resistance.begin(), resistance.end(), output_resistance.begin());
    if (overlap != nullptr)
      overlap->leave();
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
  BatchOverlapProbe overlap;
  AffineBatch a{ .ocv = { 4.0, 4.0 },
                 .resistance = { 0.1, 0.1 },
                 .overlap = &overlap };
  AffineBatch b{ .ocv = { 4.0, 4.0 },
                 .resistance = { 0.1, 0.1 },
                 .overlap = &overlap };
  const std::array<core::TheveninBatchView, 2> batches{
    core::TheveninBatchView::bind(a, 2), core::TheveninBatchView::bind(b, 2)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches, 2) == Status::Success);
  REQUIRE(solver.solve(1.0) == Status::Success);
  REQUIRE(a.calls == solver.diagnostics().iterations);
  REQUIRE(b.calls == solver.diagnostics().iterations);
  REQUIRE(overlap.max_active >= 2);
}

TEST_CASE("parallel Thevenin configuration rejects aliased batch objects",
          "[core][pack][thevenin][alias][P9-B36]")
{
  const auto topology = compile(core::parallel(std::vector{
    core::cell({ .archetype = "a" }),
    core::cell({ .archetype = "b" }) }));
  AffineBatch shared{ .ocv = { 4.0 }, .resistance = { 0.1 } };
  const std::array<core::TheveninBatchView, 2> views{
    core::TheveninBatchView::bind(shared, 1),
    core::TheveninBatchView::bind(shared, 1)
  };
  core::PackSolver solver;
  CHECK(solver.configure(topology, views, 2) == Status::Invalid_parameters);
}

TEST_CASE("production batch execution is bit-repeatable across worker counts",
          "[core][pack][thevenin][determinism][P9-B32]")
{
  std::vector<core::PackNode> cells;
  std::vector<AffineBatch> batches;
  cells.reserve(7);
  batches.reserve(7);
  for (int index = 0; index < 7; ++index) {
    cells.push_back(core::cell(
      { .archetype = "worker-" + std::to_string(index) }));
    batches.push_back(
      { .ocv = { 3.8 + 0.05 * index },
        .resistance = { 0.08 + 0.01 * index } });
  }
  const auto topology = compile(core::parallel(std::move(cells)));
  std::vector<core::TheveninBatchView> views;
  views.reserve(batches.size());
  for (auto &batch : batches)
    views.push_back(core::TheveninBatchView::bind(batch, 1));

  core::PackSolution reference;
  for (const unsigned workers : { 1U, 2U, 7U }) {
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views, workers) == Status::Success);
    REQUIRE(solver.batchWorkerCount() == workers);
    REQUIRE(solver.solve(3.0, core::PackSolveMode::ladder) == Status::Success);
    if (reference.cell_current.empty()) {
      reference = solver.solution();
      continue;
    }
    REQUIRE(solver.solution().cell_current.size()
            == reference.cell_current.size());
    REQUIRE(solver.solution().node_voltage.size()
            == reference.node_voltage.size());
    for (std::size_t index = 0; index < reference.cell_current.size(); ++index)
      CHECK(std::bit_cast<std::uint64_t>(solver.solution().cell_current[index])
            == std::bit_cast<std::uint64_t>(reference.cell_current[index]));
    for (std::size_t index = 0; index < reference.node_voltage.size(); ++index)
      CHECK(std::bit_cast<std::uint64_t>(solver.solution().node_voltage[index])
            == std::bit_cast<std::uint64_t>(reference.node_voltage[index]));
    CHECK(std::bit_cast<std::uint64_t>(solver.solution().terminal_voltage)
          == std::bit_cast<std::uint64_t>(reference.terminal_voltage));
  }
}

TEST_CASE("pack solver rejects invalid modes and malformed compiled netlists atomically",
          "[core][pack][solver][validation][P9-G4]")
{
  const auto topology = compile(core::cell({ .archetype = "affine" }));
  AffineBatch batch{ .ocv = { 4.0 }, .resistance = { 0.1 } };
  const std::array<core::TheveninBatchView, 1> batches{
    core::TheveninBatchView::bind(batch, 1)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(1.0) == Status::Success);
  const auto expected_current = solver.solution().cell_current;
  const auto expected_voltage = solver.solution().node_voltage;
  const auto expected_terminal = solver.solution().terminal_voltage;

  const auto require_unchanged = [&] {
    REQUIRE(solver.solution().cell_current == expected_current);
    REQUIRE(solver.solution().node_voltage == expected_voltage);
    REQUIRE(solver.solution().terminal_voltage == expected_terminal);
  };

  REQUIRE(solver.solve(1.0, static_cast<core::PackSolveMode>(255))
          == Status::Invalid_parameters);
  require_unchanged();

  const auto reject_reconfigure = [&](core::CompiledPackTopology malformed) {
    REQUIRE(solver.configure(malformed, batches) == Status::Invalid_parameters);
    require_unchanged();
    REQUIRE(solver.solve(1.0) == Status::Success);
    require_unchanged();
  };

  auto bad_terminal = topology;
  bad_terminal.electrical.terminal_positive = bad_terminal.electrical.node_count;
  reject_reconfigure(std::move(bad_terminal));

  auto bad_node = topology;
  bad_node.electrical.branches[0].node_positive = bad_node.electrical.node_count;
  reject_reconfigure(std::move(bad_node));

  auto bad_cell = topology;
  bad_cell.electrical.branches[0].cell = static_cast<std::uint32_t>(bad_cell.cells.size());
  reject_reconfigure(std::move(bad_cell));

  auto bad_kind = topology;
  bad_kind.electrical.branches[0].kind = static_cast<core::ElectricalBranchKind>(255);
  reject_reconfigure(std::move(bad_kind));

  auto bad_ladder = topology;
  bad_ladder.electrical.ladder_cells[0] = static_cast<std::uint32_t>(bad_ladder.cells.size());
  reject_reconfigure(std::move(bad_ladder));

  auto impossible_node_count = topology;
  impossible_node_count.electrical.node_count =
    static_cast<std::uint32_t>(std::numeric_limits<int>::max()) + 1U;
  reject_reconfigure(std::move(impossible_node_count));

  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(std::all_of(solver.solution().cell_current.begin(),
                      solver.solution().cell_current.end(),
                      [](double value) { return value == 0.0; }));
  REQUIRE(std::all_of(solver.solution().node_voltage.begin(),
                      solver.solution().node_voltage.end(),
                      [](double value) { return value == 0.0; }));
  REQUIRE(solver.solution().terminal_voltage == 0.0);
}

TEST_CASE("ladder terminal overflow cannot publish a finite-current trial",
          "[core][pack][solver][finite][P9-G4]")
{
  const auto topology = compile(core::series(
    2, core::cell({ .archetype = "affine" })));
  AffineBatch batch{ .ocv = { 4.0, 4.0 }, .resistance = { 1.0, 1.0 } };
  const std::array<core::TheveninBatchView, 1> batches{
    core::TheveninBatchView::bind(batch, 2)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(0.0, core::PackSolveMode::ladder) == Status::Success);
  const auto expected_current = solver.solution().cell_current;
  const auto expected_voltage = solver.solution().node_voltage;
  const auto expected_terminal = solver.solution().terminal_voltage;

  batch.ocv = { 1e308, 1e308 };
  REQUIRE(solver.solve(0.0, core::PackSolveMode::ladder, 1e-12, 2)
          != Status::Success);
  REQUIRE(solver.solution().cell_current == expected_current);
  REQUIRE(solver.solution().node_voltage == expected_voltage);
  REQUIRE(solver.solution().terminal_voltage == expected_terminal);
}

TEST_CASE("all solver modes reject extreme derived values without publishing a trial",
          "[core][pack][solver][finite][P9-G4]")
{
  const auto topology = compile(core::cell({ .archetype = "affine" }));

  const auto exercise = [&](core::PackSolveMode mode,
                            double overflow_ocv,
                            double overflow_resistance,
                            double applied_current) {
    AffineBatch batch{ .ocv = { 4.0 }, .resistance = { 1.0 } };
    const std::array<core::TheveninBatchView, 1> batches{
      core::TheveninBatchView::bind(batch, 1)
    };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, batches) == Status::Success);
    REQUIRE(solver.solve(1.0, mode, 1e-12, 4) == Status::Success);
    const auto expected_current = solver.solution().cell_current;
    const auto expected_voltage = solver.solution().node_voltage;
    const auto expected_terminal = solver.solution().terminal_voltage;

    batch.ocv[0] = overflow_ocv;
    batch.resistance[0] = overflow_resistance;
    REQUIRE(solver.solve(applied_current, mode, 1e-12, 2) != Status::Success);
    REQUIRE(solver.solution().cell_current == expected_current);
    REQUIRE(solver.solution().node_voltage == expected_voltage);
    REQUIRE(solver.solution().terminal_voltage == expected_terminal);
  };

  exercise(core::PackSolveMode::sparse_newton, 0.0, 1e308, 2.0);
  exercise(core::PackSolveMode::ladder, 0.0, 1e308, 2.0);
  exercise(core::PackSolveMode::relaxation, 0.0, 1e308, 2.0);
}

TEST_CASE("pack solve diagnostics are reset when changing modes",
          "[core][pack][solver][diagnostics][P9-G4]")
{
  const auto topology = compile(core::cell({ .archetype = "affine" }));
  AffineBatch batch{ .ocv = { 4.0 }, .resistance = { 1.0 } };
  const std::array<core::TheveninBatchView, 1> batches{
    core::TheveninBatchView::bind(batch, 1)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  REQUIRE(solver.solve(1.0, core::PackSolveMode::sparse_newton, 2.0, 2)
          == Status::Success);
  REQUIRE(solver.diagnostics().residual_norm > 0.0);

  REQUIRE(solver.solve(1.0, core::PackSolveMode::ladder, 1e-12, 2)
          == Status::Success);
  REQUIRE(solver.diagnostics().residual_norm == 0.0);
  REQUIRE(solver.diagnostics().constraint_drift == 0.0);
  REQUIRE(solver.diagnostics().constraint_bound == 0.0);
  REQUIRE(solver.diagnostics().relaxation_gain == 0.0);
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
  constexpr double derivative_step = 1e-4;
  std::array<double, lanes> plus_density{}, minus_density{};
  std::array<double, lanes> plus_voltage{}, minus_voltage{};
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    plus_density[i] = (probe_current[i] + derivative_step) / batch.electrode_area();
    minus_density[i] = (probe_current[i] - derivative_step) / batch.electrode_area();
  }
  REQUIRE(batch.terminalVoltage({ .i_app = plus_density }, plus_voltage) == Status::Success);
  REQUIRE(batch.terminalVoltage({ .i_app = minus_density }, minus_voltage) == Status::Success);
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    const double finite_difference = -(plus_voltage[i] - minus_voltage[i])
                                     / (2.0 * derivative_step);
    REQUIRE(std::abs(resistance[i] - finite_difference)
            <= 1e-7 * std::max(1.0, std::abs(finite_difference)));
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
