/**
 * @file core_ModeC_test.cpp
 * @brief Phase-4 waveform-relaxation/Baumgarte and 100k scale gates.
 */

#include "../../src/core/PackSolver.hpp"
#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <new>
#include <vector>

static std::atomic<std::size_t> allocation_count{ 0 };

#if defined(_WIN32)
#include <malloc.h>
static void *aligned_allocate(std::size_t bytes, std::size_t alignment)
{
  return _aligned_malloc(bytes, alignment);
}
static void aligned_release(void *pointer) { _aligned_free(pointer); }
#else
static void *aligned_allocate(std::size_t bytes, std::size_t alignment)
{
  return std::aligned_alloc(alignment,
                            ((bytes + alignment - 1) / alignment) * alignment);
}
static void aligned_release(void *pointer) { std::free(pointer); }
#endif

void *operator new(std::size_t bytes)
{
  ++allocation_count;
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes)
{
  ++allocation_count;
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  ++allocation_count;
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  ++allocation_count;
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete[](void *pointer, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept { aligned_release(pointer); }

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
                             output_resistance)
  {
    if (current.size() != ocv.size())
      return Status::Invalid_parameters;
    std::copy(ocv.begin(), ocv.end(), intercept.begin());
    std::copy(resistance.begin(), resistance.end(), output_resistance.begin());
    return Status::Success;
  }
};

core::CompiledPackTopology parallelTopology(int lanes)
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(lanes,
                                     core::cell({ .archetype = "affine" })) },
            topology)
          == Status::Success);
  return topology;
}

core::CompiledPackTopology compileTopology(const core::PackNode &root)
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription({ .root = root }, topology)
          == Status::Success);
  return topology;
}

AffineBatch heterogeneous(int lanes)
{
  AffineBatch batch;
  batch.ocv.reserve(static_cast<std::size_t>(lanes));
  batch.resistance.reserve(static_cast<std::size_t>(lanes));
  for (int lane = 0; lane < lanes; ++lane) {
    const double coordinate = lanes == 1
                                ? 0.0
                                : 2.0 * lane / static_cast<double>(lanes - 1) - 1.0;
    batch.ocv.push_back(4.0 + 0.05 * coordinate);
    batch.resistance.push_back(0.1 * (1.0 + 0.1 * coordinate));
  }
  return batch;
}

} // namespace

TEST_CASE("P4-G1 Mode C agrees with Mode A on heterogeneous 16p",
          "[core][pack][mode-c][P4-G1]")
{
  constexpr int lanes = 16;
  const auto topology = parallelTopology(lanes);
  auto exact_batch = heterogeneous(lanes);
  auto relaxation_batch = heterogeneous(lanes);
  const std::array<core::TheveninBatchView, 1> exact_view{
    core::TheveninBatchView::bind(exact_batch, lanes)
  };
  const std::array<core::TheveninBatchView, 1> relaxation_view{
    core::TheveninBatchView::bind(relaxation_batch, lanes)
  };
  core::PackSolver exact, relaxation;
  REQUIRE(exact.configure(topology, exact_view) == Status::Success);
  REQUIRE(relaxation.configure(topology, relaxation_view) == Status::Success);
  REQUIRE(exact.solve(160.0) == Status::Success);
  REQUIRE(relaxation.solve(160.0, core::PackSolveMode::relaxation)
          == Status::Success);
  double maximum_relative{};
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    maximum_relative = std::max(
      maximum_relative,
      std::abs(relaxation.solution().cell_current[i]
               - exact.solution().cell_current[i])
        / std::max(1e-12, std::abs(exact.solution().cell_current[i])));
  }
  REQUIRE(maximum_relative <= 1e-3);
  REQUIRE(relaxation.workspace().numericFactorizations() == 0);
}

TEST_CASE("P4-G2 one-node relaxation follows its analytic contraction",
          "[core][pack][mode-c][baumgarte][P4-G2]")
{
  constexpr int lanes = 16;
  const auto topology = parallelTopology(lanes);
  auto batch = heterogeneous(lanes);
  const std::array<core::TheveninBatchView, 1> view{
    core::TheveninBatchView::bind(batch, lanes)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, view) == Status::Success);
  constexpr double applied_current = 160.0;
  constexpr double alpha = 0.5;
  double initial_kcl = applied_current;
  double operation_scale = std::abs(applied_current);
  for (std::size_t lane = 0; lane < batch.ocv.size(); ++lane) {
    initial_kcl -= batch.ocv[lane] / batch.resistance[lane];
    operation_scale += std::abs(batch.ocv[lane] / batch.resistance[lane]);
  }
  initial_kcl = std::abs(initial_kcl);

  REQUIRE(solver.setRelaxationGain(alpha) == Status::Success);
  REQUIRE(solver.solve(applied_current,
                       core::PackSolveMode::relaxation,
                       1e-8,
                       40)
          == Status::Success);
  CAPTURE(solver.diagnostics().iterations,
          solver.diagnostics().constraint_drift,
          solver.diagnostics().constraint_bound);
  const double analytic_drift = initial_kcl
                                * std::pow(1.0 - alpha,
                                           solver.diagnostics().iterations);
  const double arithmetic_allowance = 512.0
                                      * std::numeric_limits<double>::epsilon()
                                      * operation_scale
                                      * solver.diagnostics().iterations;
  REQUIRE(solver.diagnostics().constraint_drift
          <= analytic_drift + arithmetic_allowance);
  REQUIRE(solver.diagnostics().constraint_drift <= 1e-8);
  REQUIRE(solver.diagnostics().residual_norm <= 1e-8);
  REQUIRE(solver.diagnostics().relaxation_gain == alpha);

  auto refused = topology;
  refused.electrical.index1_candidate = false;
  core::PackSolver invalid;
  REQUIRE(invalid.configure(refused, view) == Status::Success);
  REQUIRE(invalid.solve(160.0, core::PackSolveMode::relaxation)
          == Status::Invalid_parameters);
}

TEST_CASE("Mode C checks internal-node KCL, not only the terminal",
          "[core][pack][mode-c][residual][P9-G4]")
{
  constexpr int parallel_cells = 3;
  const auto topology = compileTopology(core::series(std::vector{
    core::cell({ .archetype = "affine" }),
    core::parallel(
      parallel_cells, core::cell({ .archetype = "affine" })) }));
  AffineBatch batch{ .ocv = std::vector<double>(1 + parallel_cells, 4.0),
                     .resistance = std::vector<double>(1 + parallel_cells, 1.0) };
  const std::array<core::TheveninBatchView, 1> view{
    core::TheveninBatchView::bind(batch, 1 + parallel_cells)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, view) == Status::Success);
  REQUIRE(solver.solve(0.0, core::PackSolveMode::ladder) == Status::Success);
  const auto expected_current = solver.solution().cell_current;
  const auto expected_voltage = solver.solution().node_voltage;
  const auto expected_terminal = solver.solution().terminal_voltage;

  constexpr double branch_change = 1e-3;
  batch.ocv[0] += parallel_cells * branch_change;
  for (int lane = 1; lane <= parallel_cells; ++lane)
    batch.ocv[static_cast<std::size_t>(lane)] += branch_change;
  REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
  REQUIRE(solver.solve(0.0,
                       core::PackSolveMode::relaxation,
                       1.5 * branch_change,
                       1)
          != Status::Success);
  REQUIRE(solver.diagnostics().residual_norm <= 1.5 * branch_change);
  REQUIRE(solver.diagnostics().constraint_drift
          > 1.5 * branch_change);
  REQUIRE(solver.solution().cell_current == expected_current);
  REQUIRE(solver.solution().node_voltage == expected_voltage);
  REQUIRE(solver.solution().terminal_voltage == expected_terminal);
}

TEST_CASE("Mode C cannot converge on update stagnation while KCL is violated",
          "[core][pack][mode-c][residual][P9-G4]")
{
  const auto topology = parallelTopology(1);
  AffineBatch batch{ .ocv = { 4.0 }, .resistance = { 1.0 } };
  const std::array<core::TheveninBatchView, 1> view{
    core::TheveninBatchView::bind(batch, 1)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, view) == Status::Success);
  REQUIRE(solver.solve(0.0, core::PackSolveMode::ladder) == Status::Success);
  const auto expected_current = solver.solution().cell_current;
  const auto expected_voltage = solver.solution().node_voltage;
  const auto expected_terminal = solver.solution().terminal_voltage;

  REQUIRE(solver.setRelaxationGain(1e-20) == Status::Success);
  REQUIRE(solver.solve(1.0, core::PackSolveMode::relaxation, 1e-12, 2)
          != Status::Success);
  REQUIRE(solver.diagnostics().constraint_drift > 1e-6);
  REQUIRE(solver.diagnostics().constraint_bound == 1e-12);
  REQUIRE(solver.solution().cell_current == expected_current);
  REQUIRE(solver.solution().node_voltage == expected_voltage);
  REQUIRE(solver.solution().terminal_voltage == expected_terminal);
}

TEST_CASE("P4-G3 100k-cell Mode C is sub-GB and allocation-free",
          "[core][pack][mode-c][scale][P4-G3]")
{
  constexpr int lanes = 100'000;
  const auto topology = parallelTopology(lanes);
  AffineBatch batch;
  batch.ocv.assign(lanes, 4.0);
  batch.resistance.assign(lanes, 0.1);
  const std::array<core::TheveninBatchView, 1> view{
    core::TheveninBatchView::bind(batch, lanes)
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, view) == Status::Success);
  REQUIRE(solver.solve(100'000.0, core::PackSolveMode::relaxation)
          == Status::Success);
  constexpr std::size_t state_bytes = static_cast<std::size_t>(lanes) * 240;
  REQUIRE(state_bytes < static_cast<std::size_t>(1'000'000'000));

  const auto before = allocation_count.load(std::memory_order_relaxed);
  const auto status = solver.solve(100'000.0, core::PackSolveMode::relaxation);
  const auto after = allocation_count.load(std::memory_order_relaxed);
  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
  REQUIRE(solver.workspace().numericFactorizations() == 0);
}

TEST_CASE("PAY-3 100k SPM cells complete an accepted pack advance",
          "[core][pack][mode-c][scale][PAY-3]")
{
  constexpr int lanes = 100'000;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, lanes, batch) == Status::Success);
  const auto topology = parallelTopology(lanes);
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Success);
  REQUIRE(batch.state().size() * sizeof(double)
          < static_cast<std::size_t>(1'000'000'000));
  const auto first_status = stepper.step(
    100'000.0, 0.0, 0.1, {}, core::PackSolveMode::relaxation, 1e-8);
  CAPTURE(static_cast<int>(first_status),
          stepper.diagnostics().iterations,
          stepper.diagnostics().residual_norm,
          stepper.diagnostics().constraint_drift,
          stepper.diagnostics().constraint_bound);
  REQUIRE(first_status == Status::Success);
  REQUIRE(stepper.diagnostics().constraint_drift <= 1e-8);

  const auto before = allocation_count.load(std::memory_order_relaxed);
  const auto status = stepper.step(100'000.0, 0.1, 0.1, {}, core::PackSolveMode::relaxation, 1e-8);
  const auto after = allocation_count.load(std::memory_order_relaxed);
  CAPTURE(static_cast<int>(status), stepper.diagnostics().iterations, stepper.diagnostics().constraint_drift);
  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
}
