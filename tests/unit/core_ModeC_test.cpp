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

TEST_CASE("P4-G2 Baumgarte gain bounds parallel constraint drift",
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
  REQUIRE(solver.setRelaxationGain(0.5) == Status::Success);
  REQUIRE(solver.solve(160.0, core::PackSolveMode::relaxation, 1e-8, 40)
          == Status::Success);
  CAPTURE(solver.diagnostics().iterations,
          solver.diagnostics().constraint_drift,
          solver.diagnostics().constraint_bound);
  REQUIRE(solver.diagnostics().constraint_drift
          <= solver.diagnostics().constraint_bound);
  REQUIRE(solver.diagnostics().relaxation_gain == 0.5);

  auto refused = topology;
  refused.electrical.index1_candidate = false;
  core::PackSolver invalid;
  REQUIRE(invalid.configure(refused, view) == Status::Success);
  REQUIRE(invalid.solve(160.0, core::PackSolveMode::relaxation)
          == Status::Invalid_parameters);
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
  REQUIRE(stepper.step(100'000.0, 0.0, 0.1, {}, core::PackSolveMode::relaxation, 1e-8)
          == Status::Success);

  const auto before = allocation_count.load(std::memory_order_relaxed);
  const auto status = stepper.step(100'000.0, 0.1, 0.1, {}, core::PackSolveMode::relaxation, 1e-8);
  const auto after = allocation_count.load(std::memory_order_relaxed);
  CAPTURE(static_cast<int>(status), stepper.diagnostics().iterations, stepper.diagnostics().constraint_drift);
  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
}
