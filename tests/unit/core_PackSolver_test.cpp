/**
 * @file core_PackSolver_test.cpp
 * @brief Phase-2 affine Thevenin, Mode A/B, and workspace-memory gates.
 */

#include "../../src/core/PackSolver.hpp"
#include "../../src/core/PackSolverInternal.hpp"
#include "../../src/core/EulerLegacy.hpp"
#include "../support/KokamSpmFixture.hpp"
#include "../support/RecordedBits.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64) \
  || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#include <xmmintrin.h>
#define SLIDE_TEST_HAS_X86_FTZ 1
#endif

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

struct ScriptedAffineBatch
{
  std::vector<std::vector<double>> ocv_by_call;
  std::vector<std::vector<double>> resistance_by_call;
  std::size_t calls{};

  Status linearizeThevenin(std::span<const double> current,
                           std::span<double>
                             output_ocv,
                           std::span<double>
                             output_resistance)
  {
    if (ocv_by_call.empty() || resistance_by_call.empty()
        || current.size() != ocv_by_call.front().size()
        || output_ocv.size() != current.size()
        || output_resistance.size() != current.size())
      return Status::Invalid_parameters;
    const auto index = std::min(calls, ocv_by_call.size() - 1);
    const auto resistance_index =
      std::min(calls, resistance_by_call.size() - 1);
    ++calls;
    std::copy(ocv_by_call[index].begin(),
              ocv_by_call[index].end(),
              output_ocv.begin());
    std::copy(resistance_by_call[resistance_index].begin(),
              resistance_by_call[resistance_index].end(),
              output_resistance.begin());
    return Status::Success;
  }
};

struct TracingAffineBatch
{
  std::vector<double> ocv;
  std::vector<double> resistance;
  std::vector<std::vector<double>> callback_current;

  Status linearizeThevenin(std::span<const double> current,
                           std::span<double>
                             output_ocv,
                           std::span<double>
                             output_resistance)
  {
    if (current.size() != ocv.size() || output_ocv.size() != ocv.size()
        || output_resistance.size() != resistance.size())
      return Status::Invalid_parameters;
    callback_current.emplace_back(current.begin(), current.end());
    std::copy(ocv.begin(), ocv.end(), output_ocv.begin());
    std::copy(resistance.begin(),
              resistance.end(),
              output_resistance.begin());
    return Status::Success;
  }
};

struct AffineFailureCase
{
  std::string_view name;
  core::PackNode root;
  std::vector<double> ocv;
  std::vector<double> resistance;
  double applied_current;
  core::PackSolveMode mode;
};

#if defined(SLIDE_TEST_HAS_X86_FTZ)
class ScopedFlushToZero
{
public:
  ScopedFlushToZero() noexcept : previous_{ _mm_getcsr() }
  {
    _mm_setcsr(previous_ | _MM_FLUSH_ZERO_MASK);
  }

  ~ScopedFlushToZero() { _mm_setcsr(previous_); }

  ScopedFlushToZero(const ScopedFlushToZero &) = delete;
  ScopedFlushToZero &operator=(const ScopedFlushToZero &) = delete;

private:
  unsigned int previous_;
};
#endif

core::CompiledPackTopology compile(const core::PackNode &root)
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription({ .root = root }, topology) == Status::Success);
  return topology;
}

Status solveAffineOnce(const core::PackNode &root,
                       std::vector<double>
                         ocv,
                       std::vector<double>
                         resistance,
                       double applied_current,
                       core::PackSolveMode mode,
                       double tolerance = 1e-12,
                       int max_iterations = 2)
{
  auto topology = compile(root);
  AffineBatch batch{ .ocv = std::move(ocv),
                     .resistance = std::move(resistance) };
  const std::array views{
    core::TheveninBatchView::bind(batch, static_cast<int>(batch.ocv.size()))
  };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, views) == Status::Success);
  return solver.solve(applied_current, mode, tolerance, max_iterations);
}

[[nodiscard]] bool haveSameBits(std::span<const double> actual,
                                std::span<const double> expected) noexcept
{
  return actual.size() == expected.size()
         && std::equal(actual.begin(),
                       actual.end(),
                       expected.begin(),
                       [](double lhs, double rhs) {
                         return std::bit_cast<std::uint64_t>(lhs)
                                == std::bit_cast<std::uint64_t>(rhs);
                       });
}

template <std::size_t CurrentCount, std::size_t NodeCount>
void checkExactSolution(
  const core::PackSolution &solution,
  double terminal_voltage,
  const std::array<double, CurrentCount> &cell_current,
  const std::array<double, NodeCount> &node_voltage)
{
  CHECK(std::bit_cast<std::uint64_t>(solution.terminal_voltage)
        == std::bit_cast<std::uint64_t>(terminal_voltage));
  CHECK(haveSameBits(solution.cell_current, cell_current));
  CHECK(haveSameBits(solution.node_voltage, node_voltage));
}

[[nodiscard]] std::vector<double> independentKclResidual(
  const core::CompiledElectricalNetlist &netlist,
  const core::PackSolution &solution,
  double applied_current)
{
  std::vector<double> residual(netlist.node_count);
  for (const auto &branch : netlist.branches) {
    const auto positive = static_cast<std::size_t>(branch.node_positive);
    const auto negative = static_cast<std::size_t>(branch.node_negative);
    const double branch_current =
      branch.kind == core::ElectricalBranchKind::cell
        ? -solution.cell_current[branch.cell]
        : (solution.node_voltage[positive] - solution.node_voltage[negative])
            / branch.resistance;
    residual[positive] += branch_current;
    residual[negative] -= branch_current;
  }
  residual[netlist.terminal_positive] += applied_current;
  residual[netlist.terminal_negative] -= applied_current;
  return residual;
}

test_support::RecordedBits recordedPackSolveTrace(
  const core::CompiledPackTopology &topology,
  std::vector<double>
    ocv,
  std::vector<double>
    resistance,
  double applied_current,
  core::PackSolveMode mode,
  int maximum_iterations,
  int expected_iterations)
{
  TracingAffineBatch batch{ .ocv = std::move(ocv),
                            .resistance = std::move(resistance) };
  const auto lanes = static_cast<int>(batch.ocv.size());
  const std::array views{ core::TheveninBatchView::bind(batch, lanes) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, views) == Status::Success);
  REQUIRE(solver.solve(
            applied_current, mode, 1e-12, maximum_iterations)
          == Status::Success);
  REQUIRE(solver.diagnostics().iterations == expected_iterations);

  test_support::RecordedBits recorded;
  for (const auto &frame : batch.callback_current)
    recorded.append(frame);
  recorded.append(solver.solution().cell_current);
  recorded.append(solver.solution().node_voltage);
  const std::array scalar_frame{
    solver.solution().terminal_voltage,
    solver.diagnostics().residual_norm,
    solver.diagnostics().constraint_drift,
    solver.diagnostics().constraint_bound,
    solver.diagnostics().relaxation_gain,
  };
  recorded.append(scalar_frame);
  return recorded;
}

} // namespace

TEST_CASE("Thevenin adapters reject every malformed public shape",
          "[core][pack][thevenin][validation][coverage]")
{
  AffineBatch first{ .ocv = { 4.0, 4.1 }, .resistance = { 0.1, 0.2 } };
  AffineBatch second{ .ocv = { 3.9 }, .resistance = { 0.3 } };
  const auto two_lane = core::TheveninBatchView::bind(first, 2);
  const auto one_lane = core::TheveninBatchView::bind(second, 1);
  std::array<double, 2> current{};
  std::array<double, 2> ocv{};
  std::array<double, 2> resistance{};
  core::TheveninBatchView empty_view;
  CHECK(empty_view.linearize(current, ocv, resistance)
        == Status::Invalid_parameters);
  CHECK(two_lane.linearize(std::span<const double>{ current }.first(1),
                           std::span<double>{ ocv }.first(1),
                           std::span<double>{ resistance }.first(1))
        == Status::Invalid_parameters);

  const auto topology = compile(core::parallel(
    2, core::cell({ .archetype = "a" })));
  const std::array one_view{ two_lane };
  core::PackTheveninSystem system;
  CHECK(system.configure(topology.cells, topology.batch_archetypes, {})
        == Status::Invalid_parameters);
  CHECK(system.configure({}, topology.batch_archetypes, one_view)
        == Status::Invalid_parameters);

  auto archetypes = topology.batch_archetypes;
  archetypes[0].clear();
  CHECK(system.configure(topology.cells, archetypes, one_view)
        == Status::Invalid_parameters);

  auto cells = topology.cells;
  cells[0].location.batch = 1;
  CHECK(system.configure(cells, topology.batch_archetypes, one_view)
        == Status::Invalid_parameters);
  cells = topology.cells;
  cells[0].archetype = "wrong";
  CHECK(system.configure(cells, topology.batch_archetypes, one_view)
        == Status::Invalid_parameters);
  cells = topology.cells;
  cells[1].location.lane = 0;
  CHECK(system.configure(cells, topology.batch_archetypes, one_view)
        == Status::Invalid_parameters);
  cells = topology.cells;
  cells[1].location.lane = 2;
  CHECK(system.configure(cells, topology.batch_archetypes, one_view)
        == Status::Invalid_parameters);

  const std::array invalid_view{ empty_view };
  CHECK(system.configure(topology.cells, topology.batch_archetypes, invalid_view)
        == Status::Invalid_parameters);
  const std::array wrong_lanes{ one_lane };
  CHECK(system.configure(topology.cells, topology.batch_archetypes, wrong_lanes)
        == Status::Invalid_parameters);

  auto two_archetypes = compile(core::series(std::vector{
    core::cell({ .archetype = "a" }), core::cell({ .archetype = "b" }) }));
  two_archetypes.batch_archetypes[1] = two_archetypes.batch_archetypes[0];
  const std::array two_views{ one_lane, two_lane };
  CHECK(system.configure(two_archetypes.cells,
                         two_archetypes.batch_archetypes,
                         two_views)
        == Status::Invalid_parameters);

  REQUIRE(system.configure(topology.cells, topology.batch_archetypes, one_view)
          == Status::Success);
  core::BatchExecutor executor;
  REQUIRE(executor.configure(1, 1) == Status::Success);
  CHECK(system.linearize(std::span<const double>{ current }.first(1),
                         ocv,
                         resistance,
                         executor)
        == Status::Invalid_parameters);

  first.ocv[0] = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  CHECK(system.linearize(current, ocv, resistance, executor)
        == Status::Invalid_states);
}

TEST_CASE("pack solver rejects invalid scalar controls and incompatible modes",
          "[core][pack][solver][validation][coverage]")
{
  core::PackSolver solver;
  CHECK(solver.setRelaxationGain(0.5) == Status::Invalid_parameters);
  CHECK(solver.solve(1.0) == Status::Invalid_parameters);

  auto topology = compile(core::cell({ .archetype = "affine" }));
  AffineBatch batch{ .ocv = { 4.0 }, .resistance = { 0.1 } };
  const std::array batches{ core::TheveninBatchView::bind(batch, 1) };
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  for (const double gain : { 0.0,
                             1.1,
                             std::bit_cast<double>(
                               UINT64_C(0x7ff8000000000000)) }) {
    CAPTURE(gain);
    CHECK(solver.setRelaxationGain(gain) == Status::Invalid_parameters);
  }
  const double positive_infinity =
    std::bit_cast<double>(UINT64_C(0x7ff0000000000000));
  CHECK(solver.solve(positive_infinity)
        == Status::Invalid_parameters);
  batch.resistance[0] = 2.0;
  CHECK(solver.solve(std::numeric_limits<double>::max())
        == Status::Invalid_states);
  batch.resistance[0] = 0.1;
  CHECK(solver.solve(1.0, core::PackSolveMode::sparse_newton, 0.0)
        == Status::Invalid_parameters);
  CHECK(solver.solve(1.0, core::PackSolveMode::sparse_newton, 1e-10, 0)
        == Status::Invalid_parameters);

  auto linked = compile(core::series(
    2, core::cell({ .archetype = "affine" }), { .resistance = 0.01 }));
  AffineBatch linked_batch{ .ocv = { 4.0, 4.0 },
                            .resistance = { 0.1, 0.1 } };
  const std::array linked_batches{
    core::TheveninBatchView::bind(linked_batch, 2)
  };
  REQUIRE(solver.configure(linked, linked_batches) == Status::Success);
  CHECK(solver.solve(1.0, core::PackSolveMode::ladder)
        == Status::Invalid_parameters);

  topology.electrical.index1_candidate = false;
  REQUIRE(solver.configure(topology, batches) == Status::Success);
  CHECK(solver.solve(1.0, core::PackSolveMode::relaxation)
        == Status::Invalid_parameters);
}

TEST_CASE("all pack solve modes preserve exact high-dynamic-range affine digits",
          "[core][pack][solver][oracle][MQ.2]")
{
  constexpr double high = 0x1p53;
  constexpr double applied_current = high - 8.0;
  constexpr std::array expected_current{
    high - 2.0, -1.0, high - 3.0, -high - 2.0
  };
  constexpr std::array expected_voltage{ 2.0, 0.0 };
  constexpr std::array zero_kcl{ 0.0, 0.0 };
  const std::array modes{ core::PackSolveMode::sparse_newton,
                          core::PackSolveMode::ladder,
                          core::PackSolveMode::relaxation };

  core::CompiledPackTopology topology;
  const auto compile_status = core::compilePackDescription(
    { .root = core::parallel(
        4, core::cell({ .archetype = "affine" })) },
    topology);
  AffineBatch batch{ .ocv = { high, 1.0, high - 1.0, -high },
                     .resistance = { 1.0, 1.0, 1.0, 1.0 } };
  const std::array views{ core::TheveninBatchView::bind(batch, 4) };
  std::array<bool, modes.size()> exact_kcl{};
  for (std::size_t index = 0; index < modes.size(); ++index) {
    core::PackSolver solver;
    REQUIRE((compile_status == Status::Success
             && solver.configure(topology, views) == Status::Success));
    REQUIRE(solver.solve(applied_current, modes[index], 1e-12, 4)
            == Status::Success);
    checkExactSolution(
      solver.solution(), 2.0, expected_current, expected_voltage);
    CHECK(solver.diagnostics().iterations == 2);
    exact_kcl[index] = haveSameBits(
      independentKclResidual(
        topology.electrical, solver.solution(), applied_current),
      zero_kcl);
    if (modes[index] == core::PackSolveMode::relaxation) {
      CHECK(std::bit_cast<std::uint64_t>(
              solver.diagnostics().constraint_drift)
            == std::bit_cast<std::uint64_t>(0.0));
      CHECK(std::bit_cast<std::uint64_t>(
              solver.diagnostics().relaxation_gain)
            == std::bit_cast<std::uint64_t>(1.0));
    }
  }
  CHECK(std::all_of(exact_kcl.begin(), exact_kcl.end(), [](bool exact) {
    return exact;
  }));
}

TEST_CASE("pack resistor arm has an exact independent affine trace",
          "[core][pack][solver][resistor][oracle][MQ.2]")
{
  constexpr std::array expected_current{ 1.0 };
  constexpr std::array expected_voltage{ 1.0, 0.0, 1.5 };
  constexpr std::array zero_kcl{ 0.0, 0.0, 0.0 };
  core::CompiledPackTopology topology;
  const auto compile_status = core::compilePackDescription(
    { .root = core::parallel(
        1,
        core::cell({ .archetype = "affine" }),
        { .resistance = 0.5 }) },
    topology);
  AffineBatch batch{ .ocv = { 2.0 }, .resistance = { 0.5 } };
  const std::array views{ core::TheveninBatchView::bind(batch, 1) };
  core::PackSolver solver;
  REQUIRE((compile_status == Status::Success
           && solver.configure(topology, views) == Status::Success));

  REQUIRE(solver.solve(
            1.0, core::PackSolveMode::sparse_newton, 1e-12, 4)
          == Status::Success);
  checkExactSolution(
    solver.solution(), 1.0, expected_current, expected_voltage);
  CHECK(solver.diagnostics().iterations == 2);
  CHECK(haveSameBits(
    independentKclResidual(topology.electrical, solver.solution(), 1.0),
    zero_kcl));

  REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
  REQUIRE(solver.solve(
            1.0, core::PackSolveMode::relaxation, 1e-12, 4)
          == Status::Success);
  checkExactSolution(
    solver.solution(), 1.0, expected_current, expected_voltage);
  CHECK(solver.diagnostics().iterations == 1);
  CHECK(std::bit_cast<std::uint64_t>(
          solver.diagnostics().constraint_drift)
        == std::bit_cast<std::uint64_t>(0.0));
  CHECK(std::bit_cast<std::uint64_t>(
          solver.diagnostics().constraint_bound)
        == std::bit_cast<std::uint64_t>(1e-12));
  CHECK(std::bit_cast<std::uint64_t>(
          solver.diagnostics().relaxation_gain)
        == std::bit_cast<std::uint64_t>(1.0));
  CHECK(haveSameBits(
    independentKclResidual(topology.electrical, solver.solution(), 1.0),
    zero_kcl));
}

TEST_CASE("sparse pack solve preserves its exact damped iteration trace",
          "[core][pack][solver][damping][oracle][MQ.2]")
{
  constexpr std::array expected_current{ 4.0, -4.0 };
  constexpr std::array expected_voltage{ 4.0, 0.0 };
  constexpr std::array zero_kcl{ 0.0, 0.0 };
  core::CompiledPackTopology topology;
  const auto compile_status = core::compilePackDescription(
    { .root = core::parallel(
        2, core::cell({ .archetype = "affine" })) },
    topology);
  AffineBatch batch{ .ocv = { 8.0, 0.0 },
                     .resistance = { 1.0, 1.0 } };
  const std::array views{ core::TheveninBatchView::bind(batch, 2) };
  core::PackSolver solver;
  REQUIRE((compile_status == Status::Success
           && solver.configure(topology, views) == Status::Success));
  REQUIRE(solver.solve(
            0.0, core::PackSolveMode::sparse_newton, 1e-12, 3)
          == Status::Success);
  checkExactSolution(
    solver.solution(), 4.0, expected_current, expected_voltage);
  CHECK(solver.diagnostics().iterations == 3);
  CHECK(haveSameBits(
    independentKclResidual(topology.electrical, solver.solution(), 0.0),
    zero_kcl));
}

TEST_CASE("non-dyadic pack solve paths retain callback and publication bits",
          "[core][pack][solver][oracle][recorded][MQ.2]")
{
  const auto all_mode_topology = compile(core::parallel(
    2, core::cell({ .archetype = "affine" })));
  constexpr std::array all_mode_ocv{ 0.0, 0.0 };
  constexpr std::array all_mode_resistance{ 0.1, 0.1 };
  constexpr std::array modes{ core::PackSolveMode::sparse_newton,
                              core::PackSolveMode::ladder,
                              core::PackSolveMode::relaxation };
  std::array<test_support::RecordedBits, 4> traces;
  for (std::size_t index = 0; index < modes.size(); ++index) {
    traces[index] = recordedPackSolveTrace(
      all_mode_topology,
      { all_mode_ocv.begin(), all_mode_ocv.end() },
      { all_mode_resistance.begin(), all_mode_resistance.end() },
      0.1,
      modes[index],
      4,
      2);
  }

  const auto damped_topology = compile(core::parallel(
    2, core::cell({ .archetype = "affine" })));
  traces[3] = recordedPackSolveTrace(damped_topology,
                                     { 0.5, 0.0 },
                                     { 0.1, 0.1 },
                                     0.0,
                                     core::PackSolveMode::sparse_newton,
                                     3,
                                     3);

  constexpr std::array<std::size_t, 4> expected_values{ 13, 13, 13, 15 };
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_RELEASE) && defined(SLIDE_TEST_IPO)
  constexpr std::array expected_fnv{ UINT64_C(0xe23683484d703894),
                                     UINT64_C(0xd645bb8d493fb31c),
                                     UINT64_C(0x3168649734a681a2),
                                     UINT64_C(0xd288c3ed5d07b427) };
  constexpr std::array expected_mixed{ UINT64_C(0xa9434a138b0616b0),
                                       UINT64_C(0xd7ade2e7849420f1),
                                       UINT64_C(0x5e3b17711acbdbf4),
                                       UINT64_C(0x65c29e2fca265041) };
#elif defined(SLIDE_TEST_RELEASE)
  constexpr std::array expected_fnv{ UINT64_C(0xe23683484d703894),
                                     UINT64_C(0xd645bb8d493fb31c),
                                     UINT64_C(0x3168649734a681a2),
                                     UINT64_C(0xd288c3ed5d07b427) };
  constexpr std::array expected_mixed{ UINT64_C(0xa9434a138b0616b0),
                                       UINT64_C(0xd7ade2e7849420f1),
                                       UINT64_C(0x5e3b17711acbdbf4),
                                       UINT64_C(0x65c29e2fca265041) };
#else
  constexpr std::array expected_fnv{ UINT64_C(0xe23683484d703894),
                                     UINT64_C(0xd645bb8d493fb31c),
                                     UINT64_C(0x3168649734a681a2),
                                     UINT64_C(0xd288c3ed5d07b427) };
  constexpr std::array expected_mixed{ UINT64_C(0xa9434a138b0616b0),
                                       UINT64_C(0xd7ade2e7849420f1),
                                       UINT64_C(0x5e3b17711acbdbf4),
                                       UINT64_C(0x65c29e2fca265041) };
#endif
#endif
  for (std::size_t index = 0; index < traces.size(); ++index) {
    CAPTURE(index, traces[index].values, traces[index].fnv1a, traces[index].mixed);
    REQUIRE(traces[index].values == expected_values[index]);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
    CHECK(traces[index].fnv1a == expected_fnv[index]);
    CHECK(traces[index].mixed == expected_mixed[index]);
#endif
  }
}

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
          == Status::Invalid_states);
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

TEST_CASE("all solver modes reject a flushed or subnormal conductance",
          "[core][pack][solver][finite][fast-math][coverage]")
{
  for (const auto mode : { core::PackSolveMode::sparse_newton,
                           core::PackSolveMode::ladder,
                           core::PackSolveMode::relaxation }) {
    DYNAMIC_SECTION(static_cast<int>(mode))
    {
      CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                            { 0.0 },
                            { std::numeric_limits<double>::max() },
                            1.0,
                            mode)
            != Status::Success);
    }
  }
}

TEST_CASE("relaxation rejects a flushed reciprocal before diagonal fallback",
          "[core][pack][solver][relaxation][fast-math][mutation]")
{
#if defined(SLIDE_TEST_HAS_X86_FTZ)
  const ScopedFlushToZero flush_subnormal_results;
  CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                        { 0.0 },
                        { std::numeric_limits<double>::max() },
                        1.0,
                        core::PackSolveMode::relaxation,
                        1e-12,
                        1)
        == Status::Invalid_states);
#else
  SUCCEED("x86 flush-to-zero control is unavailable on this platform");
#endif
}

TEST_CASE("source stepping scales before multiplying an extreme finite current",
          "[core][pack][solver][source-stepping][finite][P9]")
{
  const auto applied_current = std::ldexp(1.0, 1023);
  const auto branch_current = applied_current / 2.0;
  const auto topology = compile(core::parallel(
    2, core::cell({ .archetype = "affine" })));
  AffineBatch batch{ .ocv = { 0.0, 0.0 },
                     .resistance = { 1.0, 1.0 } };
  const std::array views{ core::TheveninBatchView::bind(batch, 2) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, views) == Status::Success);

  REQUIRE(solver.solve(applied_current,
                       core::PackSolveMode::ladder,
                       applied_current / 4.0,
                       1)
          == Status::Success);
  REQUIRE(solver.solution().cell_current.size() == 2);
  CHECK(solver.solution().cell_current[0] == branch_current);
  CHECK(solver.solution().cell_current[1] == branch_current);
  CHECK(solver.solution().terminal_voltage == -branch_current);
  CHECK(solver.diagnostics().source_steps == 8);
  CHECK(solver.diagnostics().iterations == 1);
  CHECK(batch.calls == 10);
}

TEST_CASE("finite trial currents with an unrepresentable delta are rejected atomically",
          "[core][pack][solver][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();
  const auto topology = compile(core::cell({ .archetype = "affine" }));
  AffineBatch batch{ .ocv = { 0.0 }, .resistance = { 1.0 } };
  const std::array views{ core::TheveninBatchView::bind(batch, 1) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, views) == Status::Success);
  REQUIRE(solver.solve(-maximum,
                       core::PackSolveMode::ladder,
                       maximum,
                       1)
          == Status::Success);
  const auto expected = solver.solution();

  REQUIRE(solver.solve(maximum,
                       core::PackSolveMode::ladder,
                       1.0,
                       1)
          == Status::Invalid_states);
  CHECK(solver.solution().cell_current == expected.cell_current);
  CHECK(solver.solution().node_voltage == expected.node_voltage);
  CHECK(solver.solution().terminal_voltage == expected.terminal_voltage);
}

TEST_CASE("sparse residual assembly rejects each finite arithmetic overflow",
          "[core][pack][solver][mode-a][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();
  const auto smallest = std::numeric_limits<double>::denorm_min();

  SECTION("cell voltage minus OCV")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { -maximum }, .resistance = { 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.solve(0.0) == Status::Success);
    batch.ocv[0] = maximum;
    CHECK(solver.solve(0.0) == Status::Invalid_states);
  }

  SECTION("branch-current quotient")
  {
    CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                          { 1.0 },
                          { smallest },
                          0.0,
                          core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("parallel residual accumulation")
  {
    CHECK(solveAffineOnce(
            core::parallel(2, core::cell({ .archetype = "affine" })),
            { -maximum, -maximum },
            { 1.0, 1.0 },
            0.0,
            core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("terminal-current accumulation")
  {
    CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                          { -maximum },
                          { 1.0 },
                          maximum,
                          core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }
}

TEST_CASE("sparse contraction rejects an overflow between finite residuals",
          "[core][pack][solver][mode-a][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();
  const auto topology = compile(core::cell({ .archetype = "scripted" }));
  ScriptedAffineBatch batch{
    .ocv_by_call = { { 0.0 }, { -maximum } },
    .resistance_by_call = { { 1.0 } }
  };
  const std::array views{ core::TheveninBatchView::bind(batch, 1) };
  core::PackSolver solver;
  REQUIRE(solver.configure(topology, views) == Status::Success);
  CHECK(solver.solve(1e-300,
                     core::PackSolveMode::sparse_newton,
                     std::numeric_limits<double>::denorm_min(),
                     2)
        == Status::Invalid_states);
  CHECK(batch.calls == 2);
}

TEST_CASE("sparse factorization rejects non-representable and singular matrices",
          "[core][pack][solver][mode-a][finite][coverage]")
{
  const auto smallest = std::numeric_limits<double>::denorm_min();
  const auto maximum = std::numeric_limits<double>::max();

  SECTION("reciprocal conductance")
  {
    CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                          { 0.0 },
                          { smallest },
                          0.0,
                          core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("parallel coefficient accumulation")
  {
    CHECK(solveAffineOnce(
            core::parallel(2, core::cell({ .archetype = "affine" })),
            { 0.0, 0.0 },
            { 1e-308, 1e-308 },
            0.0,
            core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("rounded singular pivot")
  {
    CHECK(solveAffineOnce(
            core::series(2, core::cell({ .archetype = "affine" })),
            { 0.0, 0.0 },
            { 1.0, maximum },
            0.0,
            core::PackSolveMode::sparse_newton)
          == Status::Numerical_failure);
  }
}

TEST_CASE("sparse updates reject every non-representable finite candidate",
          "[core][pack][solver][mode-a][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();

  SECTION("node correction")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { maximum }, .resistance = { 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.solve(0.0) == Status::Success);
    batch.resistance[0] = 0.5;
    CHECK(solver.solve(-maximum) == Status::Invalid_states);
  }

  SECTION("cell-current numerator")
  {
    CHECK(solveAffineOnce(core::cell({ .archetype = "affine" }),
                          { maximum },
                          { 2.0 },
                          maximum,
                          core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("cell-current quotient")
  {
    CHECK(solveAffineOnce(
            core::parallel(2, core::cell({ .archetype = "affine" })),
            { maximum / 2.0, -maximum / 2.0 },
            { 0.5, 1.0 },
            maximum,
            core::PackSolveMode::sparse_newton)
          == Status::Invalid_states);
  }

  SECTION("current change")
  {
    const auto topology = compile(core::parallel(
      2, core::cell({ .archetype = "affine" })));
    AffineBatch batch{ .ocv = { -maximum, maximum },
                       .resistance = { 1.0, 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 2) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.solve(0.0,
                         core::PackSolveMode::sparse_newton,
                         maximum,
                         1)
            == Status::Success);
    batch.ocv = { maximum, -maximum };
    CHECK(solver.solve(0.0) == Status::Invalid_states);
  }
}

TEST_CASE("solver kernels reject an unrepresentable resistor drop before publication",
          "[core][pack][solver][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();

  SECTION("sparse candidate")
  {
    auto topology = compile(core::parallel(
      4, core::cell({ .archetype = "affine" })));
    topology.electrical = {
      .node_count = 3,
      .terminal_positive = 0,
      .terminal_negative = 1,
      .branches = {
        { .node_positive = 2, .node_negative = 1, .kind = core::ElectricalBranchKind::resistor, .resistance = 0.125 },
        { .node_positive = 0, .node_negative = 1, .kind = core::ElectricalBranchKind::cell, .cell = 0 },
        { .node_positive = 0, .node_negative = 2, .kind = core::ElectricalBranchKind::cell, .cell = 1 },
        { .node_positive = 0, .node_negative = 1, .kind = core::ElectricalBranchKind::cell, .cell = 2 },
        { .node_positive = 0, .node_negative = 2, .kind = core::ElectricalBranchKind::cell, .cell = 3 },
      },
      .nodal_sparsity = { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 2 } },
      .connected = true,
      .index1_candidate = true,
      .series_parallel_ladder = false,
    };
    AffineBatch batch{ .ocv = std::vector<double>(4),
                       .resistance = std::vector<double>(4, 0.25) };
    const std::array views{ core::TheveninBatchView::bind(batch, 4) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    for (const double scale : { 0.25, 0.5, 0.75 }) {
      batch.ocv = { -maximum * (scale / 8.0),
                    maximum * (scale / 4.0),
                    -maximum * (scale / 8.0),
                    maximum * (scale / 4.0) };
      REQUIRE(solver.solve(0.0,
                           core::PackSolveMode::sparse_newton,
                           maximum,
                           1)
              == Status::Success);
    }
    const auto expected = solver.solution();
    batch.ocv = { -maximum * (5.0 / 32.0),
                  maximum * (5.0 / 16.0),
                  -maximum * (5.0 / 32.0),
                  maximum * (5.0 / 16.0) };
    REQUIRE(solver.solve(0.0,
                         core::PackSolveMode::sparse_newton,
                         maximum,
                         1)
            == Status::Invalid_states);
    CHECK(solver.solution().cell_current == expected.cell_current);
    CHECK(solver.solution().node_voltage == expected.node_voltage);
    CHECK(solver.solution().terminal_voltage == expected.terminal_voltage);
  }

  SECTION("relaxation candidate")
  {
    auto topology = compile(core::series(
      2,
      core::cell({ .archetype = "affine" }),
      { .resistance = 4.0 }));
    const auto resistor = std::find_if(
      topology.electrical.branches.begin(),
      topology.electrical.branches.end(),
      [](const auto &branch) {
        return branch.kind == core::ElectricalBranchKind::resistor;
      });
    REQUIRE(resistor != topology.electrical.branches.end());
    std::iter_swap(topology.electrical.branches.begin(), resistor);
    AffineBatch batch{ .ocv = { -maximum, -maximum },
                       .resistance = { 1.0, 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 2) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    CHECK(solver.solve(-maximum, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }
}

TEST_CASE("ladder layers reject every non-representable affine operation",
          "[core][pack][solver][mode-b][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();
  const auto smallest = std::numeric_limits<double>::denorm_min();
  const auto single = core::cell({ .archetype = "affine" });
  const std::vector<AffineFailureCase> cases{
    { "conductance reciprocal", single, { 0.0 }, { smallest }, 0.0, core::PackSolveMode::ladder },
    { "source product", single, { maximum }, { 0.5 }, 0.0, core::PackSolveMode::ladder },
    { "conductance accumulation", core::parallel(2, single), { 0.0, 0.0 }, { 1e-308, 1e-308 }, 0.0, core::PackSolveMode::ladder },
    { "source accumulation", core::parallel(2, single), { maximum, maximum }, { 1.0, 1.0 }, 0.0, core::PackSolveMode::ladder },
    { "layer numerator", single, { maximum }, { 1.0 }, -maximum, core::PackSolveMode::ladder },
    { "cell-current numerator", single, { maximum }, { 2.0 }, maximum, core::PackSolveMode::ladder },
    { "cell-current quotient", core::parallel(2, single), { maximum / 2.0, -maximum / 2.0 }, { 0.5, 1.0 }, maximum, core::PackSolveMode::ladder },
  };
  for (const auto &test : cases) {
    DYNAMIC_SECTION(test.name)
    {
      CHECK(solveAffineOnce(test.root,
                            test.ocv,
                            test.resistance,
                            test.applied_current,
                            test.mode)
            == Status::Invalid_states);
    }
  }
}

TEST_CASE("relaxation assembly rejects every non-representable affine operation",
          "[core][pack][solver][relaxation][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();
  const auto smallest = std::numeric_limits<double>::denorm_min();
  const auto single = core::cell({ .archetype = "affine" });
  const std::vector<AffineFailureCase> cases{
    { "conductance reciprocal", single, { 0.0 }, { smallest }, 0.0, core::PackSolveMode::relaxation },
    { "source product", single, { maximum }, { 0.5 }, 0.0, core::PackSolveMode::relaxation },
    { "coefficient accumulation", core::parallel(2, single), { 0.0, 0.0 }, { 1e-308, 1e-308 }, 0.0, core::PackSolveMode::relaxation },
    { "terminal-current accumulation", single, { maximum }, { 1.0 }, -maximum, core::PackSolveMode::relaxation },
  };
  for (const auto &test : cases) {
    DYNAMIC_SECTION(test.name)
    {
      CHECK(solveAffineOnce(test.root,
                            test.ocv,
                            test.resistance,
                            test.applied_current,
                            test.mode)
            == Status::Invalid_states);
    }
  }
}

TEST_CASE("relaxation updates reject every non-representable finite candidate",
          "[core][pack][solver][relaxation][finite][coverage]")
{
  const auto maximum = std::numeric_limits<double>::max();

  SECTION("correction difference")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { 0.0 }, .resistance = { 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    REQUIRE(solver.solve(-maximum,
                         core::PackSolveMode::ladder,
                         maximum,
                         1)
            == Status::Success);
    batch.ocv[0] = maximum;
    batch.resistance[0] = 2.0;
    CHECK(solver.solve(maximum, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }

  SECTION("cell-current numerator")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { maximum }, .resistance = { 2.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    CHECK(solver.solve(maximum, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }

  SECTION("cell-current quotient")
  {
    const auto topology = compile(core::parallel(
      2, core::cell({ .archetype = "affine" })));
    AffineBatch batch{ .ocv = { maximum / 2.0, -maximum / 2.0 },
                       .resistance = { 0.5, 1.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 2) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    CHECK(solver.solve(maximum, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }

  SECTION("roundoff operation-scale sum")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { maximum }, .resistance = { 2.0 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    CHECK(solver.solve(0.0, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }

  SECTION("roundoff operation-scale quotient")
  {
    const auto topology = compile(core::cell({ .archetype = "affine" }));
    AffineBatch batch{ .ocv = { maximum / 2.0 },
                       .resistance = { 0.5 } };
    const std::array views{ core::TheveninBatchView::bind(batch, 1) };
    core::PackSolver solver;
    REQUIRE(solver.configure(topology, views) == Status::Success);
    REQUIRE(solver.setRelaxationGain(1.0) == Status::Success);
    CHECK(solver.solve(0.0, core::PackSolveMode::relaxation)
          == Status::Invalid_states);
  }
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

TEST_CASE("pack roundoff diagnostics saturate conservatively",
          "[core][pack][solver][diagnostics][coverage]")
{
  const auto epsilon = std::numeric_limits<double>::epsilon();
  const auto expected = (2.0 * epsilon) / (1.0 - 2.0 * epsilon) * 3.0
                        + 8.0 * epsilon * 5.0;
  CHECK(core::detail::conservativePackRoundoffBound(2.0 * epsilon, 3.0, 5.0)
        == expected);
  CHECK(core::detail::conservativePackRoundoffBound(1.0, 1.0, 1.0)
        == std::numeric_limits<double>::max());
  CHECK(core::detail::conservativePackRoundoffBound(
          0.5, std::numeric_limits<double>::max(), 1.0)
        == std::numeric_limits<double>::max());
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
