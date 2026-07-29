/**
 * @file core_PackStepper_test.cpp
 * @brief Transactional electrical/thermal pack-step and restore gates.
 */

#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

using namespace slide;

namespace {

static_assert(std::is_same_v<
                decltype(std::declval<core::PackStepper &>().solver()),
                const core::PackSolver &>,
              "PackStepper must not expose mutable solver ownership");
static_assert(std::is_nothrow_move_constructible_v<core::PackStepper>);
static_assert(std::is_nothrow_move_assignable_v<core::PackStepper>);
static_assert(!std::is_copy_constructible_v<core::PackStepper>);
static_assert(!std::is_copy_assignable_v<core::PackStepper>);

core::SpmFactoryInput thermalKokam(double temperature)
{
  auto input = test_support::make_legacy_kokam_input(0.55, temperature, 298.0);
  input.design.thermal.density = 1626.0;
  input.design.thermal.heat_capacity = 750.0;
  input.design.thermal.volume = 1.0e-4;
  input.design.thermal.surface_area = 0.0;
  input.design.thermal.h_conv = 0.0;
  input.design.thermal.environment_temperature = 298.0;
  return input;
}

bool sameDiagnostics(const core::PackSolveDiagnostics &a,
                     const core::PackSolveDiagnostics &b)
{
  const auto same_scalar = [](double lhs, double rhs) {
    return std::memcmp(&lhs, &rhs, sizeof(lhs)) == 0;
  };
  return a.iterations == b.iterations
         && a.numeric_factorizations == b.numeric_factorizations
         && a.symbolic_factorizations == b.symbolic_factorizations
         && a.jacobian_refreshes == b.jacobian_refreshes
         && a.source_steps == b.source_steps
         && same_scalar(a.residual_norm, b.residual_norm)
         && same_scalar(a.constraint_drift, b.constraint_drift)
         && same_scalar(a.constraint_bound, b.constraint_bound)
         && same_scalar(a.relaxation_gain, b.relaxation_gain);
}

bool sameBits(std::span<const double> a, std::span<const double> b)
{
  return a.size() == b.size()
         && (a.empty()
             || std::memcmp(a.data(), b.data(), a.size_bytes()) == 0);
}

bool sameSolution(const core::PackSolution &a, const core::PackSolution &b)
{
  return sameBits(a.cell_current, b.cell_current)
         && sameBits(a.node_voltage, b.node_voltage)
         && std::memcmp(&a.terminal_voltage,
                        &b.terminal_voltage,
                        sizeof(a.terminal_voltage))
              == 0;
}

struct StepperSnapshot
{
  std::size_t checkpoint_size{};
  std::vector<double> checkpoint{};
  core::PackSolution solution{};
  core::PackSolveDiagnostics diagnostics{};
  std::vector<double> cell_heat{};
  std::vector<double> boundary_heat{};
  bool workspace_valid{};
  int workspace_age{};
  int numeric_factorizations{};
  int symbolic_factorizations{};
  unsigned workers{};
};

Status capture(const core::PackStepper &stepper, StepperSnapshot &snapshot)
{
  snapshot.checkpoint_size = stepper.checkpointSize();
  snapshot.checkpoint.resize(snapshot.checkpoint_size);
  const auto status = stepper.checkpoint(snapshot.checkpoint);
  if (status != Status::Success)
    return status;
  snapshot.solution = stepper.solution();
  snapshot.diagnostics = stepper.diagnostics();
  snapshot.cell_heat.assign(
    stepper.cellExternalHeat().begin(), stepper.cellExternalHeat().end());
  snapshot.boundary_heat.assign(
    stepper.boundaryHeat().begin(), stepper.boundaryHeat().end());
  const auto &workspace = stepper.solver().workspace();
  snapshot.workspace_valid = workspace.valid();
  snapshot.workspace_age = workspace.age();
  snapshot.numeric_factorizations = workspace.numericFactorizations();
  snapshot.symbolic_factorizations = workspace.symbolicFactorizations();
  snapshot.workers = stepper.batchWorkerCount();
  return Status::Success;
}

bool sameSnapshot(const StepperSnapshot &a, const StepperSnapshot &b)
{
  return a.checkpoint_size == b.checkpoint_size
         && sameBits(a.checkpoint, b.checkpoint)
         && sameSolution(a.solution, b.solution)
         && sameDiagnostics(a.diagnostics, b.diagnostics)
         && sameBits(a.cell_heat, b.cell_heat)
         && sameBits(a.boundary_heat, b.boundary_heat)
         && a.workspace_valid == b.workspace_valid
         && a.workspace_age == b.workspace_age
         && a.numeric_factorizations == b.numeric_factorizations
         && a.symbolic_factorizations == b.symbolic_factorizations
         && a.workers == b.workers;
}

bool defaultSolution(const core::PackSolution &solution)
{
  const double zero{};
  return solution.cell_current.empty() && solution.node_voltage.empty()
         && std::memcmp(
              &solution.terminal_voltage, &zero, sizeof(zero))
              == 0;
}

struct HeterogeneousThermalPack
{
  core::CompiledPackTopology topology{};
  core::SpmBatch cold{};
  core::SpmBatch hot{};
  std::array<core::SpmBatch *, 2> batches{ &cold, &hot };
  core::PackStepper stepper{};

  Status configure()
  {
    const auto root = core::parallel(std::vector{
      core::cell({ .archetype = "cold", .thermal = true }),
      core::cell({ .archetype = "hot", .thermal = true }) });
    auto status = core::compilePackDescription(
      { .root = root,
        .thermal_boundaries = { { "coolant" } },
        .thermal_links = { { "p00", "p01", 2.0 },
                           { "p01", "coolant", 0.5 } } },
      topology);
    if (status != Status::Success)
      return status;
    const core::SpmModelOptions options{ .nch = 5, .thermal = true };
    status = core::buildSpmBatch(thermalKokam(300.0), options, 1, cold);
    if (status != Status::Success)
      return status;
    status = core::buildSpmBatch(thermalKokam(310.0), options, 1, hot);
    if (status != Status::Success)
      return status;
    return stepper.configure(topology, batches, 2);
  }

  Status advance(double time)
  {
    constexpr std::array boundary{ 295.0 };
    return stepper.step(
      20.0, time, 0.1, boundary, core::PackSolveMode::ladder);
  }
};

struct OneCellPack
{
  core::CompiledPackTopology topology{};
  core::SpmBatch batch{};
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper{};

  Status configure()
  {
    auto status = core::compilePackDescription(
      { .root = core::cell({ .archetype = "old" }) }, topology);
    if (status != Status::Success)
      return status;
    status = core::buildSpmBatch(
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
      { .nch = 5 },
      1,
      batch);
    if (status != Status::Success)
      return status;
    return stepper.configure(topology, batches);
  }

  Status advance()
  {
    return stepper.step(3.0, 0.0, 0.05, {}, core::PackSolveMode::ladder);
  }
};

} // namespace

TEST_CASE("PackStepper move construction resets the source and preserves continuation",
          "[core][pack][move][ownership][MQ.2][S1.1]")
{
  HeterogeneousThermalPack source;
  HeterogeneousThermalPack control;
  REQUIRE(source.configure() == Status::Success);
  REQUIRE(control.configure() == Status::Success);
  REQUIRE(source.advance(0.0) == Status::Success);
  REQUIRE(control.advance(0.0) == Status::Success);
  StepperSnapshot before;
  StepperSnapshot control_before;
  REQUIRE(capture(source.stepper, before) == Status::Success);
  REQUIRE(capture(control.stepper, control_before) == Status::Success);
  CHECK(sameSnapshot(before, control_before));

  core::PackStepper destination{ std::move(source.stepper) };
  StepperSnapshot moved;
  REQUIRE(capture(destination, moved) == Status::Success);
  CHECK(moved.checkpoint_size == before.checkpoint_size);
  CHECK(sameBits(moved.checkpoint, before.checkpoint));
  CHECK(sameSolution(moved.solution, before.solution));
  CHECK(sameDiagnostics(moved.diagnostics, before.diagnostics));
  CHECK(sameBits(moved.cell_heat, before.cell_heat));
  CHECK(sameBits(moved.boundary_heat, before.boundary_heat));
  CHECK(moved.workspace_valid == before.workspace_valid);
  CHECK(moved.workspace_age == before.workspace_age);
  CHECK(moved.numeric_factorizations == before.numeric_factorizations);
  CHECK(moved.symbolic_factorizations == before.symbolic_factorizations);
  CHECK(moved.workers == before.workers);

  CHECK(source.stepper.checkpointSize() == 0);
  CHECK(source.stepper.batchWorkerCount() == 0);
  CHECK(defaultSolution(source.stepper.solution()));
  CHECK(sameDiagnostics(source.stepper.diagnostics(), {}));
  CHECK(source.stepper.cellExternalHeat().empty()
        && source.stepper.boundaryHeat().empty());
  CHECK_FALSE(source.stepper.solver().workspace().valid());
  CHECK(source.stepper.solver().workspace().age() == 0);
  CHECK(source.stepper.solver().workspace().numericFactorizations() == 0);
  CHECK(source.stepper.solver().workspace().symbolicFactorizations() == 0);

  std::vector<double> empty_checkpoint(source.stepper.checkpointSize());
  REQUIRE(source.stepper.checkpoint(empty_checkpoint)
          == Status::Invalid_parameters);
  CHECK(source.stepper.restore(empty_checkpoint) == Status::Invalid_parameters);
  CHECK(source.stepper.solveElectrical(0.0, core::PackSolveMode::ladder)
        == Status::Invalid_parameters);
  CHECK(source.stepper.step(0.0, 0.1, 0.1) == Status::Invalid_parameters);
  CHECK(source.stepper.stepExponential(0.0, 0.1, 0.1)
        == Status::Invalid_parameters);

  REQUIRE(destination.step(
            20.0,
            0.1,
            0.1,
            std::array{ 295.0 },
            core::PackSolveMode::ladder)
          == Status::Success);
  REQUIRE(control.advance(0.1) == Status::Success);
  StepperSnapshot continued;
  StepperSnapshot control_continued;
  REQUIRE(capture(destination, continued) == Status::Success);
  REQUIRE(capture(control.stepper, control_continued) == Status::Success);
  CHECK(sameSnapshot(continued, control_continued));
}

TEST_CASE("PackStepper move assignment replaces ownership without touching old arenas",
          "[core][pack][move][ownership][MQ.2][S1.1]")
{
  HeterogeneousThermalPack source;
  HeterogeneousThermalPack control;
  REQUIRE(source.configure() == Status::Success);
  REQUIRE(control.configure() == Status::Success);
  REQUIRE(source.advance(0.0) == Status::Success);
  REQUIRE(control.advance(0.0) == Status::Success);
  StepperSnapshot before;
  StepperSnapshot control_before;
  REQUIRE(capture(source.stepper, before) == Status::Success);
  REQUIRE(capture(control.stepper, control_before) == Status::Success);

  OneCellPack old;
  REQUIRE(old.configure() == Status::Success);
  REQUIRE(old.advance() == Status::Success);
  const std::vector<double> old_arena(
    old.batch.state().raw().begin(), old.batch.state().raw().end());

  old.stepper = std::move(source.stepper);
  StepperSnapshot moved;
  REQUIRE(capture(old.stepper, moved) == Status::Success);
  CHECK(moved.checkpoint_size == before.checkpoint_size);
  CHECK(sameBits(moved.checkpoint, before.checkpoint));
  CHECK(sameSolution(moved.solution, before.solution));
  CHECK(sameDiagnostics(moved.diagnostics, before.diagnostics));
  CHECK(sameBits(moved.cell_heat, before.cell_heat));
  CHECK(sameBits(moved.boundary_heat, before.boundary_heat));
  CHECK(moved.workspace_valid == before.workspace_valid);
  CHECK(moved.workspace_age == before.workspace_age);
  CHECK(moved.numeric_factorizations == before.numeric_factorizations);
  CHECK(moved.symbolic_factorizations == before.symbolic_factorizations);
  CHECK(moved.workers == before.workers);
  CHECK(sameBits(old.batch.state().raw(), old_arena));

  CHECK(source.stepper.checkpointSize() == 0);
  CHECK(source.stepper.batchWorkerCount() == 0);
  CHECK(defaultSolution(source.stepper.solution()));
  CHECK(sameDiagnostics(source.stepper.diagnostics(), {}));
  CHECK(source.stepper.cellExternalHeat().empty()
        && source.stepper.boundaryHeat().empty());
  CHECK_FALSE(source.stepper.solver().workspace().valid());
  CHECK(source.stepper.solver().workspace().age() == 0);
  CHECK(source.stepper.solver().workspace().numericFactorizations() == 0);
  CHECK(source.stepper.solver().workspace().symbolicFactorizations() == 0);

  std::vector<double> empty_checkpoint(source.stepper.checkpointSize());
  REQUIRE(source.stepper.checkpoint(empty_checkpoint)
          == Status::Invalid_parameters);
  CHECK(source.stepper.restore(empty_checkpoint) == Status::Invalid_parameters);
  CHECK(source.stepper.solveElectrical(0.0, core::PackSolveMode::ladder)
        == Status::Invalid_parameters);
  CHECK(source.stepper.step(0.0, 0.1, 0.1) == Status::Invalid_parameters);
  CHECK(source.stepper.stepExponential(0.0, 0.1, 0.1)
        == Status::Invalid_parameters);

  old.stepper = std::move(old.stepper);
  StepperSnapshot self_moved;
  REQUIRE(capture(old.stepper, self_moved) == Status::Success);
  CHECK(sameSnapshot(self_moved, moved));
  REQUIRE(source.stepper.configure(old.topology, old.batches)
          == Status::Success);

  REQUIRE(old.stepper.step(
            20.0,
            0.1,
            0.1,
            std::array{ 295.0 },
            core::PackSolveMode::ladder)
          == Status::Success);
  REQUIRE(control.advance(0.1) == Status::Success);
  StepperSnapshot continued;
  StepperSnapshot control_continued;
  REQUIRE(capture(old.stepper, continued) == Status::Success);
  REQUIRE(capture(control.stepper, control_continued) == Status::Success);
  CHECK(sameSnapshot(continued, control_continued));
  CHECK(sameBits(old.batch.state().raw(), old_arena));
}

TEST_CASE("compiled pack step couples thermal batches and restore invalidates the solve",
          "[core][pack][thermal][rollback][P2-G1]")
{
  const auto root = core::parallel(std::vector{
    core::cell({ .archetype = "cold", .thermal = true }),
    core::cell({ .archetype = "hot", .thermal = true }) });
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = root,
              .thermal_links = { { "p00", "p01", 2.0 } } },
            topology)
          == Status::Success);

  core::SpmBatch cold, hot;
  const core::SpmModelOptions options{ .nch = 5, .thermal = true };
  REQUIRE(core::buildSpmBatch(thermalKokam(300.0), options, 1, cold)
          == Status::Success);
  REQUIRE(core::buildSpmBatch(thermalKokam(310.0), options, 1, hot)
          == Status::Success);
  std::array<core::SpmBatch *, 2> batches{ &cold, &hot };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches, 2) == Status::Success);
  REQUIRE(stepper.batchWorkerCount() == 2);

  std::vector<double> initial(stepper.checkpointSize());
  REQUIRE(stepper.checkpoint(initial) == Status::Success);
  REQUIRE(stepper.step(20.0, 0.0, 0.1) == Status::Success);
  REQUIRE(stepper.cellExternalHeat()[0] == 20.0);
  REQUIRE(stepper.cellExternalHeat()[1] == -20.0);
  REQUIRE(stepper.cellExternalHeat()[0] + stepper.cellExternalHeat()[1] == 0.0);

  const auto first_current = stepper.solution().cell_current;
  std::vector<double> first_accepted(stepper.checkpointSize());
  REQUIRE(stepper.checkpoint(first_accepted) == Status::Success);
  REQUIRE(stepper.restore(initial) == Status::Success);
  REQUIRE_FALSE(stepper.solver().workspace().valid());
  REQUIRE(stepper.step(20.0, 0.0, 0.1) == Status::Success);
  REQUIRE(stepper.solution().cell_current == first_current);
  std::vector<double> repeated(stepper.checkpointSize());
  REQUIRE(stepper.checkpoint(repeated) == Status::Success);
  REQUIRE(std::memcmp(first_accepted.data(), repeated.data(), repeated.size() * sizeof(double))
          == 0);
}

TEST_CASE("substeps are full-dt advances under one frozen pack solve",
          "[core][pack][stepper][substeps][oracle][MQ.2][S1]")
{
  constexpr core::real_t applied_current = 8.0;
  constexpr core::real_t dt = 0.125;
  constexpr core::real_t current_tolerance = 1e-10;
  constexpr int substeps = 4;

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::cell({ .archetype = "spm" }) }, topology)
          == Status::Success);

  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  const core::SpmModelOptions options{ .nch = 5 };
  core::SpmBatch batched;
  core::SpmBatch repeated;
  REQUIRE(core::buildSpmBatch(input, options, 1, batched)
          == Status::Success);
  REQUIRE(core::buildSpmBatch(input, options, 1, repeated)
          == Status::Success);

  std::array<core::SpmBatch *, 1> batched_view{ &batched };
  core::PackStepper batched_stepper;
  core::EulerLegacy repeated_stepper;
  REQUIRE(batched_stepper.configure(topology, batched_view)
          == Status::Success);
  REQUIRE(repeated_stepper.configure(repeated) == Status::Success);

  REQUIRE(batched_stepper.step(applied_current,
                               0.0,
                               dt,
                               {},
                               core::PackSolveMode::ladder,
                               current_tolerance,
                               substeps)
          == Status::Success);
  const auto frozen_current = batched_stepper.solution().cell_current;
  REQUIRE(frozen_current.size() == 1);
  const std::array frozen_current_density{
    frozen_current[0] / repeated.electrode_area()
  };
  for (int step = 0; step < substeps; ++step)
    REQUIRE(repeated_stepper.step(repeated,
                                  frozen_current_density,
                                  static_cast<core::real_t>(step) * dt,
                                  dt)
            == Status::Success);

  REQUIRE(std::abs(frozen_current[0] - applied_current)
          <= current_tolerance);
  REQUIRE(batched.state().at(batched.layout().elapsed_time, 0, 0)
          == 0.5);
  REQUIRE(repeated.state().at(repeated.layout().elapsed_time, 0, 0)
          == 0.5);
  const auto batched_state = batched.state().raw();
  const auto repeated_state = repeated.state().raw();
  REQUIRE(batched_state.size() == repeated_state.size());
  REQUIRE(std::memcmp(batched_state.data(),
                      repeated_state.data(),
                      batched_state.size_bytes())
          == 0);
}

TEST_CASE("pack checkpoints serialize batch slices in caller order",
          "[core][pack][checkpoint][layout][oracle][MQ.2][S1.1]")
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(std::vector{
                core::cell({ .archetype = "first" }),
                core::cell({ .archetype = "second", .thermal = true }) }) },
            topology)
          == Status::Success);

  core::SpmBatch first;
  core::SpmBatch second;
  REQUIRE(core::buildSpmBatch(
            test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
            { .nch = 5 },
            1,
            first)
          == Status::Success);
  REQUIRE(core::buildSpmBatch(
            thermalKokam(303.0),
            { .nch = 12, .thermal = true },
            1,
            second)
          == Status::Success);
  std::array<core::SpmBatch *, 2> batches{ &first, &second };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Success);

  auto first_state = first.state().raw();
  auto second_state = second.state().raw();
  REQUIRE(first_state.size() != second_state.size());
  REQUIRE(stepper.checkpointSize()
          == first_state.size() + second_state.size());
  for (std::size_t i = 0; i < first_state.size(); ++i)
    first_state[i] = 1'000.0 + static_cast<double>(i);
  for (std::size_t i = 0; i < second_state.size(); ++i)
    second_state[i] = 2'000.0 + static_cast<double>(i);
  std::vector<double> expected;
  expected.reserve(first_state.size() + second_state.size());
  expected.insert(expected.end(), first_state.begin(), first_state.end());
  expected.insert(expected.end(), second_state.begin(), second_state.end());

  std::vector<double> serialized(stepper.checkpointSize());
  REQUIRE(stepper.checkpoint(serialized) == Status::Success);
  CHECK(sameBits(serialized, expected));

  std::vector<double> wire(stepper.checkpointSize());
  for (std::size_t i = 0; i < wire.size(); ++i)
    wire[i] = 10'000.0 + static_cast<double>(i);
  REQUIRE(stepper.restore(wire) == Status::Success);
  CHECK(sameBits(first.state().raw(),
                 std::span<const double>{ wire }.first(first_state.size())));
  CHECK(sameBits(
    second.state().raw(),
    std::span<const double>{ wire }.subspan(first_state.size())));
}

TEST_CASE("heterogeneous thermal substeps hold the initial assembly frozen",
          "[core][pack][thermal][substeps][oracle][MQ.2][S1.1]")
{
  constexpr double cold_initial = 300.0;
  constexpr double hot_initial = 310.0;
  constexpr double conductance = 2.0;
  constexpr double dt = 0.25;
  constexpr int substeps = 4;
  constexpr double capacity = 1626.0 * 750.0 * 1.0e-4;
  constexpr double heat = conductance * (hot_initial - cold_initial);
  constexpr double temperature_change =
    static_cast<double>(substeps) * dt * heat / capacity;

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(
                2, core::cell({ .archetype = "thermal", .thermal = true })),
              .thermal_links = { { "s00", "s01", conductance } } },
            topology)
          == Status::Success);

  auto input = thermalKokam(cold_initial);
  const core::OCVCurve zero_ocv{
    .stoichiometry = { 0.0, 1.0 },
    .value = { 0.0, 0.0 }
  };
  for (const auto domain : core::domains)
    core::domain_value(input.design.electrode, domain).active_material.ocv =
      zero_ocv;
  input.total_entropic_coefficient = {};
  input.negative_entropic_coefficient = {};
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(
            input, { .nch = 5, .thermal = true }, 2, batch)
          == Status::Success);
  batch.state().at(batch.layout().spm.temperature, 0, 1) = hot_initial;
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Success);
  REQUIRE(stepper.step(0.0,
                       0.0,
                       dt,
                       {},
                       core::PackSolveMode::ladder,
                       1e-10,
                       substeps)
          == Status::Success);

  CHECK(stepper.solution().cell_current[0] == 0.0
        && stepper.solution().cell_current[1] == 0.0);
  CHECK(stepper.cellExternalHeat()[0] == heat
        && stepper.cellExternalHeat()[1] == -heat);
  CHECK(batch.state().at(
          batch.layout().thermal.generated_heat_energy, 0, 0)
          == 0.0
        && batch.state().at(
             batch.layout().thermal.generated_heat_energy, 0, 1)
             == 0.0);
  CHECK(std::abs(batch.state().at(
                   batch.layout().spm.temperature, 0, 0)
                 - (cold_initial + temperature_change))
        <= 1e-10);
  CHECK(std::abs(batch.state().at(
                   batch.layout().spm.temperature, 0, 1)
                 - (hot_initial - temperature_change))
        <= 1e-10);
}

TEST_CASE("real multi-archetype pack steps are bit-repeatable across worker counts",
          "[core][pack][thread-pool][determinism][P9-B32]")
{
  struct Result
  {
    std::vector<double> state{};
    core::PackSolution solution{};
    core::PackSolveDiagnostics diagnostics{};
    std::vector<double> cell_heat{};
    std::vector<double> boundary_heat{};
  };

  const auto run = [&](unsigned workers) {
    const auto root = core::parallel(std::vector{
      core::cell({ .archetype = "cold", .thermal = true }),
      core::cell({ .archetype = "hot", .thermal = true }) });
    core::CompiledPackTopology topology;
    REQUIRE(core::compilePackDescription(
              { .root = root,
                .thermal_boundaries = { { "coolant" } },
                .thermal_links = { { "p00", "p01", 2.0 },
                                   { "p01", "coolant", 0.5 } } },
              topology)
            == Status::Success);
    core::SpmBatch cold, hot;
    const core::SpmModelOptions options{ .nch = 5, .thermal = true };
    REQUIRE(core::buildSpmBatch(thermalKokam(300.0), options, 1, cold)
            == Status::Success);
    REQUIRE(core::buildSpmBatch(thermalKokam(310.0), options, 1, hot)
            == Status::Success);
    std::array<core::SpmBatch *, 2> batches{ &cold, &hot };
    core::PackStepper stepper;
    REQUIRE(stepper.configure(topology, batches, workers) == Status::Success);
    REQUIRE(stepper.batchWorkerCount() == workers);
    constexpr std::array boundary{ 295.0 };
    for (int step = 0; step < 4; ++step)
      REQUIRE(stepper.step(20.0,
                           static_cast<double>(step) * 0.1,
                           0.1,
                           boundary,
                           core::PackSolveMode::ladder)
              == Status::Success);

    Result result;
    result.state.insert(
      result.state.end(), cold.state().raw().begin(), cold.state().raw().end());
    result.state.insert(
      result.state.end(), hot.state().raw().begin(), hot.state().raw().end());
    result.solution = stepper.solution();
    result.diagnostics = stepper.diagnostics();
    result.cell_heat.assign(
      stepper.cellExternalHeat().begin(), stepper.cellExternalHeat().end());
    result.boundary_heat.assign(
      stepper.boundaryHeat().begin(), stepper.boundaryHeat().end());
    return result;
  };

  const auto serial = run(1);
  const auto parallel = run(2);
  REQUIRE(parallel.state.size() == serial.state.size());
  CHECK(std::memcmp(parallel.state.data(),
                    serial.state.data(),
                    serial.state.size() * sizeof(double))
        == 0);
  REQUIRE(parallel.solution.cell_current.size()
          == serial.solution.cell_current.size());
  CHECK(std::memcmp(parallel.solution.cell_current.data(),
                    serial.solution.cell_current.data(),
                    serial.solution.cell_current.size() * sizeof(double))
        == 0);
  REQUIRE(parallel.solution.node_voltage.size()
          == serial.solution.node_voltage.size());
  CHECK(std::memcmp(parallel.solution.node_voltage.data(),
                    serial.solution.node_voltage.data(),
                    serial.solution.node_voltage.size() * sizeof(double))
        == 0);
  CHECK(std::memcmp(&parallel.solution.terminal_voltage,
                    &serial.solution.terminal_voltage,
                    sizeof(double))
        == 0);
  CHECK(sameDiagnostics(parallel.diagnostics, serial.diagnostics));
  REQUIRE(parallel.cell_heat.size() == serial.cell_heat.size());
  CHECK(std::memcmp(parallel.cell_heat.data(),
                    serial.cell_heat.data(),
                    serial.cell_heat.size() * sizeof(double))
        == 0);
  REQUIRE(parallel.boundary_heat.size() == serial.boundary_heat.size());
  CHECK(std::memcmp(parallel.boundary_heat.data(),
                    serial.boundary_heat.data(),
                    serial.boundary_heat.size() * sizeof(double))
        == 0);
}

TEST_CASE("one archetype cannot mix thermal and isothermal lanes",
          "[core][pack][compile][validation]")
{
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(std::vector{
                core::cell({ .archetype = "spm", .thermal = true }),
                core::cell({ .archetype = "spm", .thermal = false }) }) },
            topology)
          == Status::Invalid_parameters);
}

TEST_CASE("pack stepper rejects malformed public calls before touching batch state",
          "[core][pack][validation][coverage]")
{
  core::PackStepper stepper;
  core::CompiledPackTopology empty_topology;
  CHECK(stepper.configure(empty_topology, {}) == Status::Invalid_parameters);

  std::array<double, 1> state{};
  CHECK(stepper.checkpoint(state) == Status::Invalid_parameters);
  CHECK(stepper.restore(state) == Status::Invalid_parameters);
  CHECK(stepper.step(1.0, 0.0, 1.0) == Status::Invalid_parameters);

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription({ .root = core::cell() }, topology)
          == Status::Success);
  topology.cells[0].location.batch = 1;
  std::array<core::SpmBatch *, 1> batches{};
  CHECK(stepper.configure(topology, batches) == Status::Invalid_parameters);
}

TEST_CASE("parallel pack step configuration rejects aliased batch arenas",
          "[core][pack][thread-pool][alias][P9-B36]")
{
  core::SpmBatch shared, middle;
  REQUIRE(core::buildSpmBatch(
            test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
            { .nch = 5 },
            1,
            shared)
          == Status::Success);
  REQUIRE(core::buildSpmBatch(
            test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
            { .nch = 5 },
            1,
            middle)
          == Status::Success);
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(std::vector{
                core::cell({ .archetype = "a" }),
                core::cell({ .archetype = "b" }),
                core::cell({ .archetype = "c" }) }) },
            topology)
          == Status::Success);
  std::array<core::SpmBatch *, 3> batches{ &shared, &middle, &shared };
  core::PackStepper stepper;
  CHECK(stepper.configure(topology, batches, 2)
        == Status::Invalid_parameters);
}

TEST_CASE("compiled repeated ladder bricks preserve their trusted lane period",
          "[core][pack][periodic][mode-b]")
{
  constexpr int lanes = 6;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, lanes, batch) == Status::Success);
  auto collector = batch.state().row(
    batch.layout().spm.current_collector_resistance.row_begin);
  for (int lane = 0; lane < lanes; ++lane)
    collector[static_cast<std::size_t>(lane)] *= lane % 2 == 0 ? 0.9 : 1.1;

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(3,
                                   core::parallel(2,
                                                  core::cell({ .archetype = "spm" }))) },
            topology)
          == Status::Success);
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Success);
  REQUIRE(batch.trustedLanePeriod() == 2);
  REQUIRE(stepper.step(32.0, 0.0, 1.0, {}, core::PackSolveMode::ladder, 1e-10, 10)
          == Status::Success);
  for (int row = 0; row < batch.state().n_rows(); ++row) {
    const auto values = batch.state().row(row);
    REQUIRE(values[0] == values[2]);
    REQUIRE(values[0] == values[4]);
    REQUIRE(values[1] == values[3]);
    REQUIRE(values[1] == values[5]);
  }
  REQUIRE(stepper.stepExponential(32.0, 10.0, 25.0, {}, core::PackSolveMode::ladder)
          == Status::Success);
}

TEST_CASE("failed pack configuration preserves the caller lane period",
          "[core][pack][periodic][validation][P9-G4]")
{
  constexpr int lanes = 6;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(
            thermalKokam(298.0), { .nch = 5, .thermal = true }, lanes, batch)
          == Status::Success);
  auto collector = batch.state().row(
    batch.layout().spm.current_collector_resistance.row_begin);
  for (int lane = 0; lane < lanes; ++lane)
    collector[static_cast<std::size_t>(lane)] *= lane % 2 == 0 ? 0.9 : 1.1;
  REQUIRE(batch.trustedLanePeriod() == lanes);

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(
                3,
                core::parallel(2, core::cell({ .archetype = "spm" }))) },
            topology)
          == Status::Success);
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Invalid_parameters);
  REQUIRE(batch.trustedLanePeriod() == lanes);
}

TEST_CASE("failed later pack batch restores solver publication and heat diagnostics",
          "[core][pack][rollback][thermal][P9-G4]")
{
  for (const bool exponential : { false, true }) {
    CAPTURE(exponential);
    const auto root = core::series(std::vector{
      core::cell({ .archetype = "a-isothermal" }),
      core::cell({ .archetype = "b-thermal", .thermal = true }) });
    core::CompiledPackTopology topology;
    REQUIRE(core::compilePackDescription(
              { .root = root,
                .thermal_boundaries = { { "coolant" } },
                .thermal_links = { { "s01", "coolant", 1.0 } } },
              topology)
            == Status::Success);

    core::SpmBatch first, second;
    REQUIRE(core::buildSpmBatch(
              test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
              { .nch = 5 },
              1,
              first)
            == Status::Success);
    auto near_empty = thermalKokam(298.0);
    near_empty.initial_soc = 1e-4;
    REQUIRE(core::buildSpmBatch(
              near_empty, { .nch = 5, .thermal = true }, 1, second)
            == Status::Success);
    std::array<core::SpmBatch *, 2> batches{ &first, &second };
    core::PackStepper stepper;
    REQUIRE(stepper.configure(topology, batches, 2) == Status::Success);
    REQUIRE(stepper.batchWorkerCount() == 2);

    constexpr std::array baseline_boundary{ 298.0 };
    REQUIRE(stepper.step(
              0.0, 0.0, 1e-3, baseline_boundary, core::PackSolveMode::ladder)
            == Status::Success);
    const std::vector first_state(first.state().raw().begin(),
                                  first.state().raw().end());
    const std::vector second_state(second.state().raw().begin(),
                                   second.state().raw().end());
    const auto expected_solution = stepper.solution();
    const auto expected_diagnostics = stepper.diagnostics();
    const std::vector expected_cell_heat(stepper.cellExternalHeat().begin(),
                                         stepper.cellExternalHeat().end());
    const std::vector expected_boundary_heat(stepper.boundaryHeat().begin(),
                                             stepper.boundaryHeat().end());

    constexpr std::array rejected_boundary{ 310.0 };
    const auto status = exponential
                          ? stepper.stepExponential(
                              16.0,
                              1e-3,
                              1000.0,
                              rejected_boundary,
                              core::PackSolveMode::sparse_newton)
                          : stepper.step(
                              16.0,
                              1e-3,
                              1000.0,
                              rejected_boundary,
                              core::PackSolveMode::sparse_newton);
    REQUIRE(status != Status::Success);
    CHECK(std::equal(first.state().raw().begin(),
                     first.state().raw().end(),
                     first_state.begin()));
    CHECK(std::equal(second.state().raw().begin(),
                     second.state().raw().end(),
                     second_state.begin()));
    CHECK(stepper.solution().cell_current == expected_solution.cell_current);
    CHECK(stepper.solution().node_voltage == expected_solution.node_voltage);
    CHECK(stepper.solution().terminal_voltage == expected_solution.terminal_voltage);
    CHECK(sameDiagnostics(stepper.diagnostics(), expected_diagnostics));
    CHECK(std::equal(stepper.cellExternalHeat().begin(),
                     stepper.cellExternalHeat().end(),
                     expected_cell_heat.begin()));
    CHECK(std::equal(stepper.boundaryHeat().begin(),
                     stepper.boundaryHeat().end(),
                     expected_boundary_heat.begin()));
  }
}
