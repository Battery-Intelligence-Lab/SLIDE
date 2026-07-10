/**
 * @file core_PackStepper_test.cpp
 * @brief Transactional electrical/thermal pack-step and restore gates.
 */

#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
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
  return a.iterations == b.iterations
         && a.numeric_factorizations == b.numeric_factorizations
         && a.symbolic_factorizations == b.symbolic_factorizations
         && a.jacobian_refreshes == b.jacobian_refreshes
         && a.source_steps == b.source_steps
         && a.residual_norm == b.residual_norm
         && a.constraint_drift == b.constraint_drift
         && a.constraint_bound == b.constraint_bound
         && a.relaxation_gain == b.relaxation_gain;
}

} // namespace

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
  REQUIRE(stepper.configure(topology, batches) == Status::Success);

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
    REQUIRE(stepper.configure(topology, batches) == Status::Success);

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
