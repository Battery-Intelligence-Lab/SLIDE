/**
 * @file core_PackStepper_test.cpp
 * @brief Transactional electrical/thermal pack-step and restore gates.
 */

#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstring>
#include <vector>

using namespace slide;

namespace {

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
