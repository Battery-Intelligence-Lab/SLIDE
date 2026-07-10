/**
 * @file core_Experiment_test.cpp
 * @brief Phase-5 experiment grammar, cycler, and event gates.
 */

#include "../../src/core/Experiment.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace slide;

namespace {

core::SpmBatch makeBatch(double soc = 0.55)
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(
    soc, settings::T_ENV, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
  return batch;
}

double terminalVoltage(core::SpmBatch &batch, double current)
{
  const std::array density{ current / batch.electrode_area() };
  std::array<double, 1> voltage{};
  REQUIRE(batch.terminalVoltage({ .i_app = density }, voltage) == Status::Success);
  return voltage[0];
}

} // namespace

TEST_CASE("P5-G1 documented experiment strings compile atomically",
          "[core][experiment][grammar][P5-G1]")
{
  const std::vector<std::string> source{
    "Discharge at 1 C for 1 h or until 2.7 V",
    "Charge at C/2 until 4.2 V",
    "Hold at 4.2 V until C/20",
    "Rest for 30 min",
    "Discharge at 2 W for 10 s",
    "Run US06 (A)",
    "Charge at 500 mA for 30 seconds * 2",
  };
  core::Experiment experiment;
  core::ParseDiagnostic diagnostic;
  const auto parse_status = core::Experiment::parse(source, experiment, diagnostic);
  CAPTURE(diagnostic.step, diagnostic.message);
  REQUIRE(parse_status == Status::Success);
  REQUIRE(diagnostic.message.empty());
  REQUIRE(experiment.segments.size() == 8);

  const auto &cc = experiment.segments[0];
  CHECK(cc.mode == core::ControlMode::current);
  CHECK(cc.direction == core::Direction::discharge);
  CHECK(cc.value_is_c_rate);
  CHECK(cc.value == 1.0);
  CHECK(cc.duration == 3600.0);
  CHECK(cc.voltage_limit == 2.7);
  CHECK(cc.scheduled_start == -1.0);
  CHECK(cc.sample_period == -1.0);

  const auto &charge = experiment.segments[1];
  CHECK(charge.mode == core::ControlMode::current);
  CHECK(charge.direction == core::Direction::charge);
  CHECK(charge.value == 0.5);
  CHECK(charge.voltage_limit == 4.2);

  const auto &hold = experiment.segments[2];
  CHECK(hold.mode == core::ControlMode::voltage);
  CHECK(hold.value == 4.2);
  CHECK(hold.cutoff_is_c_rate);
  CHECK(hold.current_cutoff == 0.05);
  CHECK(experiment.segments[3].duration == 1800.0);
  CHECK(experiment.segments[4].mode == core::ControlMode::power);
  CHECK(experiment.segments[5].drive_cycle == "us06");
  CHECK(experiment.segments[6].value == 0.5);
  CHECK(experiment.segments[7].source == source.back());

  core::Experiment unchanged;
  unchanged.segments.push_back({ .mode = core::ControlMode::rest, .duration = 7.0 });
  for (const std::string malformed : {
         "Fly at 1 C for 1 hour",
         "Charge at bananas until 4.2 V",
         "Rest until tomorrow",
         "Hold at 4.2 V",
         "Discharge at -1 A for 1 s",
         "Run US06 (W)",
         "Rest for 1 s * 999999999999999999999999",
       }) {
    core::ParseDiagnostic failure;
    CHECK(core::Experiment::parse(std::array{ malformed }, unchanged, failure)
          == Status::Invalid_parameters);
    CHECK_FALSE(failure.message.empty());
    REQUIRE(unchanged.segments.size() == 1);
    CHECK(unchanged.segments[0].duration == 7.0);
  }
}

TEST_CASE("P5-G1 cycler executes power rest and drive-cycle controls",
          "[core][experiment][cycler][P5-G1]")
{
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  REQUIRE(cycler.registerDriveCycle({ .name = "US06",
                                      .time = { 0.0, 0.5, 1.0 },
                                      .current = { 1.0, -1.0, 0.0 } })
          == Status::Success);
  CHECK(cycler.registerDriveCycle({ .name = "us06",
                                    .time = { 0.0, 1.0 },
                                    .current = { 0.0, 0.0 } })
        == Status::Invalid_parameters);

  const std::vector<std::string> source{
    "Discharge at 2 W for 1 s",
    "Rest for 1 s",
    "Run US06 (A)",
  };
  core::Experiment experiment;
  core::ParseDiagnostic diagnostic;
  REQUIRE(core::Experiment::parse(source, experiment, diagnostic) == Status::Success);
  core::ExperimentSolution solution;
  REQUIRE(cycler.run(experiment, 0.25, solution) == Status::Success);
  REQUIRE(solution.reason == core::TerminationReason::final_time);
  REQUIRE(solution.status == Status::Success);
  REQUIRE(solution.time.size() == 13);
  REQUIRE(solution.sample_segment.size() == solution.time.size());
  CHECK(solution.voltage[0] * solution.current[0]
        == Catch::Approx(2.0).epsilon(1e-11));
  CHECK(solution.time.back() == Catch::Approx(3.0).margin(1e-14));
  CHECK(solution.voltage[1] * solution.current[1]
        == Catch::Approx(2.0).epsilon(1e-11));
  for (std::size_t i = 5; i <= 8; ++i)
    CHECK(solution.current[i] == 0.0);
  CHECK(solution.current[9] == Catch::Approx(1.0));
  CHECK(solution.current.back() == Catch::Approx(-0.5));
}

TEST_CASE("P5-G1 voltage termination is located on the event root",
          "[core][experiment][event][P5-G1]")
{
  auto probe = makeBatch();
  constexpr double current = 8.0;
  const double before = terminalVoltage(probe, current);
  core::EulerLegacy stepper{ probe };
  const std::array density{ current / probe.electrode_area() };
  REQUIRE(stepper.step(probe, density, 0.0, 1.0) == Status::Success);
  const double after = terminalVoltage(probe, current);
  REQUIRE(after < before);
  const double limit = 0.5 * (before + after);

  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch, core::CyclerIntegrator::euler_legacy)
          == Status::Success);
  core::Experiment experiment;
  experiment.segments.push_back({ .mode = core::ControlMode::current,
                                  .direction = core::Direction::discharge,
                                  .value = current,
                                  .voltage_limit = limit });
  core::ExperimentSolution solution;
  REQUIRE(cycler.run(experiment, 1.0, solution) == Status::Success);
  REQUIRE(solution.reason == core::TerminationReason::event);
  REQUIRE(solution.time.back() > 0.0);
  REQUIRE(solution.time.back() < 1.0);
  CHECK(solution.voltage.back() == Catch::Approx(limit).margin(2e-12));

  auto immediate_batch = makeBatch();
  core::CyclerV2 immediate_cycler;
  REQUIRE(immediate_cycler.configure(
            immediate_batch, core::CyclerIntegrator::euler_legacy)
          == Status::Success);
  core::Experiment immediate;
  immediate.segments.push_back({ .mode = core::ControlMode::current,
                                 .direction = core::Direction::discharge,
                                 .value = current,
                                 .current_cutoff = current + 1.0 });
  core::ExperimentSolution immediate_solution;
  REQUIRE(immediate_cycler.run(immediate, 1.0, immediate_solution)
          == Status::Success);
  CHECK(immediate_solution.reason == core::TerminationReason::event);
  CHECK(immediate_solution.time.back() == 0.0);
  CHECK(immediate_solution.current.back() == current);
}

TEST_CASE("P5-G1 Euler CC-CV sequence remains in the legacy Cycler parity band",
          "[core][experiment][cycler][parity][P5-G1]")
{
  // Compare the current actually applied over the final interval. CyclerV2 also
  // reports the next algebraic CV current, while legacy leaves the applied value
  // in Cell_SPM::I(). Registered bands: 0.2 microvolt, 20 microampere, 0.2 uAh.
  constexpr double dt = 1.0;
  constexpr double cc_current = -8.0;
  constexpr double cc_duration = 30.0;
  constexpr double cv_duration = 30.0;

  Cell_SPM target_probe;
  REQUIRE(test_support::initialize_legacy_kokam(target_probe, 0.5, cc_current)
          == Status::Success);
  const double legacy_initial_voltage = target_probe.V();
  for (int step = 0; step < static_cast<int>(cc_duration); ++step)
    target_probe.timeStep_CC(dt);
  const double target_voltage = target_probe.V();

  Cell_SPM legacy;
  REQUIRE(test_support::initialize_legacy_kokam(legacy, 0.5, cc_current)
          == Status::Success);
  Cycler legacy_cycler{ &legacy };
  ThroughputData cc_throughput{}, cv_throughput{};
  REQUIRE(legacy_cycler.CC(cc_current, 4.2, cc_duration, dt, 0, cc_throughput)
          == Status::ReachedTimeLimit);
  REQUIRE(legacy_cycler.CV(target_voltage, 0.0, cv_duration, dt, 0, cv_throughput)
          == Status::ReachedTimeLimit);

  auto batch = makeBatch(0.5);
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch, core::CyclerIntegrator::euler_legacy)
          == Status::Success);
  core::Experiment experiment;
  experiment.segments = {
    { .mode = core::ControlMode::current,
      .direction = core::Direction::charge,
      .value = std::abs(cc_current),
      .duration = cc_duration },
    { .mode = core::ControlMode::voltage,
      .value = target_voltage,
      .duration = cv_duration },
  };
  core::ExperimentSolution solution;
  REQUIRE(cycler.run(experiment, dt, solution) == Status::Success);
  REQUIRE(solution.reason == core::TerminationReason::final_time);

  REQUIRE(solution.current.size() >= 2);
  const double applied_current = solution.current[solution.current.size() - 2];
  const double applied_voltage = terminalVoltage(batch, applied_current);
  const double voltage_error = std::abs(applied_voltage - legacy.V());
  const double current_error = std::abs(applied_current - legacy.I());
  const double charge_error = std::abs(
    batch.state().at(batch.layout().charge_throughput, 0, 0)
    - legacy.getStateObj().Ah());
  CAPTURE(target_voltage, voltage_error, current_error, legacy.V(), legacy.I(), applied_voltage, applied_current, solution.current.back(), charge_error, legacy_initial_voltage);
  CHECK(voltage_error <= 2e-7);
  CHECK(current_error <= 2e-5);
  CHECK(charge_error <= 2e-7);
  CHECK(batch.state().at(batch.layout().elapsed_time, 0, 0)
        == Catch::Approx(legacy.getStateObj().time()).margin(1e-12));
}

TEST_CASE("custom explicit implicit differential controls and terminations execute",
          "[core][experiment][custom]")
{
  auto explicit_batch = makeBatch();
  core::CyclerV2 explicit_cycler;
  REQUIRE(explicit_cycler.configure(explicit_batch) == Status::Success);
  core::Experiment custom;
  custom.segments.push_back(
    { .mode = core::ControlMode::custom_explicit,
      .direction = core::Direction::discharge,
      .duration = 20.0,
      .custom_control = [](const core::ExperimentVariables &) { return 2.5; } });
  core::ExperimentSolution custom_solution;
  REQUIRE(explicit_cycler.run(custom, 10.0, custom_solution)
          == Status::Success);

  auto standard_batch = makeBatch();
  core::CyclerV2 standard_cycler;
  REQUIRE(standard_cycler.configure(standard_batch) == Status::Success);
  core::Experiment standard;
  standard.segments.push_back({ .mode = core::ControlMode::current,
                                .direction = core::Direction::discharge,
                                .value = 2.5,
                                .duration = 20.0 });
  core::ExperimentSolution standard_solution;
  REQUIRE(standard_cycler.run(standard, 10.0, standard_solution)
          == Status::Success);
  CHECK(custom_solution.time == standard_solution.time);
  CHECK(custom_solution.current == standard_solution.current);
  CHECK(custom_solution.voltage == standard_solution.voltage);

  auto implicit_batch = makeBatch();
  core::CyclerV2 implicit_cycler;
  REQUIRE(implicit_cycler.configure(implicit_batch) == Status::Success);
  core::Experiment implicit;
  implicit.segments.push_back(
    { .mode = core::ControlMode::custom_implicit,
      .direction = core::Direction::charge,
      .duration = 1.0,
      .custom_control = [](const core::ExperimentVariables &variables) {
        return variables.voltage - 3.8;
      } });
  core::ExperimentSolution implicit_solution;
  REQUIRE(implicit_cycler.run(implicit, 1.0, implicit_solution)
          == Status::Success);
  CHECK(implicit_solution.voltage.front()
        == Catch::Approx(3.8).margin(2e-10));
  CHECK(implicit_solution.current.front() < 0.0);

  auto differential_batch = makeBatch();
  core::CyclerV2 differential_cycler;
  REQUIRE(differential_cycler.configure(differential_batch) == Status::Success);
  core::Experiment differential;
  differential.segments.push_back(
    { .mode = core::ControlMode::custom_differential,
      .direction = core::Direction::discharge,
      .duration = 2.0,
      .custom_control = [](const core::ExperimentVariables &) { return 1.0; } });
  core::ExperimentSolution differential_solution;
  REQUIRE(differential_cycler.run(differential, 1.0, differential_solution)
          == Status::Success);
  REQUIRE(differential_solution.current.size() == 3);
  CHECK(differential_solution.current.front() == 1.0);
  CHECK(differential_solution.current.back() == 2.0);

  auto event_batch = makeBatch();
  core::CyclerV2 event_cycler;
  REQUIRE(event_cycler.configure(event_batch) == Status::Success);
  core::Experiment event;
  event.segments.push_back(
    { .mode = core::ControlMode::custom_explicit,
      .direction = core::Direction::discharge,
      .duration = 20.0,
      .custom_control = [](const core::ExperimentVariables &) { return 1.0; },
      .custom_terminations = { { .name = "five-second event",
                                 .indicator = [](const core::ExperimentVariables &variables) {
                                   return 5.0 - variables.local_time;
                                 } } } });
  core::ExperimentSolution event_solution;
  REQUIRE(event_cycler.run(event, 10.0, event_solution) == Status::Success);
  CHECK(event_solution.reason == core::TerminationReason::event);
  CHECK(event_solution.termination_name == "five-second event");
  CHECK(event_solution.time.back() == Catch::Approx(5.0).margin(1e-11));
}

TEST_CASE("scheduled starts cut steps and insert exact rest gaps",
          "[core][experiment][start-time]")
{
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  core::Experiment experiment;
  experiment.segments = {
    { .mode = core::ControlMode::rest,
      .duration = 3600.0,
      .scheduled_start = 0.0 },
    { .mode = core::ControlMode::rest,
      .duration = 600.0,
      .scheduled_start = 1800.0 },
    { .mode = core::ControlMode::rest,
      .duration = 1800.0,
      .scheduled_start = 3600.0 },
    { .mode = core::ControlMode::rest, .duration = 3600.0 },
  };
  core::ExperimentSolution solution;
  REQUIRE(cycler.run(experiment, 600.0, solution) == Status::Success);
  CHECK(solution.time.back() == 9000.0);
  CHECK(std::ranges::find(solution.time, 3000.0) != solution.time.end());
  CHECK(std::ranges::find(solution.time, 3600.0) != solution.time.end());
  CHECK(std::ranges::all_of(solution.current,
                            [](double current) { return current == 0.0; }));

  experiment.segments.front().scheduled_start = -1.0;
  core::ExperimentSolution invalid;
  CHECK(cycler.run(experiment, 600.0, invalid) == Status::Invalid_parameters);
}
