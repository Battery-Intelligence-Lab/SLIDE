/**
 * @file core_Experiment_test.cpp
 * @brief Phase-5 experiment grammar, cycler, and event gates.
 */

#include "../../src/core/Experiment.hpp"
#include "../support/KokamSpmFixture.hpp"
#include "../support/RecordedBits.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <iterator>
#include <numeric>
#include <stdexcept>
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

void recordParsedSegment(test_support::RecordedBits &recorded,
                         const core::ExperimentSegment &segment)
{
  // Keep enums, flags, and vector cardinalities in the numerical trace.  Each
  // segment is a separate frame, so moving a field or a repeated segment is
  // distinguishable even when the flattened values happen to be identical.
  const std::array fields{
    static_cast<double>(static_cast<int>(segment.mode)),
    static_cast<double>(static_cast<int>(segment.direction)),
    segment.value,
    segment.value_is_c_rate ? 1.0 : 0.0,
    segment.duration,
    segment.voltage_limit,
    segment.current_cutoff,
    segment.cutoff_is_c_rate ? 1.0 : 0.0,
    segment.custom_control ? 1.0 : 0.0,
    static_cast<double>(segment.custom_terminations.size()),
    segment.scheduled_start,
    segment.sample_period,
  };
  recorded.append(fields);
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

  // 9C-3 pre-refactor parser fixture.  Numeric and integer-valued metadata are
  // fingerprinted independently of text so a parser change cannot hide a
  // source-order or normalisation regression behind the same floating trace.
  test_support::RecordedBits parsed_bits;
  std::vector<std::string> parsed_sources;
  std::vector<std::string> parsed_drive_cycles;
  for (const auto &segment : experiment.segments) {
    recordParsedSegment(parsed_bits, segment);
    parsed_sources.push_back(segment.source);
    parsed_drive_cycles.push_back(segment.drive_cycle);
  }
  const std::vector<std::string> expected_sources{
    source[0],
    source[1],
    source[2],
    source[3],
    source[4],
    source[5],
    source[6],
    source[6],
  };
  const std::vector<std::string> expected_drive_cycles{
    "",
    "",
    "",
    "",
    "",
    "us06",
    "",
    "",
  };
  CHECK(parsed_sources == expected_sources);
  CHECK(parsed_drive_cycles == expected_drive_cycles);
  CAPTURE(parsed_bits.values, parsed_bits.fnv1a, parsed_bits.mixed);
  REQUIRE(parsed_bits.values == 96);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
  // Filled from the pre-split production TU in the three supported capture
  // configurations before any Experiment.cpp ownership is moved.
#if defined(SLIDE_TEST_RELEASE) && defined(SLIDE_TEST_IPO)
  constexpr auto expected_parser_fnv = UINT64_C(0x1e404a3be7365df5);
  constexpr auto expected_parser_mixed = UINT64_C(0x94a28984eed0159b);
#elif defined(SLIDE_TEST_RELEASE)
  constexpr auto expected_parser_fnv = UINT64_C(0x1e404a3be7365df5);
  constexpr auto expected_parser_mixed = UINT64_C(0x94a28984eed0159b);
#else
  constexpr auto expected_parser_fnv = UINT64_C(0x1e404a3be7365df5);
  constexpr auto expected_parser_mixed = UINT64_C(0x94a28984eed0159b);
#endif
  CHECK(parsed_bits.fnv1a == expected_parser_fnv);
  CHECK(parsed_bits.mixed == expected_parser_mixed);
#endif

  core::Experiment unchanged;
  unchanged.segments.push_back({ .mode = core::ControlMode::rest, .duration = 7.0 });
  for (const std::string malformed : {
         "Fly at 1 C for 1 hour",
         "Charge at bananas until 4.2 V",
         "Rest until tomorrow",
         "Hold at 4.2 V",
         "Hold at 4.2 V until 3.8 V",
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

  // The numeric token itself is finite, but conversion from hours to seconds
  // overflows. Parser success here would feed an infinite execution horizon.
  const std::string overflowing_duration =
    "Rest for " + std::string(306, '9') + " hours";
  core::ParseDiagnostic overflow;
  CHECK(core::Experiment::parse(
          std::array{ overflowing_duration }, unchanged, overflow)
        == Status::Invalid_parameters);
  CHECK_FALSE(overflow.message.empty());
  REQUIRE(unchanged.segments.size() == 1);
  CHECK(unchanged.segments[0].duration == 7.0);

  // Each token is below the legacy one-million repetition ceiling, but the
  // aggregate expansion is an untrusted allocation request. The parser's
  // documented short-experiment surface is capped at 10,000 segments.
  core::Experiment bounded;
  bounded.segments.push_back({ .mode = core::ControlMode::rest, .duration = 9.0 });
  core::ParseDiagnostic expansion;
  const std::vector<std::string> expansion_bomb{
    "Rest for 1 s * 6000",
    "Rest for 1 s * 6000",
  };
  CHECK(core::Experiment::parse(expansion_bomb, bounded, expansion)
        == Status::Invalid_parameters);
  CHECK(expansion.message.find("expanded segment limit") != std::string::npos);
  REQUIRE(bounded.segments.size() == 1);
  CHECK(bounded.segments[0].duration == 9.0);

  core::Experiment retained;
  retained.segments.push_back({ .mode = core::ControlMode::rest, .duration = 11.0 });
  core::ParseDiagnostic retained_diagnostic;
  const std::string retained_bomb = "Rest for 1 s"
                                    + std::string(500, ' ')
                                    + "* 10000";
  CHECK(core::Experiment::parse(
          std::array{ retained_bomb }, retained, retained_diagnostic)
        == Status::Invalid_parameters);
  CHECK(retained_diagnostic.message.find("expanded text limit")
        != std::string::npos);
  CHECK(retained.segments.size() == 1);
  if (retained.segments.size() == 1)
    CHECK(retained.segments[0].duration == 11.0);

  core::Experiment named;
  named.segments.push_back({ .mode = core::ControlMode::rest, .duration = 13.0 });
  core::ParseDiagnostic name_diagnostic;
  const std::string long_drive = "Run " + std::string(1025, 'x') + " (A)";
  CHECK(core::Experiment::parse(
          std::array{ long_drive }, named, name_diagnostic)
        == Status::Invalid_parameters);
  CHECK(name_diagnostic.message.find("drive-cycle name") != std::string::npos);
  CHECK(named.segments.size() == 1);
  if (named.segments.size() == 1)
    CHECK(named.segments[0].duration == 13.0);

  core::Experiment missing_drive_name = named;
  core::ParseDiagnostic missing_name_diagnostic;
  CHECK(core::Experiment::parse(
          std::array{ std::string{ "Run (A)" } },
          missing_drive_name,
          missing_name_diagnostic)
        == Status::Invalid_parameters);
  CHECK_FALSE(missing_name_diagnostic.message.empty());
  CHECK(missing_drive_name.segments.size() == 1);
  if (missing_drive_name.segments.size() == 1)
    CHECK(missing_drive_name.segments[0].duration == 13.0);
}

TEST_CASE("9C-3 manually constructed cycler trace retains its pre-split bits",
          "[core][experiment][cycler][9C-3][recorded]")
{
  // This fixture deliberately bypasses Experiment::parse: correlated parser
  // and runner changes must not be able to bless one another.  The sign-changing
  // nonuniform drive samples distinguish interpolation and exact breakpoint
  // ownership; the final current segment crosses a dyadic quarter-second root
  // inside a half-second step, forcing event rollback and bisection.
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  REQUIRE(cycler.registerDriveCycle(
            { .name = "9c3-breakpoints",
              .time = { 0.0, 0.25, 0.5 },
              .current = { 0.25, -0.75, 0.5 } })
          == Status::Success);

  core::Experiment experiment;
  experiment.segments = {
    { .mode = core::ControlMode::power,
      .direction = core::Direction::discharge,
      .value = 1.25,
      .duration = 0.5,
      .sample_period = 0.25 },
    { .mode = core::ControlMode::rest,
      .duration = 0.25,
      .sample_period = 0.25 },
    { .mode = core::ControlMode::drive_cycle,
      .drive_cycle = "9c3-breakpoints",
      .sample_period = 0.125 },
    { .mode = core::ControlMode::current,
      .direction = core::Direction::discharge,
      .value = 0.4,
      .duration = 1.0,
      .custom_terminations = {
        { .name = "quarter-step event",
          .indicator = [](const core::ExperimentVariables &variables) {
            return 0.25 - variables.local_time;
          } },
      },
      .sample_period = 0.5 },
  };

  test_support::RecordedBits run_bits;
  run_bits.append(batch.state().raw());
  core::ExperimentSolution solution;
  const auto run_status = cycler.run(experiment, 0.5, solution);
  REQUIRE(run_status == Status::Success);
  REQUIRE(solution.status == Status::Success);
  REQUIRE(solution.reason == core::TerminationReason::event);
  REQUIRE(solution.segment == 3);
  REQUIRE(solution.termination_name == "quarter-step event");
  REQUIRE(solution.time.size() == 9);
  REQUIRE(solution.time.back() == 1.5);

  run_bits.append(solution.time);
  run_bits.append(solution.voltage);
  run_bits.append(solution.current);
  std::vector<double> recorded_segments;
  recorded_segments.reserve(solution.sample_segment.size());
  std::ranges::transform(solution.sample_segment,
                         std::back_inserter(recorded_segments),
                         [](std::size_t segment) {
                           return static_cast<double>(segment);
                         });
  run_bits.append(recorded_segments);
  const std::array metadata{
    static_cast<double>(static_cast<int>(solution.reason)),
    static_cast<double>(static_cast<int>(solution.status)),
    static_cast<double>(solution.segment),
  };
  run_bits.append(metadata);
  run_bits.append(batch.state().raw());

  CAPTURE(run_bits.values, run_bits.fnv1a, run_bits.mixed);
  REQUIRE(run_bits.values == 503);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_RELEASE) && defined(SLIDE_TEST_IPO)
  constexpr auto expected_run_fnv = UINT64_C(0x9d97787b3976978b);
  constexpr auto expected_run_mixed = UINT64_C(0x4c7bcd0565164aa2);
#elif defined(SLIDE_TEST_RELEASE)
  constexpr auto expected_run_fnv = UINT64_C(0xe9616b6ce3131bb6);
  constexpr auto expected_run_mixed = UINT64_C(0xb72cecb529e7c3c0);
#else
  constexpr auto expected_run_fnv = UINT64_C(0xc779e41caac8339c);
  constexpr auto expected_run_mixed = UINT64_C(0x32dec619a41e0324);
#endif
  CHECK(run_bits.fnv1a == expected_run_fnv);
  CHECK(run_bits.mixed == expected_run_mixed);
#endif
}

TEST_CASE("Cycler validates direct segment metadata before stepping",
          "[core][experiment][validation][P9]")
{
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  auto invalid_integrator_batch = makeBatch();
  core::CyclerV2 invalid_integrator_cycler;
  CHECK(invalid_integrator_cycler.configure(
          invalid_integrator_batch,
          static_cast<core::CyclerIntegrator>(255))
        == Status::Invalid_parameters);
  const std::vector<double> initial(batch.state().raw().begin(),
                                    batch.state().raw().end());
  core::ExperimentSolution sentinel;
  sentinel.time = { 123.0 };

  const double infinity = std::bit_cast<double>(UINT64_C(0x7ff0000000000000));
  const double nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  for (const core::ExperimentSegment invalid : {
         core::ExperimentSegment{ .mode = core::ControlMode::current,
                                  .direction = core::Direction::discharge,
                                  .value = 1.0,
                                  .duration = infinity,
                                  .current_cutoff = 2.0 },
         core::ExperimentSegment{ .mode = core::ControlMode::current,
                                  .direction = core::Direction::discharge,
                                  .value = nan,
                                  .duration = 1.0 },
         core::ExperimentSegment{ .mode = core::ControlMode::rest,
                                  .duration = 1.0,
                                  .sample_period = infinity },
         core::ExperimentSegment{ .mode = core::ControlMode::rest,
                                  .duration = 1.0,
                                  .scheduled_start = infinity },
         core::ExperimentSegment{ .mode = static_cast<core::ControlMode>(255),
                                  .duration = 1.0 },
       }) {
    core::Experiment experiment;
    experiment.segments.push_back(invalid);
    core::ExperimentSolution output = sentinel;
    CHECK(cycler.run(experiment, 1.0, output) == Status::Invalid_parameters);
    CHECK(output.time == sentinel.time);
    REQUIRE(batch.state().raw().size() == initial.size());
    CHECK(std::memcmp(batch.state().raw().data(), initial.data(), initial.size() * sizeof(double))
          == 0);
  }
}

TEST_CASE("Cycler rolls back a step when post-advance event evaluation fails",
          "[core][experiment][rollback][P9]")
{
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  const std::vector<double> initial(batch.state().raw().begin(),
                                    batch.state().raw().end());

  core::Experiment experiment;
  experiment.segments.push_back(
    { .mode = core::ControlMode::custom_explicit,
      .direction = core::Direction::discharge,
      .duration = 2.0,
      .custom_control = [](const core::ExperimentVariables &) { return 1.0; },
      .custom_terminations = {
        { .name = "injected callback failure",
          .indicator = [](const core::ExperimentVariables &variables) {
            if (variables.local_time > 0.0)
              throw std::runtime_error("injected post-advance failure");
            return 1.0;
          } } } });
  core::ExperimentSolution output;
  CHECK(cycler.run(experiment, 1.0, output) == Status::Invalid_parameters);
  CHECK(output.reason == core::TerminationReason::error);
  REQUIRE(batch.state().raw().size() == initial.size());
  CHECK(std::memcmp(batch.state().raw().data(), initial.data(), initial.size() * sizeof(double))
        == 0);
}

TEST_CASE("event bisection selects the first of two roots and exact breakpoints",
          "[core][experiment][event][P9]")
{
  auto first_batch = makeBatch();
  core::CyclerV2 first_cycler;
  REQUIRE(first_cycler.configure(first_batch) == Status::Success);
  core::Experiment first;
  first.segments.push_back(
    { .mode = core::ControlMode::rest,
      .duration = 1.0,
      .custom_terminations = {
        { .name = "quarter",
          .indicator = [](const core::ExperimentVariables &variables) {
            return 0.25 - variables.local_time;
          } },
        { .name = "three quarters", .indicator = [](const core::ExperimentVariables &variables) {
           return 0.75 - variables.local_time;
         } },
      } });
  core::ExperimentSolution first_solution;
  REQUIRE(first_cycler.run(first, 1.0, first_solution) == Status::Success);
  CHECK(first_solution.reason == core::TerminationReason::event);
  CHECK(first_solution.termination_name == "quarter");
  CHECK(first_solution.time.back() == Catch::Approx(0.25).margin(1e-12));

  auto breakpoint_batch = makeBatch();
  core::CyclerV2 breakpoint_cycler;
  REQUIRE(breakpoint_cycler.configure(breakpoint_batch) == Status::Success);
  core::Experiment breakpoint;
  breakpoint.segments.push_back(
    { .mode = core::ControlMode::rest,
      .duration = 2.0,
      .custom_terminations = {
        { .name = "breakpoint",
          .indicator = [](const core::ExperimentVariables &variables) {
            return 1.0 - variables.local_time;
          } },
      } });
  core::ExperimentSolution breakpoint_solution;
  REQUIRE(breakpoint_cycler.run(breakpoint, 1.0, breakpoint_solution)
          == Status::Success);
  CHECK(breakpoint_solution.reason == core::TerminationReason::event);
  CHECK(breakpoint_solution.termination_name == "breakpoint");
  CHECK(breakpoint_solution.time.back() == 1.0);
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
  // in Cell_SPM::I(). AUD-2 found that commit 09e8ec5 introduced the following
  // empirical regression envelopes with the implementation; they were not
  // registered before the decisive run. The 0.2-uAh scale is consistent with
  // 20 uA sustained for 30 s (= 0.1667 uAh), but this test compares only the
  // final applied current, not a trajectory-wide current bound.
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

TEST_CASE("Experiment public limits and drive-cycle tables reject exact boundaries",
          "[core][experiment][validation][coverage]")
{
  core::Experiment sentinel;
  sentinel.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 7.0 });
  core::ParseDiagnostic diagnostic;

  const std::vector<std::string> empty;
  CHECK(core::Experiment::parse(empty, sentinel, diagnostic)
        == Status::Invalid_parameters);
  CHECK(sentinel.segments.size() == 1);

  const std::vector<std::string> too_many(10'001, "Rest for 1 s");
  CHECK(core::Experiment::parse(too_many, sentinel, diagnostic)
        == Status::Invalid_parameters);
  CHECK(sentinel.segments.size() == 1);

  const std::array oversized_step{ std::string(65'537, 'x') };
  CHECK(core::Experiment::parse(oversized_step, sentinel, diagnostic)
        == Status::Invalid_parameters);
  CHECK(sentinel.segments.size() == 1);

  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  CHECK(cycler.registerDriveCycle({}) == Status::Invalid_parameters);
  CHECK(cycler.registerDriveCycle(
          { .name = "nonmonotone",
            .time = { 0.0, 1.0, 1.0 },
            .current = { 0.0, 1.0, 2.0 } })
        == Status::Invalid_parameters);
}

TEST_CASE("Cycler failure solvers are reached through validated public segments",
          "[core][experiment][solver][coverage]")
{
  core::CyclerV2 unconfigured;
  core::Experiment rest;
  rest.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 1.0 });
  core::ExperimentSolution output;
  CHECK(unconfigured.run(rest, 1.0, output) == Status::Invalid_parameters);

  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);

  core::Experiment decreasing_schedule;
  decreasing_schedule.segments = {
    { .mode = core::ControlMode::rest,
      .duration = 1.0,
      .scheduled_start = 0.0 },
    { .mode = core::ControlMode::rest,
      .duration = 1.0,
      .scheduled_start = 2.0 },
    { .mode = core::ControlMode::rest,
      .duration = 1.0,
      .scheduled_start = 1.0 },
  };
  CHECK(cycler.run(decreasing_schedule, 1.0, output)
        == Status::Invalid_parameters);

  core::Experiment missing_cycle;
  missing_cycle.segments.push_back(
    { .mode = core::ControlMode::drive_cycle,
      .drive_cycle = "missing" });
  CHECK(cycler.run(missing_cycle, 1.0, output)
        == Status::Invalid_parameters);

  const double maximum = std::numeric_limits<double>::max();
  core::Experiment invalid_voltage_iteration;
  invalid_voltage_iteration.segments.push_back(
    { .mode = core::ControlMode::voltage,
      .value = maximum,
      .duration = 1.0 });
  CHECK(cycler.run(invalid_voltage_iteration, 1.0, output)
        == Status::Invalid_states);

  auto cycling_power_batch = makeBatch(0.001);
  core::CyclerV2 cycling_power_cycler;
  REQUIRE(cycling_power_cycler.configure(cycling_power_batch)
          == Status::Success);
  core::Experiment cycling_power;
  cycling_power.segments.push_back(
    { .mode = core::ControlMode::power,
      .direction = core::Direction::discharge,
      .value = std::exp2(17.0 / 4.0),
      .duration = 1.0 });
  CHECK(cycling_power_cycler.run(cycling_power, 1.0, output)
        == Status::Numerical_failure);

  auto cycling_voltage_batch = makeBatch(0.1);
  core::CyclerV2 cycling_voltage_cycler;
  REQUIRE(cycling_voltage_cycler.configure(cycling_voltage_batch)
          == Status::Success);
  core::Experiment cycling_voltage;
  cycling_voltage.segments = {
    { .mode = core::ControlMode::current,
      .direction = core::Direction::discharge,
      .value = std::exp2(23.0 / 4.0),
      .duration = 1e-12 },
    { .mode = core::ControlMode::voltage,
      .value = 3.5,
      .duration = 1.0 },
  };
  CHECK(cycling_voltage_cycler.run(cycling_voltage, 1e-12, output)
        == Status::Numerical_failure);

  auto singular_input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  // Make the passive ohmic term dominant so the physical maximum-power
  // tangent lies comfortably inside the model's valid concentration range.
  singular_input.initial_current_collector_resistance =
    singular_input.design.electrode_area;
  core::SpmBatch singular_power_batch;
  REQUIRE(core::buildSpmBatch(singular_input, {}, 1, singular_power_batch)
          == Status::Success);
  const auto power_derivative = [&](double target_power) {
    const double current = target_power / 3.7;
    std::array<double, 1> input{ current }, intercept{}, resistance{};
    REQUIRE(singular_power_batch.linearizeThevenin(
              input, intercept, resistance)
            == Status::Success);
    return intercept[0] - 2.0 * resistance[0] * current;
  };
  double lower_power = 0.0;
  double upper_power = 8.0;
  REQUIRE(power_derivative(lower_power) > 0.0);
  REQUIRE(power_derivative(upper_power) < 0.0);
  for (int iteration = 0; iteration < 80; ++iteration) {
    const double middle = std::midpoint(lower_power, upper_power);
    if (power_derivative(middle) > 0.0)
      lower_power = middle;
    else
      upper_power = middle;
  }
  const double singular_power =
    std::abs(power_derivative(lower_power))
        < std::abs(power_derivative(upper_power))
      ? lower_power
      : upper_power;
  REQUIRE(std::abs(power_derivative(singular_power)) <= 1e-12);
  core::CyclerV2 singular_power_cycler;
  REQUIRE(singular_power_cycler.configure(singular_power_batch)
          == Status::Success);
  core::Experiment singular_power_experiment;
  singular_power_experiment.segments.push_back(
    { .mode = core::ControlMode::power,
      .direction = core::Direction::discharge,
      .value = singular_power,
      .duration = 1.0 });
  CHECK(singular_power_cycler.run(singular_power_experiment, 1.0, output)
        == Status::Numerical_failure);

  const double nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  core::Experiment nonfinite_custom;
  nonfinite_custom.segments.push_back(
    { .mode = core::ControlMode::custom_explicit,
      .duration = 1.0,
      .custom_control = [nan](const core::ExperimentVariables &) {
        return nan;
      } });
  CHECK(cycler.run(nonfinite_custom, 1.0, output)
        == Status::Invalid_states);

  core::Experiment singular_custom;
  singular_custom.segments.push_back(
    { .mode = core::ControlMode::custom_implicit,
      .duration = 1.0,
      .custom_control = [](const core::ExperimentVariables &) {
        return 1.0;
      } });
  CHECK(cycler.run(singular_custom, 1.0, output)
        == Status::Numerical_failure);

  core::Experiment cycling_newton;
  cycling_newton.segments.push_back(
    { .mode = core::ControlMode::custom_implicit,
      .duration = 1.0,
      .custom_control = [](const core::ExperimentVariables &variables) {
        const double current = variables.current;
        // The finite-difference Newton map is exactly 0 -> 1 -> 0:
        // f(x)=2-2x around zero and f(x)=x around one.
        return current < 0.5 ? 2.0 - 2.0 * current : current;
      } });
  CHECK(cycler.run(cycling_newton, 1.0, output)
        == Status::Numerical_failure);

  core::Experiment unbounded_event;
  unbounded_event.segments.push_back(
    { .mode = core::ControlMode::rest,
      .custom_terminations = {
        { .name = "never",
          .indicator = [](const core::ExperimentVariables &) { return 1.0; } } } });
  CHECK(cycler.run(unbounded_event, 7.0 * 24.0 * 3600.0, output)
        == Status::Numerical_failure);
  CHECK(output.termination_name == "maximum step duration");
}

TEST_CASE("Cycler revalidates callback-mutable segment controls",
          "[core][experiment][callback][validation][coverage]")
{
  auto batch = makeBatch();
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  core::ExperimentSolution output;

  SECTION("a later power direction cannot be invalidated")
  {
    core::Experiment experiment;
    experiment.segments.push_back(
      { .mode = core::ControlMode::custom_explicit,
        .duration = 1.0,
        .custom_control = [&experiment](const core::ExperimentVariables &) {
          experiment.segments[1].direction = core::Direction::none;
          return 0.0;
        } });
    experiment.segments.push_back(
      { .mode = core::ControlMode::power,
        .direction = core::Direction::discharge,
        .value = 1.0,
        .duration = 1.0 });
    CHECK(cycler.run(experiment, 1.0, output)
          == Status::Invalid_parameters);
  }

  SECTION("a later custom controller cannot be removed")
  {
    core::Experiment experiment;
    experiment.segments.push_back(
      { .mode = core::ControlMode::custom_explicit,
        .duration = 1.0,
        .custom_control = [&experiment](const core::ExperimentVariables &) {
          experiment.segments[1].custom_control = {};
          return 0.0;
        } });
    experiment.segments.push_back(
      { .mode = core::ControlMode::custom_implicit,
        .duration = 1.0,
        .custom_control = [](const core::ExperimentVariables &variables) {
          return variables.current;
        } });
    CHECK(cycler.run(experiment, 1.0, output)
          == Status::Invalid_parameters);
  }

  SECTION("a later voltage event direction cannot be invalidated")
  {
    core::Experiment experiment;
    experiment.segments.push_back(
      { .mode = core::ControlMode::custom_explicit,
        .duration = 1.0,
        .custom_control = [&experiment](const core::ExperimentVariables &) {
          experiment.segments[1].direction = core::Direction::none;
          return 0.0;
        } });
    experiment.segments.push_back(
      { .mode = core::ControlMode::current,
        .direction = core::Direction::discharge,
        .value = 1.0,
        .duration = 1.0,
        .voltage_limit = 3.0 });
    CHECK(cycler.run(experiment, 1.0, output)
          == Status::Invalid_parameters);
  }

  SECTION("an active controller removed between samples is rejected")
  {
    core::Experiment experiment;
    experiment.segments.push_back(
      { .mode = core::ControlMode::custom_explicit,
        .duration = 2.0,
        .custom_control = [](const core::ExperimentVariables &) {
          return 0.0;
        },
        .custom_terminations = {
          { .name = "mutate after control evaluation",
            .indicator = [&experiment](const core::ExperimentVariables &) {
              experiment.segments.front().custom_control = {};
              return 1.0;
            } } } });
    CHECK(cycler.run(experiment, 1.0, output)
          == Status::Invalid_parameters);
  }
}
