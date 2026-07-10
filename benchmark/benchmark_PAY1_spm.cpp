/**
 * @file benchmark_PAY1_spm.cpp
 * @brief PAY-1: 10,000 Kokam SPM cells, 1C discharge, one simulated hour.
 *
 * Wall-clock output is admissible only on a quiet machine. The executable always validates
 * equal work and final-state parity so an invalid or shortened run cannot claim a speedup.
 */

#include "../src/core/EulerLegacy.hpp"
#include "../tests/support/KokamSpmFixture.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <span>
#include <string_view>
#include <vector>

using namespace slide;

namespace {

struct Options
{
  int lanes{ 10'000 };
  int steps{ 3'600 };
  int repetitions{ 3 };
};

bool parse_positive(std::string_view text, int &value)
{
  int parsed{};
  const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size() || parsed <= 0)
    return false;
  value = parsed;
  return true;
}

bool parse_options(int argc, char **argv, Options &options)
{
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 >= argc)
      return false;
    const std::string_view name{ argv[i] };
    if (name == "--lanes") {
      if (!parse_positive(argv[i + 1], options.lanes))
        return false;
    } else if (name == "--steps") {
      if (!parse_positive(argv[i + 1], options.steps))
        return false;
    } else if (name == "--repetitions") {
      if (!parse_positive(argv[i + 1], options.repetitions))
        return false;
    } else {
      return false;
    }
  }
  return true;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  const auto middle = values.size() / 2;
  return values.size() % 2 == 0 ? 0.5 * (values[middle - 1] + values[middle])
                                : values[middle];
}

template <class Function>
double seconds(Function &&function)
{
  const auto start = std::chrono::steady_clock::now();
  function();
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

int main(int argc, char **argv)
{
  Options options;
  if (!parse_options(argc, argv, options)) {
    std::fprintf(stderr,
                 "usage: benchmark_PAY1_spm [--lanes N] [--steps N] [--repetitions N]\n");
    return 2;
  }

  constexpr double initial_soc = 1.0;
  constexpr double current_A = 16.0;
  constexpr double dt = 1.0;

  Cell_SPM prototype;
  if (test_support::initialize_legacy_kokam(prototype, initial_soc, current_A)
      != Status::Success)
    return 3;
  double environment_temperature{}, reference_temperature{};
  prototype.getTemperatures(&environment_temperature, &reference_temperature);
  const auto input = test_support::make_legacy_kokam_input(
    initial_soc, prototype.T(), reference_temperature);

  core::SpmBatch batch;
  if (core::buildSpmBatch(input, {}, options.lanes, batch) != Status::Success)
    return 4;
  core::EulerLegacy stepper;
  if (stepper.configure(batch) != Status::Success)
    return 5;
  const std::vector<double> initial_arena(batch.state().raw().begin(),
                                          batch.state().raw().end());
  std::vector<Cell_SPM> legacy(static_cast<std::size_t>(options.lanes), prototype);
  const std::vector<double> current_density(
    static_cast<std::size_t>(options.lanes), current_A / batch.electrode_area());

  std::vector<double> core_seconds;
  std::vector<double> legacy_seconds;
  core_seconds.reserve(static_cast<std::size_t>(options.repetitions));
  legacy_seconds.reserve(static_cast<std::size_t>(options.repetitions));
  Status core_status = Status::Success;

  auto run_core = [&] {
    std::copy(initial_arena.begin(), initial_arena.end(), batch.state().raw().begin());
    core_seconds.push_back(seconds([&] {
      for (int step = 0; step < options.steps; ++step) {
        core_status = stepper.step(batch, current_density, step * dt, dt);
        if (core_status != Status::Success)
          return;
      }
    }));
  };
  auto run_legacy = [&] {
    std::fill(legacy.begin(), legacy.end(), prototype);
    legacy_seconds.push_back(seconds([&] {
      for (int step = 0; step < options.steps; ++step)
        for (auto &cell : legacy)
          cell.timeStep_CC(dt);
    }));
  };

  for (int repetition = 0; repetition < options.repetitions; ++repetition) {
    if (repetition % 2 == 0) {
      run_core();
      if (core_status != Status::Success)
        break;
      run_legacy();
    } else {
      run_legacy();
      run_core();
      if (core_status != Status::Success)
        break;
    }
  }
  if (core_status != Status::Success) {
    std::fprintf(stderr, "PAY-1 core run failed with status %d\n", static_cast<int>(core_status));
    return 6;
  }

  const auto &layout = batch.layout();
  double max_state_error{};
  for (const auto legacy_domain : { neg, pos }) {
    const auto domain = legacy_domain == neg ? core::Domain::neg : core::Domain::pos;
    const auto d = core::domain_index(domain);
    for (int mode = 0; mode < static_cast<int>(settings::nch); ++mode) {
      max_state_error = std::max(
        max_state_error,
        std::abs(legacy.front().getStateObj().z(static_cast<std::size_t>(mode), legacy_domain)
                 - batch.state().at(layout.spm.z[d], mode, 0)));
      max_state_error = std::max(
        max_state_error,
        std::abs(legacy.back().getStateObj().z(static_cast<std::size_t>(mode), legacy_domain)
                 - batch.state().at(layout.spm.z[d], mode, options.lanes - 1)));
    }
  }
  max_state_error = std::max(
    max_state_error,
    std::abs(legacy.front().getStateObj().Wh()
             - batch.state().at(layout.energy_throughput, 0, 0)));
  const double voltage_error = std::abs(legacy.front().V() - stepper.terminalVoltage().front());
  const double core_median = median(core_seconds);
  const double legacy_median = median(legacy_seconds);
  const double speedup = legacy_median / core_median;
  const auto [core_min, core_max] = std::minmax_element(core_seconds.begin(), core_seconds.end());
  const auto [legacy_min, legacy_max] = std::minmax_element(legacy_seconds.begin(), legacy_seconds.end());
  const double conservative_speedup = *legacy_min / *core_max;
  const double checksum = legacy.front().getStateObj().z(0, neg)
                          + batch.state().at(layout.spm.z[core::domain_index(core::Domain::neg)], 0, 0)
                          + legacy.back().V() + stepper.terminalVoltage().back();

  std::printf("PAY-1 lanes=%d steps=%d dt=%.1f repetitions=%d cell_steps=%lld\n",
              options.lanes,
              options.steps,
              dt,
              options.repetitions,
              static_cast<long long>(options.lanes) * options.steps);
  std::printf("PAY-1 core_median_s=%.9g legacy_median_s=%.9g speedup=%.6gx\n",
              core_median,
              legacy_median,
              speedup);
  std::printf("PAY-1 core_range_s=[%.9g,%.9g] legacy_range_s=[%.9g,%.9g] conservative_speedup=%.6gx\n",
              *core_min,
              *core_max,
              *legacy_min,
              *legacy_max,
              conservative_speedup);
  std::printf("PAY-1 max_state_error=%.17g voltage_error=%.17g checksum=%.17g\n",
              max_state_error,
              voltage_error,
              checksum);
  std::printf("PAY-1 verdict=%s (abort<2x, target>=5x; timing trusted only on a quiet machine)\n",
              conservative_speedup < 2.0    ? "ABORT"
              : conservative_speedup >= 5.0 ? "TARGET"
                                            : "CONTINUE");

  if (!(max_state_error <= 1e-12 && voltage_error <= 1e-12 && std::isfinite(checksum)))
    return 7;
  return conservative_speedup < 2.0 ? 8 : 0;
}
