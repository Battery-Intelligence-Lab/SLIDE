/**
 * @file benchmark_PAY4_slide_pack.cpp
 * @brief Machine-readable SLIDE half of the PAY-4 liionpack comparison.
 */

#include "../src/core/PackStepper.hpp"
#include "../src/core/ParameterSet.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;
using slide::Status;
namespace core = slide::core;

struct Options
{
  int series{ 16 };
  int parallel{ 4 };
  int steps{ 60 };
  int repetitions{ 5 };
};

bool parsePositive(const char *text, int &value)
{
  char *end{};
  const auto parsed = std::strtol(text, &end, 10);
  if (end == text || *end != '\0' || parsed <= 0 || parsed > 1'000'000)
    return false;
  value = static_cast<int>(parsed);
  return true;
}

bool parseOptions(int argc, char **argv, Options &options)
{
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 >= argc)
      return false;
    const std::string_view key{ argv[i] };
    int *destination = key == "--series"        ? &options.series
                       : key == "--parallel"    ? &options.parallel
                       : key == "--steps"       ? &options.steps
                       : key == "--repetitions" ? &options.repetitions
                                                : nullptr;
    if (destination == nullptr || !parsePositive(argv[i + 1], *destination))
      return false;
  }
  return true;
}

double median(std::vector<double> values)
{
  std::ranges::sort(values);
  const auto middle = values.size() / 2;
  if (values.size() % 2 != 0)
    return values[middle];
  return 0.5 * (values[middle - 1] + values[middle]);
}

template <typename Function>
double seconds(Function &&function)
{
  const auto begin = clock_type::now();
  function();
  return std::chrono::duration<double>(clock_type::now() - begin).count();
}

} // namespace

int main(int argc, char **argv)
{
  Options options;
  if (!parseOptions(argc, argv, options)) {
    std::fprintf(stderr,
                 "usage: benchmark_PAY4_slide_pack [--series N] [--parallel N] "
                 "[--steps N] [--repetitions N]\n");
    return 2;
  }

  const auto lanes = options.series * options.parallel;
  core::SpmBatch batch;
  core::CompiledPackTopology topology;
  core::PackStepper stepper;
  Status status = Status::Success;
  const double setup_seconds = seconds([&] {
    core::ParameterSet parameters;
    core::SpmFactoryInput input;
    status = core::ParameterSet::chen2020(parameters);
    if (status == Status::Success)
      status = parameters.set("Contact resistance [Ohm]", 1e-4, "PAY-4");
    if (status == Status::Success)
      status = parameters.toSpmInput(input);
    input.initial_soc = 0.8;
    if (status == Status::Success)
      status = core::buildSpmBatch(input, { .nch = 12 }, lanes, batch);
    if (status == Status::Success)
      status = core::compilePackDescription(
        { .root = core::series(options.series,
                               core::parallel(options.parallel,
                                              core::cell({ .archetype = "spm" }))) },
        topology);
    std::array<core::SpmBatch *, 1> batches{ &batch };
    if (status == Status::Success)
      status = stepper.configure(topology, batches);
  });
  if (status != Status::Success)
    return 3;

  std::vector<double> checkpoint(stepper.checkpointSize());
  if (stepper.checkpoint(checkpoint) != Status::Success)
    return 4;

  std::vector<double> solve_times;
  solve_times.reserve(static_cast<std::size_t>(options.repetitions));
  const double pack_current = 5.0 * options.parallel;
  for (int repetition = 0; repetition < options.repetitions; ++repetition) {
    if (stepper.restore(checkpoint) != Status::Success)
      return 5;
    solve_times.push_back(seconds([&] {
      for (int step = 0; step < options.steps && status == Status::Success; ++step)
        status = stepper.stepExponential(pack_current,
                                         step * 10.0,
                                         10.0,
                                         {},
                                         core::PackSolveMode::ladder,
                                         1e-10);
    }));
    if (status != Status::Success)
      return 6;
  }

  const auto [minimum, maximum] = std::ranges::minmax(solve_times);
  std::printf(
    "{\"tool\":\"slide\",\"series\":%d,\"parallel\":%d,\"cells\":%d,"
    "\"steps\":%d,\"repetitions\":%d,\"setup_s\":%.17g,\"solve_median_s\":%.17g,"
    "\"solve_min_s\":%.17g,\"solve_max_s\":%.17g,\"terminal_voltage_V\":%.17g}\n",
    options.series,
    options.parallel,
    lanes,
    options.steps,
    options.repetitions,
    setup_seconds,
    median(solve_times),
    minimum,
    maximum,
    stepper.solution().terminal_voltage);
  return 0;
}
