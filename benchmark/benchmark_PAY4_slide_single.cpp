/**
 * @file benchmark_PAY4_slide_single.cpp
 * @brief Machine-readable SLIDE half of the PAY-4 PyBaMM comparison.
 */

#include "../src/core/Experiment.hpp"
#include "../src/core/ParameterSet.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;
using slide::Status;
namespace core = slide::core;

double median(std::vector<double> values)
{
  std::ranges::sort(values);
  const auto middle = values.size() / 2;
  return values.size() % 2 != 0 ? values[middle]
                                : 0.5 * (values[middle - 1] + values[middle]);
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
  int repetitions = 7;
  if (argc == 3 && std::string_view{ argv[1] } == "--repetitions") {
    char *end{};
    const long parsed = std::strtol(argv[2], &end, 10);
    if (end == argv[2] || *end != '\0' || parsed <= 0 || parsed > 1'000)
      return 2;
    repetitions = static_cast<int>(parsed);
  } else if (argc != 1) {
    std::fprintf(stderr, "usage: benchmark_PAY4_slide_single [--repetitions N]\n");
    return 2;
  }

  const std::vector<std::string> instructions{
    "Discharge at 1 C until 2.5 V",
    "Charge at 1 C until 4.2 V",
    "Hold at 4.2 V until C/20",
  };
  core::SpmBatch batch;
  core::Experiment experiment;
  core::CyclerV2 cycler;
  Status status = Status::Success;
  const double setup_seconds = seconds([&] {
    core::ParameterSet parameters;
    core::SpmFactoryInput input;
    core::ParseDiagnostic diagnostic;
    status = core::ParameterSet::chen2020(parameters);
    if (status == Status::Success)
      status = parameters.toSpmInput(input);
    input.initial_soc = 1.0;
    if (status == Status::Success)
      status = core::buildSpmBatch(input, { .nch = 12 }, 1, batch);
    if (status == Status::Success)
      status = core::Experiment::parse(instructions, experiment, diagnostic);
    if (status == Status::Success)
      status = cycler.configure(batch, core::CyclerIntegrator::exponential);
  });
  if (status != Status::Success)
    return 3;

  const std::vector<double> initial(batch.state().raw().begin(),
                                    batch.state().raw().end());
  std::vector<double> solve_times;
  solve_times.reserve(static_cast<std::size_t>(repetitions));
  core::ExperimentSolution solution;
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    std::ranges::copy(initial, batch.state().raw().begin());
    solve_times.push_back(seconds([&] {
      status = cycler.run(experiment, 10.0, solution);
    }));
    if (status != Status::Success)
      return 4;
  }

  const auto [minimum, maximum] = std::ranges::minmax(solve_times);
  std::printf(
    "{\"tool\":\"slide\",\"case\":\"single_spm_1C_discharge_cccv\","
    "\"repetitions\":%d,\"setup_s\":%.17g,\"solve_median_s\":%.17g,"
    "\"solve_min_s\":%.17g,\"solve_max_s\":%.17g,\"samples\":%zu,"
    "\"duration_s\":%.17g,\"terminal_voltage_V\":%.17g}\n",
    repetitions,
    setup_seconds,
    median(solve_times),
    minimum,
    maximum,
    solution.time.size(),
    solution.time.empty() ? 0.0 : solution.time.back(),
    solution.voltage.empty() ? 0.0 : solution.voltage.back());
  return 0;
}
