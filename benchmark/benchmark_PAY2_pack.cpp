/**
 * @file benchmark_PAY2_pack.cpp
 * @brief PAY-2: heterogeneous 16s4p compiled pack versus legacy nested modules.
 */

#include "../src/core/PackStepper.hpp"
#include "../tests/support/KokamSpmFixture.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string_view>
#include <vector>

using namespace slide;

namespace {

struct Options
{
  int steps{ 1'800 };
  int repetitions{ 3 };
};

bool parsePositive(std::string_view text, int &value)
{
  int parsed{};
  const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size() || parsed <= 0)
    return false;
  value = parsed;
  return true;
}

bool parseOptions(int argc, char **argv, Options &options)
{
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 >= argc)
      return false;
    const std::string_view name{ argv[i] };
    if (name == "--steps") {
      if (!parsePositive(argv[i + 1], options.steps))
        return false;
    } else if (name == "--repetitions") {
      if (!parsePositive(argv[i + 1], options.repetitions))
        return false;
    } else {
      return false;
    }
  }
  return true;
}

template <class Function>
double seconds(Function &&function)
{
  const auto start = std::chrono::steady_clock::now();
  function();
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  const auto middle = values.size() / 2;
  return values.size() % 2 == 0 ? 0.5 * (values[middle - 1] + values[middle])
                                : values[middle];
}

constexpr int series_count = 16;
constexpr int parallel_count = 4;
constexpr int lanes = series_count * parallel_count;
constexpr double branch_current = 16.0;
constexpr double pack_current = branch_current * parallel_count;
constexpr double dt = 1.0;

double laneSoc(int lane)
{
  const int parallel_lane = lane % parallel_count;
  return 0.90 + 0.03 * (parallel_lane - 1.5);
}

struct LegacyPack
{
  Deep_ptr<Module_s> root{};
  std::array<Cell_SPM *, lanes> cells{};
};

LegacyPack makeLegacyPack()
{
  std::array<Deep_ptr<StorageUnit>, series_count> modules;
  for (int s = 0; s < series_count; ++s) {
    std::array<Deep_ptr<StorageUnit>, parallel_count> cells;
    for (int p = 0; p < parallel_count; ++p) {
      auto cell = make<Cell_SPM>();
      if (test_support::initialize_legacy_kokam(
            *cell, laneSoc(s * parallel_count + p), 0.0)
          != Status::Success)
        return {};
      cells[static_cast<std::size_t>(p)] = std::move(cell);
    }
    auto module = make<Module_p>("p", 288.0, false, false, parallel_count, 1, 1);
    module->setSUs(cells, false, false);
    modules[static_cast<std::size_t>(s)] = std::move(module);
  }
  LegacyPack pack;
  pack.root = make<Module_s>("s", 288.0, false, false, lanes, 1, 1);
  pack.root->setSUs(modules, false, false);
  pack.root->setBlockDegAndTherm(true);
  for (int s = 0; s < series_count; ++s) {
    auto *module = dynamic_cast<Module_p *>(pack.root->getSUs()[static_cast<std::size_t>(s)].get());
    for (int p = 0; p < parallel_count; ++p)
      pack.cells[static_cast<std::size_t>(s * parallel_count + p)] =
        dynamic_cast<Cell_SPM *>(module->getSUs()[static_cast<std::size_t>(p)].get());
  }
  return pack;
}

bool initializeCore(core::SpmBatch &batch)
{
  Cell_SPM prototype;
  double environment{}, reference{};
  prototype.getTemperatures(&environment, &reference);
  const auto base = test_support::make_legacy_kokam_input(0.9, prototype.T(), reference);
  if (core::buildSpmBatch(base, {}, lanes, batch) != Status::Success)
    return false;
  for (int lane = 0; lane < lanes; ++lane) {
    core::SpmBatch single;
    auto input = base;
    input.initial_soc = laneSoc(lane);
    if (core::buildSpmBatch(input, {}, 1, single) != Status::Success)
      return false;
    for (int row = 0; row < batch.state().n_rows(); ++row)
      batch.state().row(row)[static_cast<std::size_t>(lane)] = single.state().row(row).front();
  }
  return true;
}

} // namespace

int main(int argc, char **argv)
{
  Options options;
  if (!parseOptions(argc, argv, options)) {
    std::fprintf(stderr,
                 "usage: benchmark_PAY2_pack [--steps N] [--repetitions N]\n");
    return 2;
  }

  core::SpmBatch batch;
  if (!initializeCore(batch))
    return 3;
  core::CompiledPackTopology topology;
  if (core::compilePackDescription(
        { .root = core::series(series_count,
                               core::parallel(parallel_count,
                                              core::cell({ .archetype = "spm" }))) },
        topology)
      != Status::Success)
    return 4;
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper compiled;
  if (compiled.configure(topology, batches) != Status::Success)
    return 5;
  std::vector<double> initial(compiled.checkpointSize());
  if (compiled.checkpoint(initial) != Status::Success)
    return 6;

  std::vector<double> core_times, legacy_times;
  core_times.reserve(static_cast<std::size_t>(options.repetitions));
  legacy_times.reserve(static_cast<std::size_t>(options.repetitions));
  Status core_status = Status::Success;
  Status legacy_status = Status::Success;
  LegacyPack final_legacy;

  auto runCore = [&] {
    core_status = compiled.restore(initial);
    core_times.push_back(seconds([&] {
      for (int step = 0; step < options.steps && core_status == Status::Success;) {
        const int substeps = std::min(10, options.steps - step);
        core_status = compiled.step(pack_current, step * dt, dt, {}, core::PackSolveMode::ladder, 1e-6, substeps);
        step += substeps;
      }
    }));
  };
  auto runLegacy = [&] {
    auto legacy = makeLegacyPack();
    if (!legacy.root) {
      legacy_status = Status::Invalid_states;
      return;
    }
    legacy_times.push_back(seconds([&] {
      for (int step = 0; step < options.steps && legacy_status == Status::Success;) {
        const int substeps = std::min(10, options.steps - step);
        legacy_status = legacy.root->setCurrent(pack_current, false, false);
        if (legacy_status == Status::Success)
          legacy.root->timeStep_CC(dt, substeps);
        step += substeps;
      }
    }));
    final_legacy = std::move(legacy);
  };

  for (int repetition = 0; repetition < options.repetitions; ++repetition) {
    if (repetition % 2 == 0) {
      runCore();
      runLegacy();
    } else {
      runLegacy();
      runCore();
    }
    if (core_status != Status::Success || legacy_status != Status::Success)
      break;
  }
  if (core_status != Status::Success || legacy_status != Status::Success)
    return 7;

  if (compiled.solveElectrical(pack_current, core::PackSolveMode::ladder) != Status::Success
      || final_legacy.root->setCurrent(pack_current, false, false) != Status::Success)
    return 8;
  double max_current_error{};
  double max_state_error{};
  const auto &layout = batch.layout();
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    max_current_error = std::max(max_current_error,
                                 std::abs(compiled.solution().cell_current[i]
                                          - final_legacy.cells[i]->I()));
    for (const auto legacy_domain : { neg, pos }) {
      const auto domain = legacy_domain == neg ? core::Domain::neg : core::Domain::pos;
      const auto d = core::domain_index(domain);
      for (int mode = 0; mode < static_cast<int>(settings::nch); ++mode)
        max_state_error = std::max(
          max_state_error,
          std::abs(batch.state().at(layout.spm.z[d], mode, lane)
                   - final_legacy.cells[i]->getStateObj().z(
                     static_cast<std::size_t>(mode), legacy_domain)));
    }
  }
  const double voltage_error = std::abs(compiled.solution().terminal_voltage
                                        - final_legacy.root->V());
  const double core_median = median(core_times);
  const double legacy_median = median(legacy_times);
  const double speedup = legacy_median / core_median;
  const auto [core_min, core_max] = std::minmax_element(core_times.begin(), core_times.end());
  const auto [legacy_min, legacy_max] = std::minmax_element(legacy_times.begin(), legacy_times.end());
  const double conservative = *legacy_min / *core_max;

  std::printf("PAY-2 pack=16s4p steps=%d repetitions=%d cell_steps=%lld\n",
              options.steps,
              options.repetitions,
              static_cast<long long>(lanes) * options.steps);
  std::printf("PAY-2 core_median_s=%.9g legacy_median_s=%.9g speedup=%.6gx\n",
              core_median,
              legacy_median,
              speedup);
  std::printf("PAY-2 core_range_s=[%.9g,%.9g] legacy_range_s=[%.9g,%.9g] conservative_speedup=%.6gx\n",
              *core_min,
              *core_max,
              *legacy_min,
              *legacy_max,
              conservative);
  std::printf("PAY-2 max_state_error=%.17g current_error=%.17g voltage_error=%.17g\n",
              max_state_error,
              max_current_error,
              voltage_error);
  std::printf("PAY-2 verdict=%s (abort<3x, target>=10x; timing trusted only on a quiet machine)\n",
              conservative < 3.0     ? "ABORT"
              : conservative >= 10.0 ? "TARGET"
                                     : "CONTINUE");

  if (!(max_state_error <= 5e-8 && max_current_error <= 1e-5
        && voltage_error <= 1e-6))
    return 9;
  return conservative < 3.0 ? 10 : 0;
}
