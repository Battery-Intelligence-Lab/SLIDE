/**
 * @file benchmark_PAY5_cuda.cpp
 * @brief P8 GPU payoff: 100k heterogeneous lanes, exact-modal CPU versus CUDA.
 */

#include "../src/core/CudaSpmBatch.hpp"
#include "../src/core/ExponentialModal.hpp"
#include "../src/core/ParameterSet.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <span>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;
namespace core = slide::core;

template <class Function>
double seconds(Function &&function)
{
  const auto begin = clock_type::now();
  function();
  return std::chrono::duration<double>(clock_type::now() - begin).count();
}

double median(std::vector<double> values)
{
  std::ranges::sort(values);
  return values[values.size() / 2];
}

void applySpread(core::SpmBatch &batch, const core::SpmFactoryInput &input)
{
  auto &state = batch.state();
  const auto &layout = batch.layout().spm;
  for (int lane = 0; lane < batch.n_lanes(); ++lane) {
    const double fraction = static_cast<double>(lane)
                            / static_cast<double>(batch.n_lanes() - 1);
    // A symmetric 0.5C-throughput discharge/charge keeps every benchmark
    // lane in the valid OCV domain while D and contact resistance remain
    // heterogeneous across the full 100k ensemble.
    const double soc = 0.65 + 0.25 * fraction;
    for (const core::Domain domain : core::domains) {
      const auto d = core::domain_index(domain);
      const auto &material = input.design.electrode[d].active_material;
      const double base = material.x_0
                          + input.initial_soc * (material.x_100 - material.x_0);
      const double varied = material.x_0
                            + soc * (material.x_100 - material.x_0);
      for (int mode = 0; mode < layout.z[d].rows; ++mode)
        state.at(layout.z[d], mode, lane) *= varied / base;
      const double diffusion_scale = 0.95 + 0.10
                                             * static_cast<double>((lane * 17
                                                                    + static_cast<int>(d) * 13)
                                                                   % 101)
                                             / 100.0;
      state.at(layout.diffusion_coefficient[d], 0, lane) =
        material.D_s.reference_value * diffusion_scale;
    }
    const double resistance_scale = 0.95 + 0.10
                                            * static_cast<double>((lane * 29) % 97)
                                            / 96.0;
    state.at(layout.current_collector_resistance, 0, lane) =
      input.initial_current_collector_resistance * resistance_scale;
  }
}

} // namespace

int main()
{
  constexpr int lanes = 100'000;
  constexpr int steps = 360;
  constexpr int repetitions = 3;
  constexpr double dt = 10.0;
  slide::Status status = slide::Status::Success;
  core::ParameterSet parameters;
  core::SpmFactoryInput input;
  status = core::ParameterSet::chen2020(parameters);
  if (status == slide::Status::Success)
    status = parameters.toSpmInput(input);
  input.initial_soc = 0.55;
  const core::SpmModelOptions options{ .nch = 12 };

  core::SpmBatch cpu;
  const double cpu_setup = seconds([&] {
    if (status == slide::Status::Success)
      status = core::buildSpmBatch(input, options, lanes, cpu);
  });
  if (status != slide::Status::Success)
    return 2;
  applySpread(cpu, input);

  core::CudaSpmBatch gpu;
  const double gpu_setup = seconds([&] {
    status = gpu.build(input, options, lanes);
  });
  if (status != slide::Status::Success)
    return 3;
  applySpread(gpu.hostBatch(), input);
  const double h2d = seconds([&] {
    status = gpu.uploadState();
    if (status == slide::Status::Success)
      status = gpu.synchronize();
  });
  if (status != slide::Status::Success)
    return 4;

  const double current = input.design.capacity_Ah;
  std::vector<double> discharge_density(static_cast<std::size_t>(lanes),
                                        current / input.design.electrode_area);
  std::vector<double> charge_density(static_cast<std::size_t>(lanes),
                                     -current / input.design.electrode_area);
  std::vector<double> initial(cpu.state().raw().begin(), cpu.state().raw().end());
  core::ExponentialModal cpu_stepper{ cpu };

  // Warm both math/runtime paths, then restore identical bytes before timing.
  status = cpu_stepper.step(cpu, discharge_density, 0.0, dt);
  std::memcpy(cpu.state().raw().data(), initial.data(),
              initial.size() * sizeof(double));
  if (status == slide::Status::Success)
    status = gpu.checkpoint();
  if (status == slide::Status::Success)
    status = gpu.step(discharge_density, 0.0, dt);
  if (status == slide::Status::Success)
    status = gpu.synchronize();
  if (status == slide::Status::Success)
    status = gpu.restore();
  if (status == slide::Status::Success)
    status = gpu.synchronize();
  if (status != slide::Status::Success)
    return 5;

  std::vector<double> cpu_times, gpu_times;
  cpu_times.reserve(repetitions);
  gpu_times.reserve(repetitions);
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    auto run_cpu = [&] {
      std::memcpy(cpu.state().raw().data(), initial.data(),
                  initial.size() * sizeof(double));
      const double elapsed = seconds([&] {
        for (int step = 0; step < steps && status == slide::Status::Success; ++step) {
          const auto &density = step < steps / 2 ? discharge_density
                                                 : charge_density;
          status = cpu_stepper.step(cpu, density, step * dt, dt);
        }
      });
      cpu_times.push_back(elapsed);
    };
    auto run_gpu = [&] {
      status = gpu.restore();
      if (status == slide::Status::Success)
        status = gpu.synchronize();
      const double elapsed = seconds([&] {
        for (int step = 0; step < steps && status == slide::Status::Success; ++step) {
          const auto &density = step < steps / 2 ? discharge_density
                                                 : charge_density;
          status = gpu.step(density, step * dt, dt);
        }
        if (status == slide::Status::Success)
          status = gpu.synchronize();
      });
      gpu_times.push_back(elapsed);
    };
    if (repetition % 2 == 0) {
      run_gpu();
      run_cpu();
    } else {
      run_cpu();
      run_gpu();
    }
    if (status != slide::Status::Success)
      return 6;
  }

  // Make both final states correspond to a full run after the last restore.
  std::memcpy(cpu.state().raw().data(), initial.data(),
              initial.size() * sizeof(double));
  for (int step = 0; step < steps && status == slide::Status::Success; ++step) {
    const auto &density = step < steps / 2 ? discharge_density
                                           : charge_density;
    status = cpu_stepper.step(cpu, density, step * dt, dt);
  }
  const double d2h = seconds([&] { status = gpu.downloadState(); });
  if (status != slide::Status::Success)
    return 7;
  double max_state_error{};
  double max_scaled_state_error{};
  for (std::size_t i = 0; i < cpu.state().size(); ++i) {
    const double error = std::abs(cpu.state().raw()[i]
                                  - gpu.hostBatch().state().raw()[i]);
    max_state_error = std::max(max_state_error, error);
    max_scaled_state_error = std::max(
      max_scaled_state_error,
      error / (1e-12 + 2e-10 * std::max(std::abs(cpu.state().raw()[i]),
                                        std::abs(gpu.hostBatch().state().raw()[i]))));
  }
  double max_voltage_error{};
  for (int lane = 0; lane < lanes; ++lane)
    max_voltage_error = std::max(
      max_voltage_error,
      std::abs(cpu_stepper.terminalVoltage()[static_cast<std::size_t>(lane)]
               - gpu.terminalVoltage()[static_cast<std::size_t>(lane)]));

  const double cpu_median = median(cpu_times);
  const double gpu_median = median(gpu_times);
  const double speedup = cpu_median / gpu_median;
  std::printf(
    "{\"cells\":%d,\"steps\":%d,\"repetitions\":%d,"
    "\"cpu_setup_s\":%.17g,\"gpu_setup_s\":%.17g,"
    "\"h2d_s\":%.17g,\"cpu_solve_median_s\":%.17g,"
    "\"gpu_solve_median_s\":%.17g,\"d2h_s\":%.17g,"
    "\"speedup\":%.17g,\"max_state_error\":%.17g,"
    "\"max_scaled_state_error\":%.17g,\"max_voltage_error_V\":%.17g,"
    "\"device_bytes\":%zu,\"device_allocations\":%zu,"
    "\"device_wide_syncs\":%zu}\n",
    lanes, steps, repetitions, cpu_setup, gpu_setup, h2d, cpu_median,
    gpu_median, d2h, speedup, max_state_error, max_scaled_state_error,
    max_voltage_error, gpu.deviceArenaBytes(), gpu.deviceAllocationCount(),
    gpu.deviceWideSynchronizationCount());
  return max_scaled_state_error <= 1.0 && max_voltage_error <= 2e-6
             && speedup >= 2.0
           ? 0
           : 8;
}
