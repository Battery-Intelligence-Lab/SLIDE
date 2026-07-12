/**
 * @file core_CudaSpmBatch_test.cpp
 * @brief P8-G2 CUDA correctness, host-coupling, conservation, and rollback gates.
 */

#include "../../src/core/CudaSpmBatch.hpp"
#include "../../src/core/ExponentialModal.hpp"
#include "../../src/core/PackSolver.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../../src/core/SpectralModel.hpp"
#include "../support/RecordedBits.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <span>
#include <vector>

using namespace slide;

namespace {

core::SpmFactoryInput chen2020(double soc = 0.55)
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  input.initial_soc = soc;
  return input;
}

void applySpread(core::SpmBatch &batch,
                 const core::SpmFactoryInput &input,
                 double minimum_soc,
                 double maximum_soc)
{
  auto &state = batch.state();
  const auto &layout = batch.layout().spm;
  for (int lane = 0; lane < batch.n_lanes(); ++lane) {
    const double fraction = batch.n_lanes() == 1
                              ? 0.0
                              : static_cast<double>(lane)
                                  / static_cast<double>(batch.n_lanes() - 1);
    const double soc = minimum_soc + fraction * (maximum_soc - minimum_soc);
    for (const core::Domain domain : core::domains) {
      const auto d = core::domain_index(domain);
      const auto &material = input.design.electrode[d].active_material;
      const double base_stoichiometry = material.x_0
                                        + input.initial_soc
                                            * (material.x_100 - material.x_0);
      const double varied_stoichiometry = material.x_0
                                          + soc * (material.x_100 - material.x_0);
      const double scale = varied_stoichiometry / base_stoichiometry;
      for (int mode = 0; mode < layout.z[d].rows; ++mode)
        state.at(layout.z[d], mode, lane) *= scale;
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

double maximumStateError(std::span<const double> expected,
                         std::span<const double> actual)
{
  REQUIRE(expected.size() == actual.size());
  double worst{};
  for (std::size_t i = 0; i < expected.size(); ++i)
    worst = std::max(worst, std::abs(expected[i] - actual[i]));
  return worst;
}

void requireStateBand(std::span<const double> cpu,
                      std::span<const double> gpu)
{
  REQUIRE(cpu.size() == gpu.size());
  for (std::size_t i = 0; i < cpu.size(); ++i) {
    const double tolerance = 1e-12
                             + 2e-10 * std::max(std::abs(cpu[i]),
                                                std::abs(gpu[i]));
    CAPTURE(i, cpu[i], gpu[i], tolerance);
    REQUIRE(std::abs(cpu[i] - gpu[i]) <= tolerance);
  }
}

} // namespace

TEST_CASE("PC-10 CUDA batch retains its backend-local pre-refactor scalar-kernel bits",
          "[core][cuda][PC-10][recorded]")
{
  REQUIRE(core::CudaSpmBatch::available());
  constexpr int lanes = 4;
  auto input = chen2020();
  core::CudaSpmBatch gpu;
  REQUIRE(gpu.build(input, { .nch = 12 }, lanes) == Status::Success);
  applySpread(gpu.hostBatch(), input, 0.31, 0.79);

  test_support::RecordedBits recorded;
  recorded.append(std::span<const double>{ gpu.hostBatch().state().raw() });
  REQUIRE(gpu.uploadState() == Status::Success);
  const double one_c_density = input.design.capacity_Ah
                               / input.design.electrode_area;
  const std::array density{ -0.45 * one_c_density,
                            0.0,
                            0.35 * one_c_density,
                            0.85 * one_c_density };
  double time{};
  for (const double dt : { 1e-5, 7.25, 19.0 }) {
    REQUIRE(gpu.step(density, time, dt) == Status::Success);
    time += dt;
    REQUIRE(gpu.downloadState() == Status::Success);
    recorded.append(std::span<const double>{ gpu.hostBatch().state().raw() });
    recorded.append(gpu.terminalVoltage());
  }
  CAPTURE(recorded.values, recorded.fnv1a, recorded.mixed);
  REQUIRE(recorded.values == 1388);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
  CHECK(recorded.fnv1a == UINT64_C(0xa59c34d1685ddb19));
  CHECK(recorded.mixed == UINT64_C(0x2e88cc95b54eb713));
#endif
}

TEST_CASE("P8-G2 CUDA matches CPU for 10003 heterogeneous lanes and rolls back",
          "[core][cuda][P8-G2]")
{
  REQUIRE(core::CudaSpmBatch::available());
  constexpr int lanes = 10'003;
  constexpr int nch = 12;
  constexpr int steps = 60;
  constexpr double dt = 10.0;
  auto input = chen2020();
  const core::SpmModelOptions options{ .nch = nch };
  core::SpmBatch cpu;
  REQUIRE(core::buildSpmBatch(input, options, lanes, cpu) == Status::Success);
  core::CudaSpmBatch gpu;
  REQUIRE(gpu.build(input, options, lanes) == Status::Success);
  applySpread(cpu, input, 0.2, 0.9);
  applySpread(gpu.hostBatch(), input, 0.2, 0.9);
  REQUIRE(gpu.uploadState() == Status::Success);
  REQUIRE(gpu.deviceAllocationCount() == 1);
  REQUIRE(gpu.deviceWideSynchronizationCount() == 0);

  const double current = input.design.capacity_Ah;
  std::vector<double> density(static_cast<std::size_t>(lanes),
                              current / input.design.electrode_area);
  core::ExponentialModal cpu_stepper{ cpu };
  for (int step = 0; step < steps; ++step) {
    REQUIRE(cpu_stepper.step(cpu, density, step * dt, dt) == Status::Success);
    REQUIRE(gpu.step(density, step * dt, dt) == Status::Success);
  }
  REQUIRE(gpu.synchronize() == Status::Success);
  REQUIRE(gpu.downloadState() == Status::Success);
  requireStateBand(cpu.state().raw(), gpu.hostBatch().state().raw());
  double worst_voltage{};
  for (int lane = 0; lane < lanes; ++lane)
    worst_voltage = std::max(
      worst_voltage,
      std::abs(cpu_stepper.terminalVoltage()[static_cast<std::size_t>(lane)]
               - gpu.terminalVoltage()[static_cast<std::size_t>(lane)]));
  CHECK(worst_voltage <= 2e-6);

  core::PerDomain<double> radius{};
  for (const auto domain : core::domains)
    radius[core::domain_index(domain)] =
      input.design.electrode[core::domain_index(domain)].particle_radius;
  core::CompiledSpectralModel<nch> spectral;
  REQUIRE(core::compileSpectralModel<nch>(radius, spectral) == Status::Success);
  const auto neg = core::domain_index(core::Domain::neg);
  const auto &negative = input.design.electrode[neg];
  const double specific_area = 3.0 * negative.active_fraction
                               / negative.particle_radius;
  core::SpmBatch initial;
  REQUIRE(core::buildSpmBatch(input, options, lanes, initial) == Status::Success);
  applySpread(initial, input, 0.2, 0.9);
  double worst_inventory{};
  for (int lane = 0; lane < lanes; ++lane) {
    const double initial_mode = initial.state().at(
      initial.layout().spm.z[neg], spectral.zero_mode[neg], lane);
    const double final_mode = gpu.hostBatch().state().at(
      gpu.hostBatch().layout().spm.z[neg], spectral.zero_mode[neg], lane);
    const double stored = (final_mode - initial_mode)
                          * input.design.electrode_area * specific_area
                          * 96487.0 * negative.thickness
                          / (3600.0
                             * spectral.B[neg][static_cast<std::size_t>(
                               spectral.zero_mode[neg])]);
    worst_inventory = std::max(
      worst_inventory,
      std::abs(current * steps * dt / 3600.0 - stored));
  }
  CHECK(worst_inventory <= 1e-9);

  std::vector<double> accepted(gpu.hostBatch().state().raw().begin(),
                               gpu.hostBatch().state().raw().end());
  REQUIRE(gpu.checkpoint() == Status::Success);
  REQUIRE(gpu.step(density, steps * dt, dt) == Status::Success);
  REQUIRE(gpu.restore() == Status::Success);
  REQUIRE(gpu.downloadState() == Status::Success);
  CHECK(std::memcmp(accepted.data(),
                    gpu.hostBatch().state().raw().data(),
                    accepted.size() * sizeof(double))
        == 0);
  CHECK(gpu.deviceWideSynchronizationCount() == 0);
}

TEST_CASE("P8-G2 host PackSolver coupling stays inside voltage/current bands",
          "[core][cuda][pack][P8-G2]")
{
  constexpr int series = 16;
  constexpr int parallel = 4;
  constexpr int lanes = series * parallel;
  constexpr int steps = 60;
  constexpr double dt = 10.0;
  auto input = chen2020(0.8);
  const core::SpmModelOptions options{ .nch = 12 };
  core::SpmBatch cpu;
  REQUIRE(core::buildSpmBatch(input, options, lanes, cpu) == Status::Success);
  core::CudaSpmBatch gpu;
  REQUIRE(gpu.build(input, options, lanes) == Status::Success);
  applySpread(cpu, input, 0.79, 0.81);
  applySpread(gpu.hostBatch(), input, 0.79, 0.81);
  REQUIRE(gpu.uploadState() == Status::Success);

  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(series,
                                   core::parallel(parallel,
                                                  core::cell({ .archetype = "spm" }))) },
            topology)
          == Status::Success);
  auto cpu_view = core::TheveninBatchView::bind(cpu, lanes);
  auto gpu_view = core::TheveninBatchView::bind(gpu.hostBatch(), lanes);
  core::PackSolver cpu_solver, gpu_solver;
  REQUIRE(cpu_solver.configure(topology, std::span{ &cpu_view, 1 })
          == Status::Success);
  REQUIRE(gpu_solver.configure(topology, std::span{ &gpu_view, 1 })
          == Status::Success);
  core::ExponentialModal cpu_stepper{ cpu };
  std::vector<double> cpu_density(lanes), gpu_density(lanes);
  double worst_pack_voltage{};
  double worst_branch_current{};
  const double pack_current = input.design.capacity_Ah * parallel;
  for (int step = 0; step < steps; ++step) {
    REQUIRE(cpu_solver.solve(pack_current, core::PackSolveMode::ladder)
            == Status::Success);
    REQUIRE(gpu_solver.solve(pack_current, core::PackSolveMode::ladder)
            == Status::Success);
    worst_pack_voltage = std::max(
      worst_pack_voltage,
      std::abs(cpu_solver.solution().terminal_voltage
               - gpu_solver.solution().terminal_voltage));
    for (int lane = 0; lane < lanes; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      worst_branch_current = std::max(
        worst_branch_current,
        std::abs(cpu_solver.solution().cell_current[i]
                 - gpu_solver.solution().cell_current[i]));
      cpu_density[i] = cpu_solver.solution().cell_current[i]
                       / input.design.electrode_area;
      gpu_density[i] = gpu_solver.solution().cell_current[i]
                       / input.design.electrode_area;
    }
    REQUIRE(cpu_stepper.step(cpu, cpu_density, step * dt, dt)
            == Status::Success);
    REQUIRE(gpu.step(gpu_density, step * dt, dt) == Status::Success);
    REQUIRE(gpu.synchronize() == Status::Success);
    REQUIRE(gpu.downloadState() == Status::Success);
  }
  CAPTURE(worst_pack_voltage, worst_branch_current,
          maximumStateError(cpu.state().raw(), gpu.hostBatch().state().raw()));
  CHECK(worst_pack_voltage <= 2e-6);
  CHECK(worst_branch_current <= 2e-6);
  requireStateBand(cpu.state().raw(), gpu.hostBatch().state().raw());
}

TEST_CASE("P8-G3 CUDA pinned side-stream snapshots decode bitwise",
          "[core][cuda][recorder][P8-G3]")
{
  auto input = chen2020(0.7);
  core::CudaSpmBatch gpu;
  REQUIRE(gpu.build(input, { .nch = 12 }, 2) == Status::Success);
  const auto path = std::filesystem::temp_directory_path()
                    / "slide_cuda_async.slcmp";
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  core::CudaAsyncRecorder recorder;
  REQUIRE(recorder.configure(
            gpu,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  CHECK(recorder.usesPinnedMemory());
  CHECK(recorder.usesNonDefaultStream());

  const std::array<double, 2> density{
    input.design.capacity_Ah / input.design.electrode_area,
    -0.5 * input.design.capacity_Ah / input.design.electrode_area
  };
  constexpr int snapshots = 6;
  std::vector<double> expected;
  expected.reserve(static_cast<std::size_t>(snapshots)
                   * gpu.hostBatch().state().size());
  for (int step = 0; step < snapshots; ++step) {
    REQUIRE(gpu.step(density, step * 10.0, 10.0) == Status::Success);
    REQUIRE(recorder.enqueue(static_cast<std::uint64_t>(step + 1))
            == Status::Success);
    REQUIRE(gpu.downloadState() == Status::Success);
    expected.insert(expected.end(),
                    gpu.hostBatch().state().raw().begin(),
                    gpu.hostBatch().state().raw().end());
  }
  REQUIRE(recorder.finish() == Status::Success);
  CHECK(recorder.snapshotsWritten() == snapshots);
  CHECK(recorder.thinnedSnapshots() == 0);
  CHECK(gpu.deviceWideSynchronizationCount() == 0);

  core::CompressedRecording decoded;
  REQUIRE(decoded.open(path) == Status::Success);
  REQUIRE(decoded.size() == snapshots);
  const auto state_values = gpu.hostBatch().state().size();
  for (int snapshot = 0; snapshot < snapshots; ++snapshot) {
    const auto actual = decoded.snapshot(static_cast<std::size_t>(snapshot));
    CHECK(actual.accepted_step == static_cast<std::uint64_t>(snapshot + 1));
    CHECK(std::equal(actual.current_density.begin(),
                     actual.current_density.end(),
                     density.begin(),
                     density.end()));
    const auto offset = static_cast<std::size_t>(snapshot) * state_values;
    CHECK(std::equal(actual.state.begin(),
                     actual.state.end(),
                     expected.begin() + static_cast<std::ptrdiff_t>(offset),
                     expected.begin()
                       + static_cast<std::ptrdiff_t>(offset + state_values)));
  }
  std::filesystem::remove(path, ignored);
}
