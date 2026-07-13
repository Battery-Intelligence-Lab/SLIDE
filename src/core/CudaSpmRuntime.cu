/**
 * @file CudaSpmRuntime.cu
 * @brief One-thread-per-lane fused exact-modal CUDA kernel and runtime ownership.
 */

#include "CudaSpmData.hpp"
#include "SpmScalarKernels.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <vector>

namespace slide::core::cuda_detail {
namespace {

constexpr std::size_t alignment = 256;

std::size_t alignUp(std::size_t value)
{
  return (value + alignment - 1) & ~(alignment - 1);
}

__device__ double curveEval(const double *storage, CurveRef curve, double query)
{
  const double *x = storage + curve.x_offset;
  const double *y = storage + curve.y_offset;
  if (query <= x[0])
    return y[0];
  if (query >= x[curve.count - 1])
    return y[curve.count - 1];
  int left = 0;
  int right = curve.count - 1;
  while (right - left > 1) {
    const int middle = left + (right - left) / 2;
    if (query >= x[middle])
      left = middle;
    else
      right = middle;
  }
  return spm_scalar::linearInterpolate(
    query, x[left], x[left + 1], y[left], y[left + 1]);
}

__device__ double &at(double *state, int stride, int row, int lane)
{
  return state[static_cast<std::size_t>(row) * stride + lane];
}

__device__ double get(const double *state, int stride, int row, int lane)
{
  return state[static_cast<std::size_t>(row) * stride + lane];
}

__global__ void exactModalKernel(KernelParams params,
                                 double *state,
                                 const double *current_density,
                                 double *terminal_voltage,
                                 const double *curves,
                                 int *status,
                                 double dt)
{
  const int lane = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (lane >= params.lanes)
    return;

  const double temperature = get(state, params.stride, params.temperature_row, lane);
  double effective_diffusivity[2]{};
  double molar_flux[2]{};
  for (int domain = 0; domain < 2; ++domain) {
    const double arrhenius = spm_scalar::arrheniusFactor(
      params.reference_temperature, temperature, params.Rg);
    effective_diffusivity[domain] = spm_scalar::activatedValue(
      get(state, params.stride, params.diffusion_row[domain], lane),
      params.diffusion_activation[domain],
      arrhenius);
    const double denominator = spm_scalar::fluxDenominator(
      get(state, params.stride, params.specific_area_row[domain], lane),
      params.n,
      params.F,
      get(state, params.stride, params.thickness_row[domain], lane));
    const double sign = domain == 0 ? 1.0 : -1.0;
    molar_flux[domain] = spm_scalar::molarFlux(
      sign, current_density[lane], denominator);

    for (int mode = 0; mode < params.nch; ++mode) {
      double &z = at(state, params.stride, params.z_row[domain] + mode, lane);
      SLIDE_SPM_ADVANCE_MODAL_CUDA(z,
                                   effective_diffusivity[domain],
                                   params.A[domain][mode],
                                   dt,
                                   params.B[domain][mode],
                                   molar_flux[domain]);
    }
  }

  double surface_stoichiometry[2]{};
  double overpotential[2]{};
  double electrode_ocv[2]{};
  const double arrhenius = spm_scalar::arrheniusFactor(
    params.reference_temperature, temperature, params.Rg);
  for (int domain = 0; domain < 2; ++domain) {
    double concentration{};
    for (int mode = 0; mode < params.nch; ++mode)
      concentration += params.surface_C[domain][mode]
                       * get(state, params.stride,
                             params.z_row[domain] + mode, lane);
    concentration = spm_scalar::concentrationOutput(
      concentration,
      params.surface_D[domain],
      molar_flux[domain],
      effective_diffusivity[domain]);
    surface_stoichiometry[domain] = spm_scalar::surfaceStoichiometry(
      concentration, params.electrode[domain].cs_max);
    if (!(surface_stoichiometry[domain] > 0.0
          && surface_stoichiometry[domain] < 1.0)
        || !isfinite(surface_stoichiometry[domain])) {
      atomicExch(status, 1);
      terminal_voltage[lane] = nan("");
      return;
    }
    const double reaction_rate = spm_scalar::activatedValue(
      params.electrode[domain].reaction_rate_ref,
      params.electrode[domain].reaction_activation,
      arrhenius);
    const double exchange = spm_scalar::exchangeCurrent(
      reaction_rate,
      params.n,
      params.F,
      params.electrolyte_concentration,
      concentration,
      params.electrode[domain].cs_max);
    const double area = get(state, params.stride,
                            params.specific_area_row[domain], lane);
    const double thickness = get(state, params.stride,
                                 params.thickness_row[domain], lane);
    const double sign = domain == 0 ? 1.0 : -1.0;
    const double argument = spm_scalar::activationArgument(
      sign, current_density[lane], area, thickness, exchange);
    overpotential[domain] = spm_scalar::activationOverpotential(
      temperature, params.Rg, params.n, params.F, argument);
    electrode_ocv[domain] = curveEval(curves,
                                      params.electrode_ocv[domain],
                                      surface_stoichiometry[domain]);
  }

  constexpr int neg = 0;
  constexpr int pos = 1;
  const double entropic = curveEval(curves,
                                    params.total_entropic,
                                    surface_stoichiometry[pos]);
  const double ocv = spm_scalar::cellOpenCircuitVoltage(
    electrode_ocv[neg],
    electrode_ocv[pos],
    temperature,
    params.reference_temperature,
    entropic);
  const double area_neg = spm_scalar::activeArea(
    get(state, params.stride, params.specific_area_row[neg], lane),
    params.electrode_area,
    get(state, params.stride, params.thickness_row[neg], lane));
  const double area_pos = spm_scalar::activeArea(
    get(state, params.stride, params.specific_area_row[pos], lane),
    params.electrode_area,
    get(state, params.stride, params.thickness_row[pos], lane));
  const double resistance = spm_scalar::seriesResistance(
    get(state, params.stride, params.sei_thickness_row, lane),
    params.sei_resistivity_area,
    get(state, params.stride, params.specific_resistance_row[neg], lane),
    get(state, params.stride, params.specific_resistance_row[pos], lane),
    get(state, params.stride, params.collector_resistance_row, lane),
    area_neg,
    area_pos,
    params.electrode_area);
  const double current = current_density[lane] * params.electrode_area;
  const double voltage = spm_scalar::terminalVoltage(
    ocv, overpotential[neg], overpotential[pos], resistance, current);
  if (!isfinite(voltage)) {
    atomicExch(status, 1);
    terminal_voltage[lane] = nan("");
    return;
  }
  terminal_voltage[lane] = voltage;

  const double dAh = current * dt / 3600.0;
  at(state, params.stride, params.elapsed_time_row, lane) += dt;
  at(state, params.stride, params.charge_throughput_row, lane) += fabs(dAh);
  at(state, params.stride, params.energy_throughput_row, lane) +=
    fabs(dAh * voltage);
}

RuntimeResult fromCuda(cudaError_t result)
{
  return result == cudaSuccess ? RuntimeResult::success
                               : RuntimeResult::cuda_failure;
}

} // namespace

struct Runtime
{
  KernelParams params{};
  cudaStream_t stream{};
  void *allocation{};
  std::size_t allocation_bytes{};
  double *state{};
  double *saved_state{};
  double *current{};
  double *voltage{};
  double *curves{};
  int *status{};
  std::size_t state_values{};
  std::size_t curve_values{};
  std::vector<double> current_cache{};
  bool current_valid{};
  std::size_t allocation_count{};
  std::size_t device_sync_count{};
};

struct RecordingRuntime
{
  Runtime *source{};
  cudaStream_t stream{};
  void *pinned{};
  std::size_t slots{};
  std::size_t values_per_slot{};
  std::vector<cudaEvent_t> compute_ready{};
  std::vector<cudaEvent_t> copy_ready{};
};

int deviceCount() noexcept
{
  int count{};
  return cudaGetDeviceCount(&count) == cudaSuccess ? count : 0;
}

RuntimeResult create(const KernelParams &params,
                     std::span<const double> curves,
                     std::span<const double> state,
                     Runtime *&output) noexcept
{
  if (output != nullptr || params.nch <= 0 || params.nch > maximum_nch
      || params.lanes <= 0 || params.stride < params.lanes || params.rows <= 0
      || state.size() != static_cast<std::size_t>(params.rows)
                           * static_cast<std::size_t>(params.stride)
      || curves.empty())
    return RuntimeResult::invalid_parameter;
  Runtime *runtime{};
  try {
    runtime = new Runtime;
    runtime->params = params;
    runtime->state_values = state.size();
    runtime->curve_values = curves.size();
    runtime->current_cache.resize(static_cast<std::size_t>(params.lanes));
  } catch (...) {
    delete runtime;
    return RuntimeResult::cuda_failure;
  }
  if (cudaStreamCreateWithFlags(&runtime->stream, cudaStreamNonBlocking)
      != cudaSuccess) {
    delete runtime;
    return RuntimeResult::cuda_failure;
  }

  std::size_t cursor{};
  const auto reserve = [&](std::size_t bytes) {
    const auto offset = cursor;
    cursor = alignUp(cursor + bytes);
    return offset;
  };
  const auto state_offset = reserve(state.size_bytes());
  const auto saved_offset = reserve(state.size_bytes());
  const auto current_offset = reserve(static_cast<std::size_t>(params.lanes)
                                      * sizeof(double));
  const auto voltage_offset = reserve(static_cast<std::size_t>(params.lanes)
                                      * sizeof(double));
  const auto curves_offset = reserve(curves.size_bytes());
  const auto status_offset = reserve(sizeof(int));
  runtime->allocation_bytes = cursor;
  if (cudaMalloc(&runtime->allocation, runtime->allocation_bytes)
      != cudaSuccess) {
    cudaStreamDestroy(runtime->stream);
    delete runtime;
    return RuntimeResult::cuda_failure;
  }
  runtime->allocation_count = 1;
  auto *base = static_cast<unsigned char *>(runtime->allocation);
  runtime->state = reinterpret_cast<double *>(base + state_offset);
  runtime->saved_state = reinterpret_cast<double *>(base + saved_offset);
  runtime->current = reinterpret_cast<double *>(base + current_offset);
  runtime->voltage = reinterpret_cast<double *>(base + voltage_offset);
  runtime->curves = reinterpret_cast<double *>(base + curves_offset);
  runtime->status = reinterpret_cast<int *>(base + status_offset);
  if (cudaMemcpyAsync(runtime->state,
                      state.data(),
                      state.size_bytes(),
                      cudaMemcpyHostToDevice,
                      runtime->stream)
        != cudaSuccess
      || cudaMemcpyAsync(runtime->curves,
                         curves.data(),
                         curves.size_bytes(),
                         cudaMemcpyHostToDevice,
                         runtime->stream)
           != cudaSuccess
      || cudaMemsetAsync(runtime->current,
                         0,
                         static_cast<std::size_t>(params.lanes) * sizeof(double),
                         runtime->stream)
           != cudaSuccess
      || cudaMemsetAsync(runtime->status, 0, sizeof(int), runtime->stream)
           != cudaSuccess
      || cudaStreamSynchronize(runtime->stream) != cudaSuccess) {
    destroy(runtime);
    return RuntimeResult::cuda_failure;
  }
  output = runtime;
  return RuntimeResult::success;
}

void destroy(Runtime *runtime) noexcept
{
  if (runtime == nullptr)
    return;
  if (runtime->stream != nullptr)
    cudaStreamDestroy(runtime->stream);
  if (runtime->allocation != nullptr)
    cudaFree(runtime->allocation);
  delete runtime;
}

RuntimeResult uploadState(Runtime *runtime,
                          std::span<const double> state) noexcept
{
  if (runtime == nullptr || state.size() != runtime->state_values)
    return RuntimeResult::invalid_parameter;
  return fromCuda(cudaMemcpyAsync(runtime->state,
                                  state.data(),
                                  state.size_bytes(),
                                  cudaMemcpyHostToDevice,
                                  runtime->stream));
}

RuntimeResult downloadState(Runtime *runtime,
                            std::span<double> state) noexcept
{
  if (runtime == nullptr || state.size() != runtime->state_values)
    return RuntimeResult::invalid_parameter;
  if (cudaMemcpyAsync(state.data(),
                      runtime->state,
                      state.size_bytes(),
                      cudaMemcpyDeviceToHost,
                      runtime->stream)
        != cudaSuccess
      || cudaStreamSynchronize(runtime->stream) != cudaSuccess)
    return RuntimeResult::cuda_failure;
  return RuntimeResult::success;
}

RuntimeResult launchStep(Runtime *runtime,
                         std::span<const double> current_density,
                         double time,
                         double dt) noexcept
{
  if (runtime == nullptr
      || current_density.size() != static_cast<std::size_t>(runtime->params.lanes)
      || !std::isfinite(time) || !std::isfinite(dt) || dt < 0.0)
    return RuntimeResult::invalid_parameter;
  if (!runtime->current_valid
      || !std::equal(current_density.begin(),
                     current_density.end(),
                     runtime->current_cache.begin())) {
    std::copy(current_density.begin(),
              current_density.end(),
              runtime->current_cache.begin());
    runtime->current_valid = true;
    if (cudaMemcpyAsync(runtime->current,
                        current_density.data(),
                        current_density.size_bytes(),
                        cudaMemcpyHostToDevice,
                        runtime->stream)
        != cudaSuccess)
      return RuntimeResult::cuda_failure;
  }
  if (cudaMemsetAsync(runtime->status, 0, sizeof(int), runtime->stream)
      != cudaSuccess)
    return RuntimeResult::cuda_failure;
  constexpr int threads = 256;
  const int blocks = (runtime->params.lanes + threads - 1) / threads;
  exactModalKernel<<<blocks, threads, 0, runtime->stream>>>(runtime->params,
                                                            runtime->state,
                                                            runtime->current,
                                                            runtime->voltage,
                                                            runtime->curves,
                                                            runtime->status,
                                                            dt);
  return fromCuda(cudaPeekAtLastError());
}

RuntimeResult synchronize(Runtime *runtime,
                          std::span<double> terminal_voltage) noexcept
{
  if (runtime == nullptr
      || terminal_voltage.size()
           != static_cast<std::size_t>(runtime->params.lanes))
    return RuntimeResult::invalid_parameter;
  int status{};
  if (cudaMemcpyAsync(terminal_voltage.data(),
                      runtime->voltage,
                      terminal_voltage.size_bytes(),
                      cudaMemcpyDeviceToHost,
                      runtime->stream)
        != cudaSuccess
      || cudaMemcpyAsync(&status,
                         runtime->status,
                         sizeof(status),
                         cudaMemcpyDeviceToHost,
                         runtime->stream)
           != cudaSuccess
      || cudaStreamSynchronize(runtime->stream) != cudaSuccess)
    return RuntimeResult::cuda_failure;
  return status == 0 ? RuntimeResult::success
                     : RuntimeResult::invalid_state;
}

RuntimeResult checkpoint(Runtime *runtime) noexcept
{
  if (runtime == nullptr)
    return RuntimeResult::invalid_parameter;
  return fromCuda(cudaMemcpyAsync(runtime->saved_state,
                                  runtime->state,
                                  runtime->state_values * sizeof(double),
                                  cudaMemcpyDeviceToDevice,
                                  runtime->stream));
}

RuntimeResult restore(Runtime *runtime) noexcept
{
  if (runtime == nullptr)
    return RuntimeResult::invalid_parameter;
  return fromCuda(cudaMemcpyAsync(runtime->state,
                                  runtime->saved_state,
                                  runtime->state_values * sizeof(double),
                                  cudaMemcpyDeviceToDevice,
                                  runtime->stream));
}

std::size_t deviceArenaBytes(const Runtime *runtime) noexcept
{
  return runtime == nullptr ? 0 : runtime->allocation_bytes;
}

std::size_t deviceAllocationCount(const Runtime *runtime) noexcept
{
  return runtime == nullptr ? 0 : runtime->allocation_count;
}

std::size_t deviceWideSynchronizationCount(const Runtime *runtime) noexcept
{
  return runtime == nullptr ? 0 : runtime->device_sync_count;
}

RuntimeResult createRecording(Runtime *runtime,
                              std::size_t slots,
                              RecordingRuntime *&output) noexcept
{
  if (runtime == nullptr || output != nullptr || slots < 3 || slots > 1024)
    return RuntimeResult::invalid_parameter;
  const std::size_t values_per_slot = runtime->state_values
                                      + static_cast<std::size_t>(
                                        runtime->params.lanes);
  if (values_per_slot > std::numeric_limits<std::size_t>::max() / slots
      || values_per_slot * slots
           > std::numeric_limits<std::size_t>::max() / sizeof(double))
    return RuntimeResult::invalid_parameter;
  RecordingRuntime *recording{};
  try {
    recording = new RecordingRuntime;
    recording->source = runtime;
    recording->slots = slots;
    recording->values_per_slot = values_per_slot;
    recording->compute_ready.resize(slots);
    recording->copy_ready.resize(slots);
  } catch (...) {
    delete recording;
    return RuntimeResult::cuda_failure;
  }
  if (cudaStreamCreateWithFlags(&recording->stream, cudaStreamNonBlocking)
      != cudaSuccess) {
    delete recording;
    return RuntimeResult::cuda_failure;
  }
  const std::size_t bytes = values_per_slot * slots * sizeof(double);
  if (cudaHostAlloc(&recording->pinned, bytes, cudaHostAllocPortable)
      != cudaSuccess) {
    cudaStreamDestroy(recording->stream);
    delete recording;
    return RuntimeResult::cuda_failure;
  }
  for (std::size_t slot = 0; slot < slots; ++slot) {
    if (cudaEventCreateWithFlags(&recording->compute_ready[slot],
                                 cudaEventDisableTiming)
          != cudaSuccess
        || cudaEventCreateWithFlags(&recording->copy_ready[slot],
                                    cudaEventDisableTiming)
             != cudaSuccess) {
      destroyRecording(recording);
      return RuntimeResult::cuda_failure;
    }
  }
  output = recording;
  return RuntimeResult::success;
}

void destroyRecording(RecordingRuntime *recording) noexcept
{
  if (recording == nullptr)
    return;
  for (auto event : recording->compute_ready)
    if (event != nullptr)
      cudaEventDestroy(event);
  for (auto event : recording->copy_ready)
    if (event != nullptr)
      cudaEventDestroy(event);
  if (recording->stream != nullptr)
    cudaStreamDestroy(recording->stream);
  if (recording->pinned != nullptr)
    cudaFreeHost(recording->pinned);
  delete recording;
}

RuntimeResult recordSnapshot(RecordingRuntime *recording,
                             std::size_t slot) noexcept
{
  if (recording == nullptr || recording->source == nullptr
      || slot >= recording->slots)
    return RuntimeResult::invalid_parameter;
  auto &source = *recording->source;
  auto *destination = static_cast<double *>(recording->pinned)
                      + slot * recording->values_per_slot;
  auto *state_destination = destination;
  auto *current_destination = destination + source.state_values;
  if (cudaEventRecord(recording->compute_ready[slot], source.stream)
        != cudaSuccess
      || cudaStreamWaitEvent(recording->stream,
                             recording->compute_ready[slot],
                             0)
           != cudaSuccess
      || cudaMemcpyAsync(state_destination,
                         source.state,
                         source.state_values * sizeof(double),
                         cudaMemcpyDeviceToHost,
                         recording->stream)
           != cudaSuccess
      || cudaMemcpyAsync(current_destination,
                         source.current,
                         static_cast<std::size_t>(source.params.lanes)
                           * sizeof(double),
                         cudaMemcpyDeviceToHost,
                         recording->stream)
           != cudaSuccess
      || cudaEventRecord(recording->copy_ready[slot], recording->stream)
           != cudaSuccess
      || cudaStreamWaitEvent(source.stream, recording->copy_ready[slot], 0)
           != cudaSuccess)
    return RuntimeResult::cuda_failure;
  return RuntimeResult::success;
}

RuntimeResult waitSnapshot(RecordingRuntime *recording,
                           std::size_t slot) noexcept
{
  if (recording == nullptr || slot >= recording->slots)
    return RuntimeResult::invalid_parameter;
  return fromCuda(cudaEventSynchronize(recording->copy_ready[slot]));
}

std::span<double> recordingState(RecordingRuntime *recording,
                                 std::size_t slot) noexcept
{
  if (recording == nullptr || recording->source == nullptr
      || slot >= recording->slots)
    return {};
  auto *data = static_cast<double *>(recording->pinned)
               + slot * recording->values_per_slot;
  return { data, recording->source->state_values };
}

std::span<double> recordingCurrent(RecordingRuntime *recording,
                                   std::size_t slot) noexcept
{
  if (recording == nullptr || recording->source == nullptr
      || slot >= recording->slots)
    return {};
  auto *data = static_cast<double *>(recording->pinned)
               + slot * recording->values_per_slot
               + recording->source->state_values;
  return { data,
           static_cast<std::size_t>(recording->source->params.lanes) };
}

bool recordingUsesNonDefaultStream(const RecordingRuntime *recording) noexcept
{
  return recording != nullptr && recording->stream != nullptr;
}

bool recordingUsesPinnedMemory(const RecordingRuntime *recording) noexcept
{
  return recording != nullptr && recording->pinned != nullptr;
}

} // namespace slide::core::cuda_detail
