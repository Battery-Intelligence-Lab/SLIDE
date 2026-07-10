/**
 * @file CudaSpmBatch.cpp
 * @brief CUDA-free host facade and cold parameter flattening.
 */

#include "CudaSpmBatch.hpp"

#include "CudaSpmData.hpp"
#include "Numeric.hpp"
#include "SpectralModel.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <new>
#include <thread>
#include <utility>
#include <vector>

namespace slide::core {
namespace {

slide::Status mapResult(cuda_detail::RuntimeResult result)
{
  switch (result) {
  case cuda_detail::RuntimeResult::success:
    return slide::Status::Success;
  case cuda_detail::RuntimeResult::invalid_parameter:
    return slide::Status::Invalid_parameters;
  case cuda_detail::RuntimeResult::invalid_state:
    return slide::Status::Invalid_states;
  case cuda_detail::RuntimeResult::cuda_failure:
  default:
    return slide::Status::Numerical_failure;
  }
}

#if defined(SLIDE_WITH_CUDA)

cuda_detail::CurveRef appendCurve(const OCVCurve &curve,
                                  std::vector<double> &storage)
{
  const auto offset = static_cast<int>(storage.size());
  storage.insert(storage.end(), curve.stoichiometry.begin(), curve.stoichiometry.end());
  const auto y_offset = static_cast<int>(storage.size());
  storage.insert(storage.end(), curve.value.begin(), curve.value.end());
  return { .x_offset = offset,
           .y_offset = y_offset,
           .count = static_cast<int>(curve.stoichiometry.size()) };
}

cuda_detail::CurveRef appendCurveOrZero(const OCVCurve &curve,
                                        std::vector<double> &storage)
{
  if (!curve.stoichiometry.empty())
    return appendCurve(curve, storage);
  const OCVCurve zero{ .stoichiometry = { 0.0, 1.0 },
                       .value = { 0.0, 0.0 } };
  return appendCurve(zero, storage);
}

template <int NCH>
slide::Status flatten(const SpmFactoryInput &input,
                      const SpmBatch &batch,
                      cuda_detail::KernelParams &params,
                      std::vector<double> &curves)
{
  PerDomain<double> radius{};
  for (const Domain domain : domains)
    radius[domain_index(domain)] = input.design.electrode[domain_index(domain)].particle_radius;
  CompiledSpectralModel<NCH> spectral;
  const auto status = compileSpectralModel<NCH>(radius, spectral);
  if (status != slide::Status::Success)
    return status;

  params.nch = NCH;
  params.lanes = batch.n_lanes();
  params.stride = batch.state().stride();
  params.rows = batch.state().n_rows();
  const auto &layout = batch.layout();
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    params.z_row[d] = layout.spm.z[d].row_begin;
    params.diffusion_row[d] = layout.spm.diffusion_coefficient[d].row_begin;
    params.thickness_row[d] = layout.spm.electrode_thickness[d].row_begin;
    params.specific_area_row[d] = layout.spm.specific_surface_area[d].row_begin;
    params.specific_resistance_row[d] = layout.spm.specific_resistance[d].row_begin;
    params.diffusion_activation[d] = input.design.electrode[d].active_material.D_s.activation_energy;
    params.electrode[d] = {
      .cs_max = input.design.electrode[d].active_material.cs_max,
      .reaction_rate_ref = input.design.electrode[d].active_material.k_ct.reference_value,
      .reaction_activation = input.design.electrode[d].active_material.k_ct.activation_energy
    };
    for (int mode = 0; mode < NCH; ++mode) {
      params.A[d][mode] = spectral.A[d][static_cast<std::size_t>(mode)];
      params.B[d][mode] = spectral.B[d][static_cast<std::size_t>(mode)];
      params.surface_C[d][mode] = spectral.C[d][0][static_cast<std::size_t>(mode)];
    }
    params.surface_D[d] = spectral.D[d][0];
    params.electrode_ocv[d] = appendCurve(
      input.design.electrode[d].active_material.ocv, curves);
  }
  params.total_entropic = appendCurveOrZero(input.total_entropic_coefficient,
                                             curves);
  params.temperature_row = layout.spm.temperature.row_begin;
  params.sei_thickness_row = layout.spm.sei_thickness.row_begin;
  params.collector_resistance_row = layout.spm.current_collector_resistance.row_begin;
  params.elapsed_time_row = layout.elapsed_time.row_begin;
  params.charge_throughput_row = layout.charge_throughput.row_begin;
  params.energy_throughput_row = layout.energy_throughput.row_begin;
  params.electrolyte_concentration = input.design.electrolyte.concentration;
  params.reference_temperature = input.design.thermal.reference_temperature;
  params.electrode_area = input.design.electrode_area;
  params.sei_resistivity_area = input.sei_resistivity_area;
  return slide::Status::Success;
}

#endif

} // namespace

struct CudaSpmBatch::Impl
{
  SpmBatch host{};
  std::vector<real_t> voltage{};
#if defined(SLIDE_WITH_CUDA)
  cuda_detail::Runtime *runtime{};
  ~Impl() { cuda_detail::destroy(runtime); }
#endif
};

CudaSpmBatch::CudaSpmBatch() = default;
CudaSpmBatch::~CudaSpmBatch() = default;
CudaSpmBatch::CudaSpmBatch(CudaSpmBatch &&) noexcept = default;
CudaSpmBatch &CudaSpmBatch::operator=(CudaSpmBatch &&) noexcept = default;

bool CudaSpmBatch::available() noexcept
{
#if defined(SLIDE_WITH_CUDA)
  return cuda_detail::deviceCount() > 0;
#else
  return false;
#endif
}

slide::Status CudaSpmBatch::build(const SpmFactoryInput &input,
                                  const SpmModelOptions &options,
                                  int n_lanes)
{
#if !defined(SLIDE_WITH_CUDA)
  (void)input;
  (void)options;
  (void)n_lanes;
  return slide::Status::NotImplementedYet;
#else
  if (!available() || options.thermal || options.has_ageing() || n_lanes <= 0)
    return slide::Status::Invalid_parameters;
  try {
    auto candidate = std::make_unique<Impl>();
    auto status = buildSpmBatch(input, options, n_lanes, candidate->host);
    if (status != slide::Status::Success)
      return status;
    cuda_detail::KernelParams params;
    std::vector<double> curves;
    curves.reserve(input.design.electrode[0].active_material.ocv.stoichiometry.size() * 2
                   + input.design.electrode[1].active_material.ocv.stoichiometry.size() * 2
                   + std::max<std::size_t>(2, input.total_entropic_coefficient.stoichiometry.size()) * 2);
    switch (options.nch) {
    case 5:
      status = flatten<5>(input, candidate->host, params, curves);
      break;
    case 8:
      status = flatten<8>(input, candidate->host, params, curves);
      break;
    case 12:
      status = flatten<12>(input, candidate->host, params, curves);
      break;
    default:
      return slide::Status::Invalid_parameters;
    }
    if (status != slide::Status::Success)
      return status;
    const auto result = cuda_detail::create(params,
                                             curves,
                                             candidate->host.state().raw(),
                                             candidate->runtime);
    if (result != cuda_detail::RuntimeResult::success)
      return mapResult(result);
    candidate->voltage.resize(static_cast<std::size_t>(n_lanes));
    std::vector<double> zero_current(static_cast<std::size_t>(n_lanes));
    auto runtime_status = cuda_detail::launchStep(candidate->runtime,
                                                   zero_current,
                                                   0.0,
                                                   0.0);
    if (runtime_status == cuda_detail::RuntimeResult::success)
      runtime_status = cuda_detail::synchronize(candidate->runtime,
                                                 candidate->voltage);
    if (runtime_status != cuda_detail::RuntimeResult::success)
      return mapResult(runtime_status);
    impl_ = std::move(candidate);
    return slide::Status::Success;
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  }
#endif
}

bool CudaSpmBatch::valid() const noexcept
{
#if defined(SLIDE_WITH_CUDA)
  return impl_ != nullptr && impl_->host.valid() && impl_->runtime != nullptr;
#else
  return false;
#endif
}

SpmBatch &CudaSpmBatch::hostBatch()
{
  assert(impl_ != nullptr);
  return impl_->host;
}

const SpmBatch &CudaSpmBatch::hostBatch() const
{
  assert(impl_ != nullptr);
  return impl_->host;
}

slide::Status CudaSpmBatch::uploadState()
{
#if defined(SLIDE_WITH_CUDA)
  if (!valid())
    return slide::Status::Invalid_parameters;
  return mapResult(cuda_detail::uploadState(impl_->runtime,
                                             impl_->host.state().raw()));
#else
  return slide::Status::NotImplementedYet;
#endif
}

slide::Status CudaSpmBatch::downloadState()
{
#if defined(SLIDE_WITH_CUDA)
  if (!valid())
    return slide::Status::Invalid_parameters;
  return mapResult(cuda_detail::downloadState(impl_->runtime,
                                               impl_->host.state().raw()));
#else
  return slide::Status::NotImplementedYet;
#endif
}

slide::Status CudaSpmBatch::step(std::span<const real_t> current_density,
                                 real_t time,
                                 real_t dt)
{
#if defined(SLIDE_WITH_CUDA)
  if (!valid() || static_cast<int>(current_density.size()) != nLanes()
      || !is_finite(time) || !is_finite(dt) || !(dt > 0.0))
    return slide::Status::Invalid_parameters;
  for (const auto current : current_density)
    if (!is_finite(current))
      return slide::Status::Invalid_parameters;
  return mapResult(cuda_detail::launchStep(impl_->runtime,
                                            current_density,
                                            time,
                                            dt));
#else
  (void)current_density;
  (void)time;
  (void)dt;
  return slide::Status::NotImplementedYet;
#endif
}

slide::Status CudaSpmBatch::synchronize()
{
#if defined(SLIDE_WITH_CUDA)
  if (!valid())
    return slide::Status::Invalid_parameters;
  return mapResult(cuda_detail::synchronize(impl_->runtime, impl_->voltage));
#else
  return slide::Status::NotImplementedYet;
#endif
}

std::span<const real_t> CudaSpmBatch::terminalVoltage() const
{
  return impl_ == nullptr ? std::span<const real_t>{}
                          : std::span<const real_t>{ impl_->voltage };
}

slide::Status CudaSpmBatch::checkpoint()
{
#if defined(SLIDE_WITH_CUDA)
  return valid() ? mapResult(cuda_detail::checkpoint(impl_->runtime))
                 : slide::Status::Invalid_parameters;
#else
  return slide::Status::NotImplementedYet;
#endif
}

slide::Status CudaSpmBatch::restore()
{
#if defined(SLIDE_WITH_CUDA)
  return valid() ? mapResult(cuda_detail::restore(impl_->runtime))
                 : slide::Status::Invalid_parameters;
#else
  return slide::Status::NotImplementedYet;
#endif
}

int CudaSpmBatch::nLanes() const
{
  return impl_ == nullptr ? 0 : impl_->host.n_lanes();
}

std::size_t CudaSpmBatch::deviceArenaBytes() const noexcept
{
#if defined(SLIDE_WITH_CUDA)
  return valid() ? cuda_detail::deviceBytes(impl_->runtime) : 0;
#else
  return 0;
#endif
}

std::size_t CudaSpmBatch::deviceAllocationCount() const noexcept
{
#if defined(SLIDE_WITH_CUDA)
  return valid() ? cuda_detail::deviceAllocations(impl_->runtime) : 0;
#else
  return 0;
#endif
}

std::size_t CudaSpmBatch::deviceWideSynchronizationCount() const noexcept
{
#if defined(SLIDE_WITH_CUDA)
  return valid() ? cuda_detail::deviceWideSynchronizations(impl_->runtime) : 0;
#else
  return 0;
#endif
}

struct CudaAsyncRecorder::Impl
{
  enum class SlotState : unsigned char { empty,
                                         copying,
                                         draining };

  CudaSpmBatch *batch{};
  AsyncRecorder writer{};
  AsyncRecorderConfig config{};
#if defined(SLIDE_WITH_CUDA)
  cuda_detail::RecordingRuntime *runtime{};
#endif
  std::vector<SlotState> states{};
  std::vector<std::uint64_t> steps{};
  std::size_t write_slot{};
  std::size_t read_slot{};
  std::uint64_t previous_step{};
  bool has_previous_step{};
  bool closing{};
  bool finished{};
  std::mutex mutex{};
  std::condition_variable ready{};
  std::condition_variable space{};
  std::thread worker{};
  std::atomic<std::uint64_t> thinned{};
  std::atomic<int> status{ static_cast<int>(slide::Status::Success) };

  ~Impl()
  {
#if defined(SLIDE_WITH_CUDA)
    cuda_detail::destroyRecording(runtime);
#endif
  }

  slide::Status workerStatus() const
  {
    return static_cast<slide::Status>(status.load(std::memory_order_relaxed));
  }

  void fail(slide::Status failure)
  {
    int expected = static_cast<int>(slide::Status::Success);
    status.compare_exchange_strong(expected,
                                   static_cast<int>(failure),
                                   std::memory_order_relaxed);
  }

  void drainLoop()
  {
#if defined(SLIDE_WITH_CUDA)
    while (true) {
      std::unique_lock lock{ mutex };
      ready.wait(lock, [&] {
        return states[read_slot] == SlotState::copying || closing;
      });
      if (states[read_slot] != SlotState::copying) {
        if (closing)
          break;
        continue;
      }
      const std::size_t slot = read_slot;
      const std::uint64_t accepted_step = steps[slot];
      states[slot] = SlotState::draining;
      lock.unlock();

      auto result = cuda_detail::waitSnapshot(runtime, slot);
      if (result == cuda_detail::RuntimeResult::success) {
        const auto state = cuda_detail::recordingState(runtime, slot);
        const auto current = cuda_detail::recordingCurrent(runtime, slot);
        const auto &host = batch->hostBatch();
        const auto index = static_cast<std::size_t>(
                             host.layout().elapsed_time.row_begin)
                           * static_cast<std::size_t>(host.state().stride());
        if (index >= state.size()) {
          fail(slide::Status::Invalid_states);
        } else {
          const auto write_status = writer.enqueueSnapshot(accepted_step,
                                                            state[index],
                                                            current,
                                                            state);
          if (write_status != slide::Status::Success)
            fail(write_status);
        }
      } else {
        fail(mapResult(result));
      }

      lock.lock();
      states[slot] = SlotState::empty;
      read_slot = (read_slot + 1) % states.size();
      lock.unlock();
      space.notify_all();
    }
#endif
  }
};

CudaAsyncRecorder::CudaAsyncRecorder() = default;

CudaAsyncRecorder::~CudaAsyncRecorder()
{
  (void)finish();
}

slide::Status CudaAsyncRecorder::configure(
  CudaSpmBatch &batch,
  const std::filesystem::path &path,
  AsyncRecorderConfig config)
{
#if !defined(SLIDE_WITH_CUDA)
  (void)batch;
  (void)path;
  (void)config;
  return slide::Status::NotImplementedYet;
#else
  if (impl_ != nullptr || !batch.valid() || path.empty()
      || config.cadence == 0 || config.ring_slots < 3
      || config.ring_slots > 1024)
    return slide::Status::Invalid_parameters;
  try {
    auto candidate = std::make_unique<Impl>();
    candidate->batch = &batch;
    candidate->config = config;
    candidate->states.resize(config.ring_slots, Impl::SlotState::empty);
    candidate->steps.resize(config.ring_slots);
    auto result = cuda_detail::createRecording(batch.impl_->runtime,
                                                config.ring_slots,
                                                candidate->runtime);
    if (result != cuda_detail::RuntimeResult::success)
      return mapResult(result);
    auto writer_config = config;
    writer_config.cadence = 1;
    writer_config.backpressure = AsyncBackpressurePolicy::block;
    auto status = candidate->writer.configure(batch.hostBatch(),
                                               path,
                                               writer_config);
    if (status != slide::Status::Success) {
      cuda_detail::destroyRecording(candidate->runtime);
      candidate->runtime = nullptr;
      return status;
    }
    try {
      candidate->worker = std::thread([object = candidate.get()] {
        object->drainLoop();
      });
    } catch (...) {
      (void)candidate->writer.finish();
      cuda_detail::destroyRecording(candidate->runtime);
      candidate->runtime = nullptr;
      return slide::Status::Numerical_failure;
    }
    impl_ = std::move(candidate);
    return slide::Status::Success;
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  }
#endif
}

slide::Status CudaAsyncRecorder::enqueue(std::uint64_t accepted_step)
{
#if !defined(SLIDE_WITH_CUDA)
  (void)accepted_step;
  return slide::Status::NotImplementedYet;
#else
  if (impl_ == nullptr || impl_->finished)
    return slide::Status::Invalid_parameters;
  if (accepted_step % impl_->config.cadence != 0)
    return slide::Status::Success;
  std::unique_lock lock{ impl_->mutex };
  if (impl_->has_previous_step && accepted_step <= impl_->previous_step)
    return slide::Status::Invalid_parameters;
  impl_->previous_step = accepted_step;
  impl_->has_previous_step = true;
  if (impl_->workerStatus() != slide::Status::Success)
    return impl_->workerStatus();
  auto available = [&] {
    return impl_->states[impl_->write_slot] == Impl::SlotState::empty
           || impl_->closing
           || impl_->workerStatus() != slide::Status::Success;
  };
  if (!available()) {
    if (impl_->config.backpressure == AsyncBackpressurePolicy::thin) {
      impl_->thinned.fetch_add(1, std::memory_order_relaxed);
      return slide::Status::Success;
    }
    impl_->space.wait(lock, available);
  }
  if (impl_->closing || impl_->workerStatus() != slide::Status::Success)
    return impl_->workerStatus() == slide::Status::Success
             ? slide::Status::Invalid_states
             : impl_->workerStatus();
  const std::size_t slot = impl_->write_slot;
  impl_->states[slot] = Impl::SlotState::copying;
  impl_->steps[slot] = accepted_step;
  const auto result = cuda_detail::recordSnapshot(impl_->runtime, slot);
  if (result != cuda_detail::RuntimeResult::success) {
    impl_->states[slot] = Impl::SlotState::empty;
    impl_->fail(mapResult(result));
    return impl_->workerStatus();
  }
  impl_->write_slot = (impl_->write_slot + 1) % impl_->states.size();
  lock.unlock();
  impl_->ready.notify_one();
  return slide::Status::Success;
#endif
}

slide::Status CudaAsyncRecorder::finish()
{
  if (impl_ == nullptr)
    return slide::Status::Success;
  if (impl_->finished)
    return impl_->workerStatus();
  {
    const std::lock_guard lock{ impl_->mutex };
    impl_->closing = true;
  }
  impl_->ready.notify_all();
  impl_->space.notify_all();
  if (impl_->worker.joinable())
    impl_->worker.join();
  const auto writer_status = impl_->writer.finish();
  if (writer_status != slide::Status::Success)
    impl_->fail(writer_status);
#if defined(SLIDE_WITH_CUDA)
  cuda_detail::destroyRecording(impl_->runtime);
  impl_->runtime = nullptr;
#endif
  impl_->finished = true;
  return impl_->workerStatus();
}

std::uint64_t CudaAsyncRecorder::thinnedSnapshots() const
{
  return impl_ == nullptr
           ? 0
           : impl_->thinned.load(std::memory_order_relaxed)
               + impl_->writer.thinnedSnapshots();
}

std::uint64_t CudaAsyncRecorder::snapshotsWritten() const
{
  return impl_ == nullptr ? 0 : impl_->writer.snapshotsWritten();
}

bool CudaAsyncRecorder::usesPinnedMemory() const
{
#if defined(SLIDE_WITH_CUDA)
  return impl_ != nullptr
         && cuda_detail::recordingUsesPinnedMemory(impl_->runtime);
#else
  return false;
#endif
}

bool CudaAsyncRecorder::usesNonDefaultStream() const
{
#if defined(SLIDE_WITH_CUDA)
  return impl_ != nullptr
         && cuda_detail::recordingUsesNonDefaultStream(impl_->runtime);
#else
  return false;
#endif
}

} // namespace slide::core
