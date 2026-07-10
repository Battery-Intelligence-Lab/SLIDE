/**
 * @file CudaSpmBatch.cpp
 * @brief CUDA-free host facade and cold parameter flattening.
 */

#include "CudaSpmBatch.hpp"

#include "CudaSpmData.hpp"
#include "Numeric.hpp"
#include "SpectralModel.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <new>
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

} // namespace slide::core
