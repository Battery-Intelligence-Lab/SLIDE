/**
 * @file SpmFactory.cpp
 * @brief Explicit SPM composition registry and cold parameter compilation.
 */

#include "SpmFactory.hpp"

#include "SpectralModel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

namespace slide::core {

SpmBatch::SpmBatch(void *implementation,
                   EvaluateFn evaluate,
                   ObserveVoltageFn observe_voltage,
                   LinearizeTheveninFn linearize_thevenin,
                   SetLanePeriodFn set_lane_period,
                   GetLanePeriodFn get_lane_period,
                   StoreStressFn store_stress,
                   FusedEulerFn fused_euler,
                   ExponentialFn exponential,
                   DestroyFn destroy,
                   int nch,
                   real_t capacity_Ah,
                   real_t electrode_area,
                   SpmComposition composition,
                   SpmPipelineLayout layout,
                   StateArena state,
                   StateArena derivative,
                   std::vector<StateRole>
                     roles)
  : implementation_{ implementation }, evaluate_{ evaluate },
    observe_voltage_{ observe_voltage },
    linearize_thevenin_{ linearize_thevenin },
    set_lane_period_{ set_lane_period },
    get_lane_period_{ get_lane_period },
    store_stress_{ store_stress }, fused_euler_{ fused_euler },
    exponential_{ exponential }, destroy_{ destroy }, nch_{ nch },
    capacity_Ah_{ capacity_Ah }, electrode_area_{ electrode_area },
    composition_{ composition }, layout_{ layout }, state_{ std::move(state) },
    derivative_{ std::move(derivative) }, roles_{ std::move(roles) }
{}

SpmBatch::~SpmBatch()
{
  reset();
}

SpmBatch::SpmBatch(SpmBatch &&other) noexcept
  : implementation_{ std::exchange(other.implementation_, nullptr) },
    evaluate_{ std::exchange(other.evaluate_, nullptr) },
    observe_voltage_{ std::exchange(other.observe_voltage_, nullptr) },
    linearize_thevenin_{ std::exchange(other.linearize_thevenin_, nullptr) },
    set_lane_period_{ std::exchange(other.set_lane_period_, nullptr) },
    get_lane_period_{ std::exchange(other.get_lane_period_, nullptr) },
    store_stress_{ std::exchange(other.store_stress_, nullptr) },
    fused_euler_{ std::exchange(other.fused_euler_, nullptr) },
    exponential_{ std::exchange(other.exponential_, nullptr) },
    destroy_{ std::exchange(other.destroy_, nullptr) },
    nch_{ std::exchange(other.nch_, 0) },
    capacity_Ah_{ std::exchange(other.capacity_Ah_, 0.0) },
    electrode_area_{ std::exchange(other.electrode_area_, 0.0) },
    composition_{ other.composition_ },
    layout_{ other.layout_ }, state_{ std::move(other.state_) },
    derivative_{ std::move(other.derivative_) }, roles_{ std::move(other.roles_) }
{}

SpmBatch &SpmBatch::operator=(SpmBatch &&other) noexcept
{
  if (this != &other) {
    reset();
    implementation_ = std::exchange(other.implementation_, nullptr);
    evaluate_ = std::exchange(other.evaluate_, nullptr);
    observe_voltage_ = std::exchange(other.observe_voltage_, nullptr);
    linearize_thevenin_ = std::exchange(other.linearize_thevenin_, nullptr);
    set_lane_period_ = std::exchange(other.set_lane_period_, nullptr);
    get_lane_period_ = std::exchange(other.get_lane_period_, nullptr);
    store_stress_ = std::exchange(other.store_stress_, nullptr);
    fused_euler_ = std::exchange(other.fused_euler_, nullptr);
    exponential_ = std::exchange(other.exponential_, nullptr);
    destroy_ = std::exchange(other.destroy_, nullptr);
    nch_ = std::exchange(other.nch_, 0);
    capacity_Ah_ = std::exchange(other.capacity_Ah_, 0.0);
    electrode_area_ = std::exchange(other.electrode_area_, 0.0);
    composition_ = other.composition_;
    layout_ = other.layout_;
    state_ = std::move(other.state_);
    derivative_ = std::move(other.derivative_);
    roles_ = std::move(other.roles_);
  }
  return *this;
}

void SpmBatch::reset() noexcept
{
  if (implementation_ != nullptr)
    destroy_(implementation_);
  implementation_ = nullptr;
  evaluate_ = nullptr;
  observe_voltage_ = nullptr;
  linearize_thevenin_ = nullptr;
  set_lane_period_ = nullptr;
  get_lane_period_ = nullptr;
  store_stress_ = nullptr;
  fused_euler_ = nullptr;
  exponential_ = nullptr;
  destroy_ = nullptr;
  nch_ = 0;
  capacity_Ah_ = 0.0;
  electrode_area_ = 0.0;
  state_ = {};
  derivative_ = {};
  roles_.clear();
}

slide::Status SpmBatch::rhs(std::span<const real_t> trial_state,
                            std::span<real_t>
                              trial_derivative,
                            const StepCtx &ctx)
{
  if (!valid() || trial_state.size() != state_.size()
      || trial_derivative.size() != derivative_.size()
      || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt))
    return slide::Status::Invalid_parameters;
  RhsViews views{ BatchShape::from(state_) };
  views.rebind(trial_state, trial_derivative);
  return evaluate_(implementation_, views, ctx);
}

slide::Status SpmBatch::evaluate(const StepCtx &ctx)
{
  return rhs(std::span<const real_t>{ state_.raw() }, derivative_.raw(), ctx);
}

slide::Status SpmBatch::fusedEuler(const StepCtx &ctx, real_t dt,
                                   std::span<real_t> terminal_voltage)
{
  if (!valid() || fused_euler_ == nullptr || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(terminal_voltage.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt) || !is_finite(dt) || !(dt > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView state_view{ BatchShape::from(state_), state_.raw() };
  return fused_euler_(implementation_, state_view, ctx, dt, terminal_voltage);
}

slide::Status SpmBatch::exponentialStep(const StepCtx &ctx, real_t dt,
                                        std::span<real_t> terminal_voltage)
{
  if (!valid() || exponential_ == nullptr
      || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(terminal_voltage.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt)
      || !is_finite(dt) || !(dt > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView state_view{ BatchShape::from(state_), state_.raw() };
  BatchView derivative_view{ BatchShape::from(derivative_), derivative_.raw() };
  return exponential_(implementation_, state_view, derivative_view, ctx, dt, terminal_voltage);
}

slide::Status SpmBatch::terminalVoltage(const StepCtx &ctx,
                                        std::span<real_t>
                                          output)
{
  if (!valid() || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(output.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt))
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return observe_voltage_(implementation_, state_view, ctx, output);
}

slide::Status SpmBatch::terminalVoltageAt(
  std::span<const real_t> snapshot,
  std::span<const real_t>
    current_density,
  std::span<real_t>
    output)
{
  if (!valid() || observe_voltage_ == nullptr || snapshot.size() != state_.size()
      || static_cast<int>(current_density.size()) != n_lanes()
      || output.size() != current_density.size())
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_), snapshot };
  const StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = current_density };
  return observe_voltage_(implementation_, state_view, ctx, output);
}

slide::Status SpmBatch::linearizeThevenin(
  std::span<const real_t> current,
  std::span<real_t>
    intercept_ocv,
  std::span<real_t>
    resistance)
{
  if (!valid() || linearize_thevenin_ == nullptr
      || static_cast<int>(current.size()) != n_lanes()
      || intercept_ocv.size() != current.size()
      || resistance.size() != current.size())
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return linearize_thevenin_(implementation_, state_view, current, intercept_ocv, resistance);
}

slide::Status SpmBatch::setTrustedLanePeriod(int maximum_period)
{
  if (!valid() || set_lane_period_ == nullptr || maximum_period <= 0
      || n_lanes() % maximum_period != 0)
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return set_lane_period_(implementation_, state_view, maximum_period);
}

int SpmBatch::trustedLanePeriod() const
{
  return valid() && get_lane_period_ != nullptr
           ? get_lane_period_(implementation_)
           : n_lanes();
}

slide::Status SpmBatch::storeStressHistory(real_t interval)
{
  if (!valid() || !is_finite(interval) || !(interval > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView view{ BatchShape::from(state_), state_.raw() };
  store_stress_(implementation_, view, interval);
  return slide::Status::Success;
}

struct SpmBatchFactoryAccess
{
  template <class Pipeline>
  static SpmBatch make(Pipeline &&pipeline,
                       int nch,
                       real_t capacity_Ah,
                       real_t electrode_area,
                       SpmComposition composition,
                       SpmPipelineLayout layout,
                       StateArena state,
                       StateArena derivative,
                       std::vector<StateRole>
                         roles)
  {
    using Concrete = std::remove_cvref_t<Pipeline>;
    auto *implementation = new Concrete(std::forward<Pipeline>(pipeline));
    SpmBatch::FusedEulerFn fused_euler = nullptr;
    if constexpr (Concrete::supports_fused_euler)
      fused_euler = [](void *object, BatchView state_view, const StepCtx &ctx, real_t dt, std::span<real_t> terminal_voltage) {
        return static_cast<Concrete *>(object)->advanceEuler(
          state_view, ctx, dt, terminal_voltage);
      };
    return { implementation,
             [](void *object, RhsViews &views, const StepCtx &ctx) {
               return static_cast<Concrete *>(object)->evaluate(views, ctx);
             },
             [](void *object, const ConstBatchView &state_view, const StepCtx &ctx, std::span<real_t> output) {
               return static_cast<Concrete *>(object)->observeTerminalVoltage(
                 state_view, ctx, output);
             },
             [](void *object, const ConstBatchView &state_view, std::span<const real_t> current, std::span<real_t> intercept_ocv, std::span<real_t> resistance) {
               return static_cast<Concrete *>(object)->linearizeThevenin(
                 state_view, current, intercept_ocv, resistance);
             },
             [](void *object, const ConstBatchView &state_view, int maximum_period) {
               return static_cast<Concrete *>(object)->setTrustedLanePeriod(
                 state_view, maximum_period);
             },
             [](const void *object) {
               return static_cast<const Concrete *>(object)->trustedLanePeriod();
             },
             [](void *object, BatchView state_view, real_t interval) {
               static_cast<Concrete *>(object)->storeStressHistory(state_view, interval);
             },
             fused_euler,
             [](void *object, BatchView state_view, BatchView derivative_view, const StepCtx &ctx, real_t dt, std::span<real_t> terminal_voltage) {
               return static_cast<Concrete *>(object)->advanceExponential(
                 state_view, derivative_view, ctx, dt, terminal_voltage);
             },
             [](void *object) { delete static_cast<Concrete *>(object); },
             nch,
             capacity_Ah,
             electrode_area,
             composition,
             layout,
             std::move(state),
             std::move(derivative),
             std::move(roles) };
  }
};

namespace {

  bool validOptions(const SpmModelOptions &options)
  {
    return valid_optional_ageing_model_mask<4>(options.sei_model_mask)
           && valid_optional_ageing_model_mask<5>(options.surface_crack_model_mask)
           && valid_optional_ageing_model_mask<4>(options.lam_model_mask)
           && (!options.sei_porosity || options.sei_model_mask != 0)
           && (!options.surface_crack_diffusivity
               || options.surface_crack_model_mask != 0);
  }

  bool optionsNeedStress(const SpmModelOptions &options)
  {
    constexpr std::uint8_t crack_stress_models = surface_crack_model_bit(1) | surface_crack_model_bit(2);
    return (options.surface_crack_model_mask & crack_stress_models) != 0
           || (options.lam_model_mask & lam_model_bit(1)) != 0;
  }

  slide::Status buildCurveOrZero(const OCVCurve &source,
                                 IndexedPiecewiseLinear &output)
  {
    if (source.stoichiometry.empty() && source.value.empty()) {
      constexpr std::array x{ 0.0, 1.0 };
      constexpr std::array y{ 0.0, 0.0 };
      return output.build(x, y);
    }
    return output.build(source.stoichiometry, source.value);
  }

  slide::Status validateInput(const SpmFactoryInput &input,
                              const SpmModelOptions &options,
                              int n_lanes)
  {
    if (n_lanes <= 0 || !validOptions(options))
      return slide::Status::Invalid_parameters;
    const auto &design = input.design;
    const double reference_temperature = design.thermal.reference_temperature;
    if (!(is_finite(reference_temperature) && reference_temperature > 0.0
          && is_finite(design.capacity_Ah) && design.capacity_Ah > 0.0
          && is_finite(design.electrode_area) && design.electrode_area > 0.0
          && is_finite(design.electrolyte.concentration)
          && design.electrolyte.concentration > 0.0
          && is_finite(input.initial_soc) && input.initial_soc >= 0.0
          && input.initial_soc <= 1.0 && is_finite(input.initial_temperature)
          && input.initial_temperature >= 0.0
          && is_finite(input.initial_sei_thickness)
          && input.initial_sei_thickness > 0.0
          && is_finite(input.initial_lost_lithium)
          && input.initial_lost_lithium >= 0.0
          && is_finite(input.initial_crack_surface_fraction)
          && input.initial_crack_surface_fraction >= 0.0
          && is_finite(input.initial_plated_lithium_thickness)
          && input.initial_plated_lithium_thickness >= 0.0
          && is_finite(input.initial_current_collector_resistance)
          && input.initial_current_collector_resistance >= 0.0
          && is_finite(input.initial_stress_interval)
          && input.initial_stress_interval > 0.0
          && is_finite(input.sei_resistivity_area)
          && input.sei_resistivity_area >= 0.0))
      return slide::Status::Invalid_parameters;

    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      const auto &electrode = design.electrode[d];
      const auto &material = electrode.active_material;
      const double initial_stoichiometry = material.x_0 + input.initial_soc * (material.x_100 - material.x_0);
      if (!(is_finite(electrode.thickness) && electrode.thickness > 0.0
            && is_finite(electrode.porosity) && electrode.porosity >= 0.0
            && electrode.porosity < 1.0 && is_finite(electrode.active_fraction)
            && electrode.active_fraction > 0.0 && electrode.active_fraction <= 1.0
            && is_finite(electrode.particle_radius)
            && electrode.particle_radius > 0.0 && is_finite(material.cs_max)
            && material.cs_max > 0.0 && is_finite(material.x_0)
            && is_finite(material.x_100) && is_finite(initial_stoichiometry)
            && initial_stoichiometry > 0.0 && initial_stoichiometry < 1.0
            && is_finite(material.D_s.reference_value)
            && material.D_s.reference_value > 0.0
            && is_finite(material.D_s.activation_energy)
            && material.D_s.reference_temperature == reference_temperature
            && is_finite(material.k_ct.reference_value)
            && material.k_ct.reference_value >= 0.0
            && is_finite(material.k_ct.activation_energy)
            && material.k_ct.reference_temperature == reference_temperature
            && is_finite(input.initial_specific_resistance[d])
            && input.initial_specific_resistance[d] >= 0.0))
        return slide::Status::Invalid_parameters;
      if (optionsNeedStress(options)
          && !(is_finite(electrode.stress.youngs_modulus)
               && electrode.stress.youngs_modulus > 0.0
               && is_finite(electrode.stress.poisson_ratio)
               && electrode.stress.poisson_ratio < 1.0
               && is_finite(electrode.stress.partial_molar_volume)))
        return slide::Status::Invalid_parameters;
    }
    return slide::Status::Success;
  }

  template <int NCH>
  slide::Status compileParameters(const SpmFactoryInput &input,
                                  const SpmModelOptions &options,
                                  CompiledSpectralModel<NCH> &spectral,
                                  SpmPipelineParams<NCH> &output)
  {
    SpmPipelineParams<NCH> params;
    PerDomain<double> radius{};
    for (const Domain domain : domains)
      radius[domain_index(domain)] = input.design.electrode[domain_index(domain)].particle_radius;
    auto status = compileSpectralModel<NCH>(radius, spectral);
    if (status != slide::Status::Success)
      return status;

    params.diffusion.A = spectral.A;
    params.diffusion.B = spectral.B;
    auto &electrical = params.electrical;
    electrical.electrolyte_concentration = input.design.electrolyte.concentration;
    electrical.reference_temperature = input.design.thermal.reference_temperature;
    electrical.electrode_area = input.design.electrode_area;
    electrical.sei_resistivity_area = input.sei_resistivity_area;
    electrical.concentration.F = electrical.F;
    electrical.concentration.Rg = electrical.Rg;
    electrical.concentration.n = electrical.n;
    electrical.concentration.T_ref = electrical.reference_temperature;
    electrical.concentration.C = spectral.C;
    electrical.concentration.Dout = spectral.D;
    electrical.concentration.Cc = spectral.Cc;
    electrical.concentration.cc_coeff = spectral.cc_coeff;

    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      const auto &source = input.design.electrode[d];
      const auto &material = source.active_material;
      electrical.concentration.R[d] = source.particle_radius;
      electrical.concentration.D_T[d] = material.D_s.activation_energy;
      electrical.electrode[d] = { .cs_max = material.cs_max,
                                  .x_0 = material.x_0,
                                  .x_100 = material.x_100,
                                  .thickness = source.thickness,
                                  .active_fraction = source.active_fraction,
                                  .particle_radius = source.particle_radius,
                                  .specific_surface_area = 3.0 * source.active_fraction / source.particle_radius,
                                  .diffusion_ref = material.D_s.reference_value,
                                  .diffusion_activation = material.D_s.activation_energy,
                                  .reaction_rate_ref = material.k_ct.reference_value,
                                  .reaction_activation = material.k_ct.activation_energy };
      status = electrical.electrode_ocv[d].build(material.ocv.stoichiometry,
                                                 material.ocv.value);
      if (status != slide::Status::Success)
        return status;
    }
    status = buildCurveOrZero(input.total_entropic_coefficient,
                              electrical.total_entropic_coefficient);
    if (status != slide::Status::Success)
      return status;
    status = buildCurveOrZero(input.negative_entropic_coefficient,
                              electrical.negative_entropic_coefficient);
    if (status != slide::Status::Success)
      return status;

    if (options.thermal) {
      status = compileThermalLumped(input.design.thermal, params.thermal);
      if (status != slide::Status::Success)
        return status;
    }

    if (options.has_ageing()) {
      params.enable_stress = optionsNeedStress(options);
      params.enable_sei = options.sei_model_mask != 0;
      params.enable_surface_crack = options.surface_crack_model_mask != 0;
      params.enable_lam = options.lam_model_mask != 0;
      params.enable_lithium_plating = options.lithium_plating;
      if (params.enable_stress) {
        params.stress.x_inner = spectral.x_inner;
        params.stress.integration = spectral.integration;
        for (const Domain domain : domains) {
          const auto d = domain_index(domain);
          const auto &stress = input.design.electrode[d].stress;
          params.stress.partial_molar_volume[d] = stress.partial_molar_volume;
          params.stress.youngs_modulus[d] = stress.youngs_modulus;
          params.stress.poisson_ratio[d] = stress.poisson_ratio;
        }
        status = buildCurveOrZero(input.negative_laresgoiti_stress,
                                  params.stress.laresgoiti_negative);
        if (status != slide::Status::Success)
          return status;
        status = validateSpmStressParams(params.stress);
        if (status != slide::Status::Success)
          return status;
      }

      const auto neg = domain_index(Domain::neg);
      params.sei.mechanism = input.sei;
      params.sei.mechanism.model_mask = options.sei_model_mask;
      params.sei.mechanism.reduce_active_fraction = options.sei_porosity;
      params.sei.mechanism.F = electrical.F;
      params.sei.mechanism.Rg = electrical.Rg;
      params.sei.mechanism.n = electrical.n;
      params.sei.mechanism.reference_temperature = electrical.reference_temperature;
      params.sei.mechanism.electrode_area = electrical.electrode_area;
      params.sei.mechanism.negative_particle_radius = input.design.electrode[neg].particle_radius;
      params.sei.mechanism.sei_resistivity_area = input.sei_resistivity_area;
      params.sei.negative_input_map = spectral.B[neg];
      if (options.sei_model_mask != 0) {
        status = validateSeiParams(params.sei.mechanism);
        if (status != slide::Status::Success)
          return status;
      }

      params.surface_crack.mechanism = input.surface_crack;
      params.surface_crack.mechanism.model_mask = options.surface_crack_model_mask;
      params.surface_crack.mechanism.reduce_negative_diffusivity = options.surface_crack_diffusivity;
      params.surface_crack.mechanism.F = electrical.F;
      params.surface_crack.mechanism.Rg = electrical.Rg;
      params.surface_crack.mechanism.reference_temperature = electrical.reference_temperature;
      params.surface_crack.mechanism.electrode_area = electrical.electrode_area;
      params.surface_crack.mechanism.negative_cs_max = input.design.electrode[neg].active_material.cs_max;
      params.surface_crack.mechanism.sei_resistivity_area = input.sei_resistivity_area;
      const double negative_area = electrical.electrode[neg].specific_surface_area
                                   * electrical.electrode_area * electrical.electrode[neg].thickness;
      if (!(params.surface_crack.mechanism.model4_max_surface > 0.0))
        params.surface_crack.mechanism.model4_max_surface = 5.0 * input.initial_crack_surface_fraction * negative_area;
      params.surface_crack.negative_input_map = spectral.B[neg];
      if (options.surface_crack_model_mask != 0) {
        status = validateSurfaceCrackParams(params.surface_crack.mechanism);
        if (status != slide::Status::Success)
          return status;
      }

      params.lam = input.lam;
      params.lam.model_mask = options.lam_model_mask;
      params.lam.F = electrical.F;
      params.lam.Rg = electrical.Rg;
      params.lam.n = electrical.n;
      params.lam.reference_temperature = electrical.reference_temperature;
      params.lam.particle_radius = radius;
      status = params.lam.positive_ocv.build(
        input.design.electrode[domain_index(Domain::pos)].active_material.ocv.stoichiometry,
        input.design.electrode[domain_index(Domain::pos)].active_material.ocv.value);
      if (status != slide::Status::Success)
        return status;
      if (options.lam_model_mask != 0) {
        status = validateLamParams(params.lam);
        if (status != slide::Status::Success)
          return status;
      }

      params.lithium_plating.mechanism = input.lithium_plating;
      params.lithium_plating.mechanism.F = electrical.F;
      params.lithium_plating.mechanism.Rg = electrical.Rg;
      params.lithium_plating.mechanism.n = electrical.n;
      params.lithium_plating.mechanism.reference_temperature = electrical.reference_temperature;
      params.lithium_plating.mechanism.electrode_area = electrical.electrode_area;
      params.lithium_plating.mechanism.sei_resistivity_area = input.sei_resistivity_area;
      if (!options.lithium_plating)
        params.lithium_plating.mechanism.reaction_rate_ref = 0.0;
      params.lithium_plating.negative_input_map = spectral.B[neg];
      if (options.lithium_plating) {
        status = validateLithiumPlatingParams(
          params.lithium_plating.mechanism);
        if (status != slide::Status::Success)
          return status;
      }
    }

    output = std::move(params);
    return slide::Status::Success;
  }

  template <int NCH, bool WithThermal, bool WithAgeing>
  slide::Status buildRegistered(const SpmFactoryInput &input,
                                const SpmModelOptions &options,
                                int n_lanes,
                                SpmComposition composition,
                                SpmBatch &output)
  {
    using Pipeline = SpmPipeline<NCH,
                                 WithThermal,
                                 WithAgeing,
                                 WithAgeing,
                                 WithAgeing,
                                 WithAgeing>;
    CompiledSpectralModel<NCH> spectral;
    SpmPipelineParams<NCH> params;
    auto status = compileParameters(input, options, spectral, params);
    if (status != slide::Status::Success)
      return status;

    BatchBuilder builder;
    const auto layout = Pipeline::declareLayout(builder);
    std::vector<StateRole> roles(builder.roles().begin(), builder.roles().end());
    StateArena state = builder.build(n_lanes);
    StateArena derivative{ state.n_rows(), n_lanes };

    const double temperature = input.initial_temperature > 0.0
                                 ? input.initial_temperature
                                 : input.design.thermal.reference_temperature;
    for (int lane = 0; lane < n_lanes; ++lane) {
      state.at(layout.spm.temperature, 0, lane) = temperature;
      state.at(layout.spm.sei_thickness, 0, lane) = input.initial_sei_thickness;
      state.at(layout.spm.lost_lithium, 0, lane) = input.initial_lost_lithium;
      state.at(layout.spm.plated_lithium_thickness, 0, lane) = input.initial_plated_lithium_thickness;
      state.at(layout.spm.current_collector_resistance, 0, lane) = input.initial_current_collector_resistance;
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const auto &electrode = input.design.electrode[d];
        const auto &material = electrode.active_material;
        const double stoichiometry = material.x_0 + input.initial_soc * (material.x_100 - material.x_0);
        const double concentration = stoichiometry * material.cs_max;
        const int zero = spectral.zero_mode[d];
        double uniform_mode{};
        for (int node = 0; node < NCH; ++node) {
          const double radial_state = electrode.particle_radius * concentration
                                      * spectral.x_inner[static_cast<std::size_t>(node)];
          uniform_mode += spectral.state_transform[d][static_cast<std::size_t>(zero)]
                                                  [static_cast<std::size_t>(node)]
                          * radial_state;
        }
        state.at(layout.spm.z[d], zero, lane) = uniform_mode;
        state.at(layout.spm.active_fraction[d], 0, lane) = electrode.active_fraction;
        state.at(layout.spm.diffusion_coefficient[d], 0, lane) = material.D_s.reference_value;
        state.at(layout.spm.electrode_thickness[d], 0, lane) = electrode.thickness;
        state.at(layout.spm.specific_surface_area[d], 0, lane) = 3.0 * electrode.active_fraction / electrode.particle_radius;
        state.at(layout.spm.specific_resistance[d], 0, lane) = input.initial_specific_resistance[d];
      }
      const auto neg = domain_index(Domain::neg);
      state.at(layout.spm.crack_surface, 0, lane) = input.initial_crack_surface_fraction
                                                    * state.at(layout.spm.specific_surface_area[neg], 0, lane)
                                                    * input.design.electrode_area
                                                    * state.at(layout.spm.electrode_thickness[neg], 0, lane);
      if constexpr (WithAgeing)
        state.at(layout.stress_history.interval, 0, lane) = input.initial_stress_interval;
    }

    Pipeline pipeline{ std::move(params), layout, n_lanes };
    std::vector<real_t> zero_current(static_cast<std::size_t>(n_lanes));
    const StepCtx initial_ctx{ .time = 0.0, .dt = 0.0, .i_app = zero_current };
    RhsViews views{ BatchShape::from(state) };
    views.rebind(state.raw(), derivative.raw());
    status = pipeline.evaluate(views, initial_ctx);
    if (status != slide::Status::Success)
      return status;
    if constexpr (WithAgeing) {
      BatchView state_view{ BatchShape::from(state), state.raw() };
      pipeline.storeStressHistory(state_view, input.initial_stress_interval);
    }

    SpmBatch candidate = SpmBatchFactoryAccess::make(std::move(pipeline),
                                                     NCH,
                                                     input.design.capacity_Ah,
                                                     input.design.electrode_area,
                                                     composition,
                                                     layout,
                                                     std::move(state),
                                                     std::move(derivative),
                                                     std::move(roles));
    output = std::move(candidate);
    return slide::Status::Success;
  }

  template <int NCH>
  slide::Status selectComposition(const SpmFactoryInput &input,
                                  const SpmModelOptions &options,
                                  int n_lanes,
                                  SpmBatch &output)
  {
    if (options.has_ageing()) {
      if (options.thermal)
        return buildRegistered<NCH, true, true>(input,
                                                options,
                                                n_lanes,
                                                SpmComposition::thermal_ageing,
                                                output);
      return buildRegistered<NCH, false, true>(input,
                                               options,
                                               n_lanes,
                                               SpmComposition::isothermal_ageing,
                                               output);
    }
    if (options.thermal)
      return buildRegistered<NCH, true, false>(input,
                                               options,
                                               n_lanes,
                                               SpmComposition::thermal,
                                               output);
    return buildRegistered<NCH, false, false>(input,
                                              options,
                                              n_lanes,
                                              SpmComposition::isothermal,
                                              output);
  }

} // namespace

slide::Status buildSpmBatch(const SpmFactoryInput &input,
                            const SpmModelOptions &options,
                            int n_lanes,
                            SpmBatch &output)
{
  const auto status = validateInput(input, options, n_lanes);
  if (status != slide::Status::Success)
    return status;
  switch (options.nch) {
  case 5:
    return selectComposition<5>(input, options, n_lanes, output);
  case 8:
    return selectComposition<8>(input, options, n_lanes, output);
  case 12:
    return selectComposition<12>(input, options, n_lanes, output);
  default:
    return slide::Status::Invalid_parameters;
  }
}

} // namespace slide::core
