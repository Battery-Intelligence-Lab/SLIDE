/**
 * @file SpmFactory.cpp
 * @brief Explicit SPM composition registry and cold parameter compilation.
 */

#include "SpmFactory.hpp"

#include "SpectralModel.hpp"
#include "SpmPipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

namespace slide::core {


struct SpmBatchFactoryAccess
{
  template <class Pipeline>
  static SpmBatch make(Pipeline &&pipeline,
                       int nch,
                       real_t capacity_Ah,
                       real_t electrode_area,
                       SpmComposition composition,
                       SpmBatchLayout layout,
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

  /**
   * Propagate the batch-owned constants into every mechanism member that
   * declares the same concept. Call this only after copying the mechanism:
   * adding one of these names to a parameter type intentionally opts it in.
   */
  template <int NCH, class Mechanism>
  void applySharedConstants(const SpmElectricalParams<NCH> &electrical,
                            Mechanism &mechanism)
  {
    mechanism.F = electrical.F;
    mechanism.Rg = electrical.Rg;
    mechanism.reference_temperature = electrical.reference_temperature;
    if constexpr (requires { mechanism.n; })
      mechanism.n = electrical.n;
    if constexpr (requires { mechanism.electrode_area; })
      mechanism.electrode_area = electrical.electrode_area;
    if constexpr (requires { mechanism.sei_resistivity_area; })
      mechanism.sei_resistivity_area = electrical.sei_resistivity_area;
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
      const auto &negative_input_map = spectral.B[neg];
      params.sei.mechanism = input.sei;
      params.sei.mechanism.model_mask = options.sei_model_mask;
      params.sei.mechanism.reduce_active_fraction = options.sei_porosity;
      applySharedConstants(electrical, params.sei.mechanism);
      params.sei.mechanism.negative_particle_radius = input.design.electrode[neg].particle_radius;
      params.sei.negative_input_map = negative_input_map;
      if (options.sei_model_mask != 0) {
        status = validateSeiParams(params.sei.mechanism);
        if (status != slide::Status::Success)
          return status;
      }

      params.surface_crack.mechanism = input.surface_crack;
      params.surface_crack.mechanism.model_mask = options.surface_crack_model_mask;
      params.surface_crack.mechanism.reduce_negative_diffusivity = options.surface_crack_diffusivity;
      applySharedConstants(electrical, params.surface_crack.mechanism);
      params.surface_crack.mechanism.negative_cs_max = input.design.electrode[neg].active_material.cs_max;
      const double negative_area = electrical.electrode[neg].specific_surface_area
                                   * electrical.electrode_area * electrical.electrode[neg].thickness;
      if (!(params.surface_crack.mechanism.model4_max_surface > 0.0))
        params.surface_crack.mechanism.model4_max_surface = 5.0 * input.initial_crack_surface_fraction * negative_area;
      params.surface_crack.negative_input_map = negative_input_map;
      if (options.surface_crack_model_mask != 0) {
        status = validateSurfaceCrackParams(params.surface_crack.mechanism);
        if (status != slide::Status::Success)
          return status;
      }

      params.lam = input.lam;
      params.lam.model_mask = options.lam_model_mask;
      applySharedConstants(electrical, params.lam);
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
      applySharedConstants(electrical, params.lithium_plating.mechanism);
      if (!options.lithium_plating)
        params.lithium_plating.mechanism.reaction_rate_ref = 0.0;
      params.lithium_plating.negative_input_map = negative_input_map;
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
