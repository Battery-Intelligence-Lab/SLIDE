/**
 * @file SpmPipeline.hpp
 * @brief Compile-time composed, allocation-free SPM RHS pipeline.
 */

#pragma once

#include "AgeingKernel.hpp"
#include "Lam.hpp"
#include "LithiumPlating.hpp"
#include "Sei.hpp"
#include "SpmScalarKernels.hpp"
#include "SurfaceCrack.hpp"
#include "ThermalLumped.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace slide::core {

template <int NCH>
struct SpmDiffusionRhsParams
{
  PerDomain<std::array<real_t, NCH>> A{};
  PerDomain<std::array<real_t, NCH>> B{};
};

/** Add the two solid-diffusion modal systems using the shared observable reconstruction. */
template <int NCH, class Real>
void addSpmDiffusionRhs(const SpmDiffusionRhsParams<NCH> &p,
                        const BasicBatchView<const Real> &state,
                        BasicBatchView<Real>
                          derivative,
                        const SpmStateLayout &layout,
                        const BasicSpmObservables<Real> &observables)
{
  const int lanes = state.n_lanes();
  assert(derivative.n_lanes() == lanes && derivative.n_rows() == state.n_rows());
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    assert(layout.z[d].rows == NCH
           && static_cast<int>(observables.effective_diffusivity[d].size()) == lanes
           && static_cast<int>(observables.molar_flux[d].size()) == lanes);
    for (int mode = 0; mode < NCH; ++mode) {
      const Real A = p.A[d][static_cast<std::size_t>(mode)];
      const Real B = p.B[d][static_cast<std::size_t>(mode)];
      for (int lane = 0; lane < lanes; ++lane) {
        const auto i = static_cast<std::size_t>(lane);
        derivative.at(layout.z[d], mode, lane) += spm_scalar::diffusionRate(
          state.at(layout.z[d], mode, lane),
          observables.effective_diffusivity[d][i],
          A,
          B,
          observables.molar_flux[d][i]);
      }
    }
  }
}

struct SpmPipelineLayout
{
  SpmStateLayout spm{};
  ThermalLumpedLayout thermal{};
  StressHistoryLayout stress_history{};
  StateSlice elapsed_time{};
  StateSlice charge_throughput{};
  StateSlice energy_throughput{};
};

template <int NCH>
struct SpmPipelineParams
{
  bool enable_stress{ true };
  bool enable_sei{ true };
  bool enable_surface_crack{ true };
  bool enable_lam{ true };
  bool enable_lithium_plating{ true };
  SpmDiffusionRhsParams<NCH> diffusion{};
  SpmElectricalParams<NCH> electrical{};
  ThermalLumpedParams thermal{};
  SpmStressParams<NCH> stress{};
  SeiRhsParams<NCH> sei{};
  SurfaceCrackRhsParams<NCH> surface_crack{};
  LamParams lam{};
  LithiumPlatingRhsParams<NCH> lithium_plating{};
};

struct EmptyPipelineScratch
{
  explicit EmptyPipelineScratch(int) {}
};

/**
 * A fixed pipeline is one batch archetype. Optional mechanisms are compile-time branches,
 * while every hot evaluation has one observable pass, one mandatory derivative clear, and no
 * allocation. Runtime model selection therefore happens once, when the batch is constructed.
 */
template <int NCH,
          bool WithThermal,
          bool WithSei,
          bool WithSurfaceCrack,
          bool WithLam,
          bool WithLithiumPlating>
class SpmPipeline
{
public:
  static constexpr bool needs_stress = WithSurfaceCrack || WithLam;
  static constexpr bool needs_sei_scratch = WithSei || WithSurfaceCrack;
  static constexpr bool needs_full_rhs_observables = WithThermal || WithSei
                                                     || WithSurfaceCrack || WithLam
                                                     || WithLithiumPlating;
  static constexpr bool supports_fused_euler = !needs_full_rhs_observables;

  static SpmPipelineLayout declareLayout(BatchBuilder &builder)
  {
    SpmPipelineLayout layout;
    layout.spm = declareSpmState<NCH>(builder,
                                      WithThermal,
                                      WithSei || WithSurfaceCrack || WithLam
                                        || WithLithiumPlating);
    if constexpr (WithThermal)
      layout.thermal = declareThermalLumped(builder, layout.spm.temperature);
    if constexpr (needs_stress)
      layout.stress_history = declareStressHistory(builder);
    layout.elapsed_time = builder.declare({ "elapsed_time", 1, Unit::s, StateRole::cumulative });
    layout.charge_throughput = builder.declare({ "charge_throughput", 1, Unit::Ah, StateRole::cumulative });
    layout.energy_throughput = builder.declare({ "energy_throughput", 1, Unit::Wh, StateRole::cumulative });
    return layout;
  }

  SpmPipeline(SpmPipelineParams<NCH> params,
              SpmPipelineLayout layout,
              int n_lanes)
    : params_{ std::move(params) }, layout_{ layout },
      n_lanes_{ checkedLaneCount(n_lanes) },
      observables_{ n_lanes_ }, stress_{ n_lanes_ }, sei_{ n_lanes_ },
      surface_crack_{ n_lanes_ }, lam_{ n_lanes_ },
      transport_cache_{ n_lanes_ },
      single_observables_{ 1 }, single_transport_cache_{ 1 },
      thevenin_current_density_(static_cast<std::size_t>(n_lanes_)),
      slow_rate_scratch_(slowRateScratchSize(layout_, n_lanes_)),
      plating_{ n_lanes_ }
  {
    assert(n_lanes > 0);
  }

  [[nodiscard]] slide::Status evaluate(RhsViews &views, const StepCtx &ctx)
  {
    assert(views.y.n_lanes() == n_lanes_ && views.ydot.n_lanes() == n_lanes_);
    views.zero_derivative();

    auto observable_view = observables_.view();
    slide::Status status = slide::Status::Success;
    if constexpr (needs_full_rhs_observables) {
      status = computeSpmObservables(params_.electrical,
                                     views.y,
                                     layout_.spm,
                                     ctx,
                                     observable_view,
                                     &transport_cache_);
      if (status != slide::Status::Success)
        return status;
    } else {
      computeSpmTransport(params_.electrical.concentration,
                          views.y,
                          layout_.spm,
                          ctx,
                          observable_view.effective_diffusivity,
                          observable_view.molar_flux,
                          &transport_cache_);
    }

    BasicSpmStress<real_t> stress_view;
    if constexpr (needs_stress) {
      computeSpmConcentrations(params_.electrical.concentration,
                               views.y,
                               layout_.spm,
                               ctx,
                               observable_view.concentration,
                               observable_view.effective_diffusivity,
                               observable_view.molar_flux);
      if (params_.enable_stress) {
        stress_view = stress_.view();
        status = computeSpmStress(params_.stress,
                                  observable_view,
                                  n_lanes_,
                                  stress_view);
        if (status != slide::Status::Success)
          return status;
      }
    }

    addSpmDiffusionRhs(params_.diffusion,
                       views.y,
                       views.ydot,
                       layout_.spm,
                       observable_view);

    if constexpr (WithThermal) {
      status = addThermalLumpedRhs(params_.thermal,
                                   views.y,
                                   views.ydot,
                                   layout_.thermal,
                                   observable_view,
                                   ctx);
      if (status != slide::Status::Success)
        return status;
    }

    BasicSeiOutput<real_t> sei_view;
    if constexpr (needs_sei_scratch)
      sei_view = sei_.view();
    if constexpr (WithSei) {
      status = detail::evaluate_ageing_stage(
        params_.enable_sei,
        sei_view,
        [&](BasicSeiOutput<real_t> output) {
          return computeSei(params_.sei.mechanism,
                            views.y,
                            layout_.spm,
                            ctx,
                            observable_view,
                            output);
        },
        [&](const BasicSeiOutput<real_t> &output) {
          addSeiRhs(params_.sei,
                    views.y,
                    views.ydot,
                    layout_.spm,
                    output);
        });
      if (status != slide::Status::Success)
        return status;
    }
    if constexpr (WithSurfaceCrack) {
      bool sei_enabled = false;
      if constexpr (WithSei)
        sei_enabled = params_.enable_sei;
      if (!sei_enabled)
        detail::clear_ageing_fields<real_t, 2>(
          n_lanes_,
          { sei_view.side_reaction_current, sei_view.active_fraction_rate });
    }

    if constexpr (WithSurfaceCrack) {
      auto crack_view = surface_crack_.view();
      status = detail::evaluate_ageing_stage(
        params_.enable_surface_crack,
        crack_view,
        [&](BasicSurfaceCrackOutput<real_t> output) {
          return computeSurfaceCrack(params_.surface_crack.mechanism,
                                     views.y,
                                     layout_.spm,
                                     layout_.stress_history,
                                     ctx,
                                     observable_view,
                                     stress_view,
                                     output);
        },
        [&](const BasicSurfaceCrackOutput<real_t> &output) {
          addSurfaceCrackRhs(params_.surface_crack,
                             views.y,
                             views.ydot,
                             layout_.spm,
                             sei_view,
                             output);
        });
      if (status != slide::Status::Success)
        return status;
    }

    if constexpr (WithLam) {
      auto lam_view = lam_.view();
      status = detail::evaluate_ageing_stage(
        params_.enable_lam,
        lam_view,
        [&](BasicLamOutput<real_t> output) {
          return computeLam(params_.lam,
                            views.y,
                            layout_.spm,
                            layout_.stress_history,
                            ctx,
                            observable_view,
                            stress_view,
                            output);
        },
        [&](const BasicLamOutput<real_t> &output) {
          addLamRhs(params_.lam, views.ydot, layout_.spm, output);
        });
      if (status != slide::Status::Success)
        return status;
    }

    if constexpr (WithLithiumPlating) {
      auto plating_view = plating_.view();
      status = detail::evaluate_ageing_stage(
        params_.enable_lithium_plating,
        plating_view,
        [&](BasicLithiumPlatingOutput<real_t> output) {
          return computeLithiumPlating(params_.lithium_plating.mechanism,
                                       views.y,
                                       layout_.spm,
                                       ctx,
                                       observable_view,
                                       output);
        },
        [&](const BasicLithiumPlatingOutput<real_t> &output) {
          addLithiumPlatingRhs(
            params_.lithium_plating,
            views.y,
            views.ydot,
            layout_.spm,
            BasicLithiumPlatingOutput<const real_t>{
              output.side_reaction_current });
        });
      if (status != slide::Status::Success)
        return status;
    }
    return slide::Status::Success;
  }

  /** Fuse transport, diffusion RHS, and Euler update for the base isothermal archetype. */
  [[nodiscard]] slide::Status advanceEuler(BatchView state,
                                           const StepCtx &ctx,
                                           real_t dt,
                                           std::span<real_t> terminal_voltage)
    requires(supports_fused_euler)
  {
    assert(static_cast<int>(terminal_voltage.size()) == n_lanes_);
    const ConstBatchView current{ state.shape(), std::span<const real_t>{ state.raw() } };
    const int period = lanePeriod(current, ctx);
    const bool coalesced = period < n_lanes_;
    auto observable_view = period == 1 ? single_observables_.view()
                                       : observables_.view(period);
    const BatchShape evaluation_shape = coalesced
                                          ? BatchShape{ state.n_rows(), period, state.stride() }
                                          : state.shape();
    const ConstBatchView evaluation_state{ evaluation_shape,
                                           std::span<const real_t>{ state.raw() } };
    const StepCtx evaluation_ctx{ .time = ctx.time,
                                  .dt = ctx.dt,
                                  .i_app = coalesced ? ctx.i_app.first(static_cast<std::size_t>(period))
                                                     : ctx.i_app };
    auto *cache = period == 1          ? &single_transport_cache_
                  : period == n_lanes_ ? &transport_cache_
                                       : nullptr;
    computeSpmTransport(params_.electrical.concentration,
                        evaluation_state,
                        layout_.spm,
                        evaluation_ctx,
                        observable_view.effective_diffusivity,
                        observable_view.molar_flux,
                        cache);
    const int evaluated_lanes = period;
    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      for (int mode = 0; mode < NCH; ++mode) {
        const real_t A = params_.diffusion.A[d][static_cast<std::size_t>(mode)];
        const real_t B = params_.diffusion.B[d][static_cast<std::size_t>(mode)];
        auto state_row = state.row(layout_.spm.z[d].row_begin + mode);
        for (int lane = 0; lane < evaluated_lanes; ++lane) {
          const auto i = static_cast<std::size_t>(lane);
          auto &z = state_row[i];
          z += dt * spm_scalar::diffusionRate(z,
                                               observable_view.effective_diffusivity[d][i],
                                               A,
                                               B,
                                               observable_view.molar_flux[d][i]);
        }
        if (coalesced)
          for (int lane = period; lane < n_lanes_; ++lane)
            state_row[static_cast<std::size_t>(lane)] = state_row[static_cast<std::size_t>(lane % period)];
      }
    }

    const ConstBatchView accepted_state{ evaluation_shape,
                                         std::span<const real_t>{ state.raw() } };
    const StepCtx accepted_ctx{ .time = ctx.time + dt,
                                .dt = 0.0,
                                .i_app = evaluation_ctx.i_app };
    auto status = computeSpmObservables(params_.electrical,
                                        accepted_state,
                                        layout_.spm,
                                        accepted_ctx,
                                        observable_view,
                                        cache);
    if (status != slide::Status::Success)
      return status;
    if (coalesced) {
      for (int lane = 0; lane < n_lanes_; ++lane)
        terminal_voltage[static_cast<std::size_t>(lane)] =
          observable_view.terminal_voltage[static_cast<std::size_t>(lane % period)];
    } else {
      std::copy(observable_view.terminal_voltage.begin(),
                observable_view.terminal_voltage.end(),
                terminal_voltage.begin());
    }
    return status;
  }

  /** Exact diagonal modal diffusion wrapped in a symmetric slow/fast/slow split. */
  [[nodiscard]] slide::Status advanceExponential(
    BatchView state,
    BatchView derivative,
    const StepCtx &ctx,
    real_t dt,
    std::span<real_t> terminal_voltage)
  {
    if (state.n_lanes() != n_lanes_ || derivative.n_lanes() != n_lanes_
        || state.n_rows() != derivative.n_rows()
        || static_cast<int>(terminal_voltage.size()) != n_lanes_)
      return slide::Status::Invalid_parameters;

    auto evaluate_slow = [&](BatchView output) -> slide::Status {
      RhsViews views{ state.shape() };
      views.rebind(state.raw(), output.raw());
      const auto status = evaluate(views, ctx);
      if (status != slide::Status::Success)
        return status;
      const auto observable = observables_.view();
      for (int row = 0; row < output.n_rows(); ++row) {
        auto rates = output.row(row);
        for (int lane = 0; lane < n_lanes_; ++lane) {
          for (const Domain domain : domains) {
            const auto d = domain_index(domain);
            const int mode = row - layout_.spm.z[d].row_begin;
            if (mode >= 0 && mode < NCH) {
              const auto i = static_cast<std::size_t>(lane);
              rates[i] -= spm_scalar::diffusionRate(
                state.row(row)[i],
                observable.effective_diffusivity[d][i],
                params_.diffusion.A[d][static_cast<std::size_t>(mode)],
                params_.diffusion.B[d][static_cast<std::size_t>(mode)],
                observable.molar_flux[d][i]);
            }
          }
        }
      }
      return slide::Status::Success;
    };

    auto apply_slow = [&](real_t h) -> slide::Status {
      if constexpr (!needs_full_rhs_observables) {
        (void)h;
        return slide::Status::Success;
      } else {
        BatchView first_rate{ state.shape(), slow_rate_scratch_ };
        auto status = evaluate_slow(first_rate);
        if (status != slide::Status::Success)
          return status;
        for (int row = 0; row < state.n_rows(); ++row) {
          auto values = state.row(row);
          const auto first = first_rate.row(row);
          for (int lane = 0; lane < n_lanes_; ++lane)
            values[static_cast<std::size_t>(lane)] += h * first[static_cast<std::size_t>(lane)];
        }
        status = evaluate_slow(derivative);
        if (status != slide::Status::Success)
          return status;
        for (int row = 0; row < state.n_rows(); ++row) {
          auto values = state.row(row);
          const auto first = first_rate.row(row);
          const auto second = derivative.row(row);
          for (int lane = 0; lane < n_lanes_; ++lane) {
            const auto i = static_cast<std::size_t>(lane);
            values[i] += 0.5 * h * (second[i] - first[i]);
          }
        }
        return slide::Status::Success;
      }
    };

    auto status = apply_slow(0.5 * dt);
    if (status != slide::Status::Success)
      return status;

    const ConstBatchView midpoint{ state.shape(),
                                   std::span<const real_t>{ state.raw() } };
    auto observable = observables_.view();
    computeSpmTransport(params_.electrical.concentration, midpoint, layout_.spm, ctx, observable.effective_diffusivity, observable.molar_flux, &transport_cache_);
    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      for (int mode = 0; mode < NCH; ++mode) {
        const real_t eigenvalue = params_.diffusion.A[d][static_cast<std::size_t>(mode)];
        const real_t input = params_.diffusion.B[d][static_cast<std::size_t>(mode)];
        auto values = state.row(layout_.spm.z[d].row_begin + mode);
        for (int lane = 0; lane < n_lanes_; ++lane) {
          const auto i = static_cast<std::size_t>(lane);
          SLIDE_SPM_ADVANCE_MODAL_STD(values[i],
                                      observable.effective_diffusivity[d][i],
                                      eigenvalue,
                                      dt,
                                      input,
                                      observable.molar_flux[d][i]);
        }
      }
    }

    status = apply_slow(0.5 * dt);
    if (status != slide::Status::Success)
      return status;
    const ConstBatchView accepted{ state.shape(),
                                   std::span<const real_t>{ state.raw() } };
    const StepCtx accepted_ctx{ .time = ctx.time + dt,
                                .dt = 0.0,
                                .i_app = ctx.i_app };
    return observeTerminalVoltage(accepted, accepted_ctx, terminal_voltage);
  }

  [[nodiscard]] slide::Status observeTerminalVoltage(
    const ConstBatchView &state,
    const StepCtx &ctx,
    std::span<real_t> terminal_voltage)
  {
    if (state.n_lanes() != n_lanes_
        || static_cast<int>(terminal_voltage.size()) != n_lanes_)
      return slide::Status::Invalid_parameters;
    auto observable_view = observables_.view();
    const auto status = computeSpmObservables(params_.electrical,
                                              state,
                                              layout_.spm,
                                              ctx,
                                              observable_view,
                                              &transport_cache_);
    if (status != slide::Status::Success)
      return status;
    std::copy(observable_view.terminal_voltage.begin(),
              observable_view.terminal_voltage.end(),
              terminal_voltage.begin());
    return slide::Status::Success;
  }

  /** Frozen-state tangent V(I) = intercept - resistance*I for the pack solver. */
  [[nodiscard]] slide::Status linearizeThevenin(
    const ConstBatchView &state,
    std::span<const real_t> current,
    std::span<real_t> intercept_ocv,
    std::span<real_t> resistance)
  {
    if (state.n_lanes() != n_lanes_
        || static_cast<int>(current.size()) != n_lanes_
        || intercept_ocv.size() != current.size()
        || resistance.size() != current.size()
        || !(is_finite(params_.electrical.electrode_area)
             && params_.electrical.electrode_area > 0.0))
      return slide::Status::Invalid_parameters;

    const int period = lanePeriod(state, current);
    for (int lane = 0; lane < period; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      if (!is_finite(current[i]))
        return slide::Status::Invalid_parameters;
      thevenin_current_density_[i] = current[i] / params_.electrical.electrode_area;
    }
    const BatchShape evaluation_shape = period < n_lanes_
                                          ? BatchShape{ state.n_rows(), period, state.stride() }
                                          : state.shape();
    const ConstBatchView evaluation_state{ evaluation_shape, state.raw() };
    StepCtx ctx{ .time = 0.0,
                 .dt = 0.0,
                 .i_app = std::span<const real_t>{ thevenin_current_density_ }.first(
                   static_cast<std::size_t>(period)) };
    auto observable = observables_.view(period);
    auto *cache = period == n_lanes_ ? &transport_cache_ : nullptr;
    auto status = computeSpmObservables(params_.electrical, evaluation_state, layout_.spm, ctx, observable, cache);
    if (status != slide::Status::Success)
      return status;

    const auto &p = params_.electrical;
    for (int lane = 0; lane < period; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const real_t temperature = state.at(layout_.spm.temperature, 0, lane);
      PerDomain<real_t> concentration_slope{};
      PerDomain<real_t> stoichiometry_slope{};
      PerDomain<real_t> overpotential_slope{};
      PerDomain<real_t> ocv_slope{};
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const auto &electrode = p.electrode[d];
        const real_t area = state.at(layout_.spm.specific_surface_area[d], 0, lane);
        const real_t thickness = state.at(layout_.spm.electrode_thickness[d], 0, lane);
        const real_t dflux_dcurrent = static_cast<real_t>(molar_flux_sign(domain))
                                      / (p.electrode_area * area * p.n * p.F * thickness);
        concentration_slope[d] = p.concentration.Dout[d][0] * dflux_dcurrent
                                 / observable.effective_diffusivity[d][i];
        stoichiometry_slope[d] = concentration_slope[d] / electrode.cs_max;
        ocv_slope[d] = p.electrode_ocv[d].derivative(
                         observable.surface_stoichiometry[d][i])
                       * stoichiometry_slope[d];

        const real_t exchange = observable.exchange_current_density[d][i];
        const real_t cs = observable.concentration[d][i];
        const real_t exchange_log_slope = 0.5 * concentration_slope[d]
                                          * (1.0 / cs
                                             - 1.0 / (electrode.cs_max - cs));
        const real_t argument_factor = 0.5 * static_cast<real_t>(molar_flux_sign(domain))
                                       / (p.electrode_area * area * thickness);
        const real_t argument = argument_factor * current[i] / exchange;
        const real_t argument_slope = argument_factor / exchange
                                      - argument * exchange_log_slope;
        overpotential_slope[d] = 2.0 * p.Rg * temperature / (p.n * p.F)
                                 * argument_slope / std::sqrt(1.0 + argument * argument);
      }
      const auto neg = domain_index(Domain::neg);
      const auto pos = domain_index(Domain::pos);
      const real_t entropic_slope = p.total_entropic_coefficient.derivative(
                                      observable.surface_stoichiometry[pos][i])
                                    * stoichiometry_slope[pos];
      const real_t voltage_slope = ocv_slope[pos] - ocv_slope[neg]
                                   + (temperature - p.reference_temperature) * entropic_slope
                                   + overpotential_slope[pos] - overpotential_slope[neg]
                                   - observable.resistance[i];
      const real_t r = -voltage_slope;
      const real_t intercept = observable.terminal_voltage[i] + r * current[i];
      if (!is_finite(r) || !(r > 0.0) || !is_finite(intercept))
        return slide::Status::Invalid_states;
      resistance[i] = r;
      intercept_ocv[i] = intercept;
    }
    for (int lane = period; lane < n_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const auto source = static_cast<std::size_t>(lane % period);
      resistance[i] = resistance[source];
      intercept_ocv[i] = intercept_ocv[source];
    }
    return slide::Status::Success;
  }

  /** Store accepted-step stress values in the checkpointed algebraic history rows. */
  void storeStressHistory(BatchView state, real_t interval)
  {
    if constexpr (needs_stress) {
      if (!params_.enable_stress)
        return;
      assert(state.n_lanes() == n_lanes_ && interval > 0.0);
      const auto stress_view = stress_.view();
      for (int lane = 0; lane < n_lanes_; ++lane) {
        const auto i = static_cast<std::size_t>(lane);
        for (const Domain domain : domains)
          state.at(layout_.stress_history.previous_dai[domain_index(domain)], 0, lane) = stress_view.dai_maximum_hydrostatic[domain_index(domain)][i];
        state.at(layout_.stress_history.previous_laresgoiti_negative, 0, lane) = stress_view.laresgoiti_negative[i];
        state.at(layout_.stress_history.interval, 0, lane) = interval;
      }
    } else {
      (void)state;
      (void)interval;
    }
  }

  const SpmPipelineLayout &layout() const { return layout_; }
  const SpmPipelineParams<NCH> &params() const { return params_; }
  int n_lanes() const { return n_lanes_; }
  int trustedLanePeriod() const { return trusted_lane_period_; }

  [[nodiscard]] slide::Status setTrustedLanePeriod(const ConstBatchView &state,
                                                   int maximum_period)
  {
    if (state.n_lanes() != n_lanes_ || maximum_period <= 0
        || n_lanes_ % maximum_period != 0)
      return slide::Status::Invalid_parameters;
    for (int period = 1; period <= maximum_period; ++period) {
      if (maximum_period % period != 0)
        continue;
      bool equal = true;
      for (int row = 0; row < state.n_rows() && equal; ++row) {
        const auto values = state.row(row);
        for (int lane = period; lane < n_lanes_; ++lane)
          if (values[static_cast<std::size_t>(lane)]
              != values[static_cast<std::size_t>(lane % period)]) {
            equal = false;
            break;
          }
      }
      if (equal) {
        trusted_lane_period_ = period;
        return slide::Status::Success;
      }
    }
    return slide::Status::Invalid_states;
  }

private:
  static int checkedLaneCount(int lanes)
  {
    if (lanes <= 0)
      throw std::invalid_argument{ "SPM pipeline requires at least one lane" };
    constexpr int padding = 7;
    if (lanes > std::numeric_limits<int>::max() - padding)
      throw std::length_error{ "SPM pipeline padded lane count is not representable" };
    return lanes;
  }

  static std::size_t slowRateScratchSize(const SpmPipelineLayout &layout,
                                         int lanes)
  {
    if constexpr (!needs_full_rhs_observables)
      return 0;

    const int row_begin = layout.elapsed_time.row_begin;
    const int rows = layout.elapsed_time.rows;
    if (row_begin < 0 || rows <= 0)
      throw std::invalid_argument{ "SPM pipeline elapsed-time layout is invalid" };

    constexpr std::size_t extra_rows = 2;
    constexpr std::size_t maximum = std::numeric_limits<std::size_t>::max();
    const auto begin = static_cast<std::size_t>(row_begin);
    const auto count = static_cast<std::size_t>(rows);
    if (count > maximum - extra_rows
        || begin > maximum - count - extra_rows)
      throw std::length_error{ "SPM pipeline scratch row count is not representable" };
    const std::size_t row_count = begin + count + extra_rows;

    constexpr std::size_t block = 8;
    const auto lane_count = static_cast<std::size_t>(lanes);
    const std::size_t padded_lanes = ((lane_count + block - 1) / block) * block;
    if (row_count > maximum / padded_lanes)
      throw std::length_error{ "SPM pipeline scratch extent is not representable" };
    const std::size_t elements = row_count * padded_lanes;
    if (elements > maximum / sizeof(real_t))
      throw std::length_error{ "SPM pipeline scratch byte count is not representable" };
    return elements;
  }

  int lanePeriod(const ConstBatchView &state,
                 std::span<const real_t>
                   current) const
  {
    if (n_lanes_ <= 1)
      return n_lanes_;
    if (trusted_lane_period_ < n_lanes_) {
      for (int lane = trusted_lane_period_; lane < n_lanes_; ++lane)
        if (current[static_cast<std::size_t>(lane)]
            != current[static_cast<std::size_t>(lane % trusted_lane_period_)])
          return n_lanes_;
      return trusted_lane_period_;
    }
    for (int period = 1; period < n_lanes_; ++period) {
      if (n_lanes_ % period != 0)
        continue;
      bool equal = true;
      for (int lane = period; lane < n_lanes_ && equal; ++lane)
        equal = current[static_cast<std::size_t>(lane)]
                == current[static_cast<std::size_t>(lane % period)];
      for (int row = 0; row < state.n_rows() && equal; ++row) {
        const auto values = state.row(row);
        for (int lane = period; lane < n_lanes_; ++lane)
          if (values[static_cast<std::size_t>(lane)]
              != values[static_cast<std::size_t>(lane % period)]) {
            equal = false;
            break;
          }
      }
      if (equal)
        return period;
    }
    return n_lanes_;
  }

  int lanePeriod(const ConstBatchView &state, const StepCtx &ctx) const
  {
    return lanePeriod(state, ctx.i_app);
  }

  SpmPipelineParams<NCH> params_;
  SpmPipelineLayout layout_{};
  int n_lanes_{};
  int trusted_lane_period_{ n_lanes_ };
  SpmObservableScratch<NCH> observables_;
  std::conditional_t<needs_stress, SpmStressScratch<>, EmptyPipelineScratch> stress_;
  std::conditional_t<needs_sei_scratch, SeiScratch<>, EmptyPipelineScratch> sei_;
  std::conditional_t<WithSurfaceCrack, SurfaceCrackScratch<>, EmptyPipelineScratch>
    surface_crack_;
  std::conditional_t<WithLam, LamScratch<>, EmptyPipelineScratch> lam_;
  SpmTransportCache transport_cache_;
  SpmObservableScratch<NCH> single_observables_;
  SpmTransportCache single_transport_cache_;
  std::vector<real_t> thevenin_current_density_{};
  std::vector<real_t> slow_rate_scratch_{};
  std::conditional_t<WithLithiumPlating,
                     LithiumPlatingScratch<>,
                     EmptyPipelineScratch>
    plating_;
};

} // namespace slide::core
