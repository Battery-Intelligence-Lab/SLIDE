/**
 * @file SpmPipeline.hpp
 * @brief Compile-time composed, allocation-free SPM RHS pipeline.
 */

#pragma once

#include "Lam.hpp"
#include "LithiumPlating.hpp"
#include "Sei.hpp"
#include "SurfaceCrack.hpp"
#include "ThermalLumped.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <span>
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
        derivative.at(layout.z[d], mode, lane) += observables.effective_diffusivity[d][i] * A
                                                    * state.at(layout.z[d], mode, lane)
                                                  + B * observables.molar_flux[d][i];
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
    : params_{ std::move(params) }, layout_{ layout }, n_lanes_{ n_lanes },
      observables_{ n_lanes }, stress_{ n_lanes }, sei_{ n_lanes },
      surface_crack_{ n_lanes }, lam_{ n_lanes },
      transport_cache_{ n_lanes },
      single_observables_{ 1 }, single_transport_cache_{ 1 },
      thevenin_current_density_(static_cast<std::size_t>(n_lanes)),
      thevenin_voltage_(static_cast<std::size_t>(n_lanes)),
      thevenin_voltage_plus_(static_cast<std::size_t>(n_lanes)),
      thevenin_voltage_minus_(static_cast<std::size_t>(n_lanes)),
      thevenin_step_(static_cast<std::size_t>(n_lanes)),
      plating_current_(WithLithiumPlating ? static_cast<std::size_t>(n_lanes) : 0)
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
      if (params_.enable_sei) {
        status = computeSei(params_.sei.mechanism,
                            views.y,
                            layout_.spm,
                            ctx,
                            observable_view,
                            sei_view);
        if (status != slide::Status::Success)
          return status;
        addSeiRhs(params_.sei,
                  views.y,
                  views.ydot,
                  layout_.spm,
                  sei_view);
      } else if constexpr (WithSurfaceCrack) {
        std::fill(sei_view.side_reaction_current.begin(),
                  sei_view.side_reaction_current.end(),
                  real_t{});
        std::fill(sei_view.active_fraction_rate.begin(),
                  sei_view.active_fraction_rate.end(),
                  real_t{});
      }
    } else if constexpr (WithSurfaceCrack) {
      std::fill(sei_view.side_reaction_current.begin(),
                sei_view.side_reaction_current.end(),
                real_t{});
      std::fill(sei_view.active_fraction_rate.begin(),
                sei_view.active_fraction_rate.end(),
                real_t{});
    }

    if constexpr (WithSurfaceCrack) {
      if (params_.enable_surface_crack) {
        auto crack_view = surface_crack_.view();
        status = computeSurfaceCrack(params_.surface_crack.mechanism,
                                     views.y,
                                     layout_.spm,
                                     layout_.stress_history,
                                     ctx,
                                     observable_view,
                                     stress_view,
                                     crack_view);
        if (status != slide::Status::Success)
          return status;
        addSurfaceCrackRhs(params_.surface_crack,
                           views.y,
                           views.ydot,
                           layout_.spm,
                           sei_view,
                           crack_view);
      }
    }

    if constexpr (WithLam) {
      if (params_.enable_lam) {
        auto lam_view = lam_.view();
        status = computeLam(params_.lam,
                            views.y,
                            layout_.spm,
                            layout_.stress_history,
                            ctx,
                            observable_view,
                            stress_view,
                            lam_view);
        if (status != slide::Status::Success)
          return status;
        addLamRhs(params_.lam, views.ydot, layout_.spm, lam_view);
      }
    }

    if constexpr (WithLithiumPlating) {
      if (params_.enable_lithium_plating) {
        status = computeLithiumPlating(params_.lithium_plating.mechanism,
                                       views.y,
                                       layout_.spm,
                                       ctx,
                                       observable_view,
                                       std::span<real_t>{ plating_current_ });
        if (status != slide::Status::Success)
          return status;
        addLithiumPlatingRhs(params_.lithium_plating,
                             views.y,
                             views.ydot,
                             layout_.spm,
                             std::span<const real_t>{ plating_current_ });
      }
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
    const bool coalesced = lanesEqual(current, ctx);
    auto observable_view = coalesced ? single_observables_.view() : observables_.view();
    const BatchShape evaluation_shape = coalesced
                                          ? BatchShape{ state.n_rows(), 1, state.stride() }
                                          : state.shape();
    const ConstBatchView evaluation_state{ evaluation_shape,
                                           std::span<const real_t>{ state.raw() } };
    const StepCtx evaluation_ctx{ .time = ctx.time,
                                  .dt = ctx.dt,
                                  .i_app = coalesced ? ctx.i_app.first(1) : ctx.i_app };
    auto *cache = coalesced ? &single_transport_cache_ : &transport_cache_;
    computeSpmTransport(params_.electrical.concentration,
                        evaluation_state,
                        layout_.spm,
                        evaluation_ctx,
                        observable_view.effective_diffusivity,
                        observable_view.molar_flux,
                        cache);
    const int evaluated_lanes = coalesced ? 1 : n_lanes_;
    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      for (int mode = 0; mode < NCH; ++mode) {
        const real_t A = params_.diffusion.A[d][static_cast<std::size_t>(mode)];
        const real_t B = params_.diffusion.B[d][static_cast<std::size_t>(mode)];
        auto state_row = state.row(layout_.spm.z[d].row_begin + mode);
        for (int lane = 0; lane < evaluated_lanes; ++lane) {
          const auto i = static_cast<std::size_t>(lane);
          auto &z = state_row[i];
          z += dt * (observable_view.effective_diffusivity[d][i] * A * z + B * observable_view.molar_flux[d][i]);
        }
        if (coalesced)
          std::fill(state_row.begin() + 1, state_row.end(), state_row.front());
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
      std::fill(terminal_voltage.begin(), terminal_voltage.end(), observable_view.terminal_voltage.front());
    } else {
      std::copy(observable_view.terminal_voltage.begin(),
                observable_view.terminal_voltage.end(),
                terminal_voltage.begin());
    }
    return status;
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

    constexpr real_t relative_step = 6.0554544523933429e-6; // cbrt(epsilon)
    for (int lane = 0; lane < n_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      if (!is_finite(current[i]))
        return slide::Status::Invalid_parameters;
      thevenin_step_[i] = relative_step * std::max(real_t{ 1 }, std::abs(current[i]));
      thevenin_current_density_[i] = current[i] / params_.electrical.electrode_area;
    }
    StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = thevenin_current_density_ };
    auto status = observeTerminalVoltage(state, ctx, thevenin_voltage_);
    if (status != slide::Status::Success)
      return status;
    for (int lane = 0; lane < n_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      thevenin_current_density_[i] = (current[i] + thevenin_step_[i])
                                     / params_.electrical.electrode_area;
    }
    status = observeTerminalVoltage(state, ctx, thevenin_voltage_plus_);
    if (status != slide::Status::Success)
      return status;
    for (int lane = 0; lane < n_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      thevenin_current_density_[i] = (current[i] - thevenin_step_[i])
                                     / params_.electrical.electrode_area;
    }
    status = observeTerminalVoltage(state, ctx, thevenin_voltage_minus_);
    if (status != slide::Status::Success)
      return status;
    for (int lane = 0; lane < n_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const real_t tangent = (thevenin_voltage_plus_[i] - thevenin_voltage_minus_[i])
                             / (2.0 * thevenin_step_[i]);
      const real_t r = -tangent;
      const real_t intercept = thevenin_voltage_[i] + r * current[i];
      if (!is_finite(r) || !(r > 0.0) || !is_finite(intercept))
        return slide::Status::Invalid_states;
      resistance[i] = r;
      intercept_ocv[i] = intercept;
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

private:
  bool lanesEqual(const ConstBatchView &state, const StepCtx &ctx) const
  {
    if (n_lanes_ <= 1)
      return true;
    for (int lane = 1; lane < n_lanes_; ++lane)
      if (ctx.i_app[static_cast<std::size_t>(lane)] != ctx.i_app.front())
        return false;
    for (int row = 0; row < state.n_rows(); ++row) {
      const auto values = state.row(row);
      for (int lane = 1; lane < n_lanes_; ++lane)
        if (values[static_cast<std::size_t>(lane)] != values.front())
          return false;
    }
    return true;
  }

  SpmPipelineParams<NCH> params_;
  SpmPipelineLayout layout_{};
  int n_lanes_{};
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
  std::vector<real_t> thevenin_voltage_{};
  std::vector<real_t> thevenin_voltage_plus_{};
  std::vector<real_t> thevenin_voltage_minus_{};
  std::vector<real_t> thevenin_step_{};
  std::vector<real_t> plating_current_{};
};

} // namespace slide::core
