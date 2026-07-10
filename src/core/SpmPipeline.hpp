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

  static SpmPipelineLayout declareLayout(BatchBuilder &builder)
  {
    SpmPipelineLayout layout;
    layout.spm = declareSpmState<NCH>(builder);
    if constexpr (WithThermal)
      layout.thermal = declareThermalLumped(builder, layout.spm.temperature);
    if constexpr (needs_stress)
      layout.stress_history = declareStressHistory(builder);
    return layout;
  }

  SpmPipeline(SpmPipelineParams<NCH> params,
              SpmPipelineLayout layout,
              int n_lanes)
    : params_{ std::move(params) }, layout_{ layout }, n_lanes_{ n_lanes },
      observables_{ n_lanes }, stress_{ n_lanes }, sei_{ n_lanes },
      surface_crack_{ n_lanes }, lam_{ n_lanes },
      plating_current_(WithLithiumPlating ? static_cast<std::size_t>(n_lanes) : 0)
  {
    assert(n_lanes > 0);
  }

  [[nodiscard]] slide::Status evaluate(RhsViews &views, const StepCtx &ctx)
  {
    assert(views.y.n_lanes() == n_lanes_ && views.ydot.n_lanes() == n_lanes_);
    views.zero_derivative();

    auto observable_view = observables_.view();
    auto status = computeSpmObservables(params_.electrical,
                                        views.y,
                                        layout_.spm,
                                        ctx,
                                        observable_view);
    if (status != slide::Status::Success)
      return status;

    BasicSpmStress<real_t> stress_view;
    if constexpr (needs_stress) {
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
  SpmPipelineParams<NCH> params_;
  SpmPipelineLayout layout_{};
  int n_lanes_{};
  SpmObservableScratch<NCH> observables_;
  std::conditional_t<needs_stress, SpmStressScratch<>, EmptyPipelineScratch> stress_;
  std::conditional_t<needs_sei_scratch, SeiScratch<>, EmptyPipelineScratch> sei_;
  std::conditional_t<WithSurfaceCrack, SurfaceCrackScratch<>, EmptyPipelineScratch>
    surface_crack_;
  std::conditional_t<WithLam, LamScratch<>, EmptyPipelineScratch> lam_;
  std::vector<real_t> plating_current_{};
};

} // namespace slide::core
