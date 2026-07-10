/**
 * @file SpmObservables.hpp
 * @brief v4 core OBSERVABLE-RECONSTRUCTION layer for the SPM (PLAN.md §3.7/§3.12, D-10, D-23).
 *
 * "Observables derived from state, one code path shared by the RHS internals and the Recorder."
 * These are the quantities the legacy Cell_SPM recomputes on demand from its 29 states — surface
 * concentration, then (built on it) overpotential / OCV / terminal voltage. Only time/Ah/Wh are
 * path-dependent; everything here is a pure function of the arena state + operating point, so it
 * is reconstructed lazily instead of stored (D-10).
 *
 * The first stage reconstructs the complete particle concentration profile at
 * {surface, NCH interior nodes, centre}. The surface uses the Model_SPM C/D output map and the
 * centre uses its independent Cc/cc_coeff derivative identity (the historical nch!=5 bug path,
 * §2.2). Kinetics, OCV, resistance, terminal voltage and heat then compose on the surface row.
 *
 * SHAPE (mirrors SpectralDiffusion.hpp, deliberately): batch-shared cold-built params in a struct
 * of designated-initialisable named fields (D-15, no loose 15-arg call); a free-function kernel
 * over a rebindable BatchView and StepCtx (D-23); and a lane sweep over whole SoA rows (PC-3).
 * Output scratch is caller-owned and allocated once by the composed batch (PC-1).
 *
 * @date 2026-07-10
 */

#pragma once

#include "BatchView.hpp"
#include "CellDesign.hpp"
#include "CompiledCurve.hpp"
#include "SpmState.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>
#include <type_traits>
#include <vector>

namespace slide::core {

/**
 * Batch-shared parameters for the concentration observable of one SPM archetype.
 * Indexed by the explicit v4 Domain (negative first). Cold-built once (by the test/factory);
 * immutable in the hot loop. C/Dout contain the surface and interior Model_SPM output rows;
 * Cc/cc_coeff reconstruct the centre node from those outputs.
 */
template <int NCH>
struct SpmConcentrationParams
{
  //!< Physical constants — mirror src/settings/constants.hpp PhyConst (F=96487, Rg=8.314, n=1).
  real_t F{ 96487.0 }; //!< Faraday constant            [C mol⁻¹]
  real_t Rg{ 8.314 };  //!< ideal gas constant          [J mol⁻¹ K⁻¹]
  real_t n{ 1.0 };     //!< electrons in main reaction    [-]
  real_t T_ref{};      //!< Arrhenius reference temperature [K]

  PerDomain<std::array<std::array<real_t, NCH>, NCH + 1>> C{};
  PerDomain<std::array<real_t, NCH + 1>> Dout{};
  std::array<real_t, NCH + 1> Cc{}; //!< centre-node derivative map over surface+interior values
  real_t cc_coeff{};
  PerDomain<real_t> R{};   //!< particle radius [m]
  PerDomain<real_t> D_T{}; //!< electrode Arrhenius activation for D
};

/** Derived transport cache. Exact state-value keys make it restart-safe and self-invalidating. */
struct SpmTransportCache
{
  explicit SpmTransportCache(int lanes)
    : lanes_{ lanes }, valid_(static_cast<std::size_t>(2 * lanes)),
      temperature_(static_cast<std::size_t>(2 * lanes)),
      diffusion_reference_(static_cast<std::size_t>(2 * lanes)),
      specific_area_(static_cast<std::size_t>(2 * lanes)),
      thickness_(static_cast<std::size_t>(2 * lanes)),
      effective_diffusivity_(static_cast<std::size_t>(2 * lanes)),
      flux_denominator_(static_cast<std::size_t>(2 * lanes))
  {}

  int lanes() const { return lanes_; }

  std::size_t index(Domain domain, int lane) const
  {
    return domain_index(domain) * static_cast<std::size_t>(lanes_)
           + static_cast<std::size_t>(lane);
  }

  int lanes_{};
  std::vector<unsigned char> valid_{};
  std::vector<real_t> temperature_{};
  std::vector<real_t> diffusion_reference_{};
  std::vector<real_t> specific_area_{};
  std::vector<real_t> thickness_{};
  std::vector<real_t> effective_diffusivity_{};
  std::vector<real_t> flux_denominator_{};
};

template <int NCH, class Real>
void computeSpmTransportLane(const SpmConcentrationParams<NCH> &p,
                             const BasicBatchView<const Real> &state,
                             const SpmStateLayout &layout,
                             const BasicStepCtx<Real> &ctx,
                             Domain domain,
                             int lane,
                             Real &effective_diffusivity,
                             Real &molar_flux,
                             SpmTransportCache *cache)
{
  const auto d = domain_index(domain);
  const Real temperature = state.at(layout.temperature, 0, lane);
  const Real diffusion_reference = state.at(layout.diffusion_coefficient[d], 0, lane);
  const Real specific_area = state.at(layout.specific_surface_area[d], 0, lane);
  const Real thickness = state.at(layout.electrode_thickness[d], 0, lane);
  if constexpr (std::is_same_v<std::remove_cv_t<Real>, real_t>) {
    if (cache != nullptr) {
      assert(cache->lanes() == state.n_lanes());
      const auto i = cache->index(domain, lane);
      if (cache->valid_[i] != 0 && cache->temperature_[i] == temperature
          && cache->diffusion_reference_[i] == diffusion_reference
          && cache->specific_area_[i] == specific_area
          && cache->thickness_[i] == thickness) {
        effective_diffusivity = cache->effective_diffusivity_[i];
        molar_flux = static_cast<Real>(molar_flux_sign(domain)) * ctx.i_app[lane]
                     / cache->flux_denominator_[i];
        return;
      }
      using std::exp;
      const Real arrhenius = (Real{ 1 } / p.T_ref - Real{ 1 } / temperature) / p.Rg;
      effective_diffusivity = diffusion_reference * exp(p.D_T[d] * arrhenius);
      const Real denominator = specific_area * p.n * p.F * thickness;
      molar_flux = static_cast<Real>(molar_flux_sign(domain)) * ctx.i_app[lane]
                   / denominator;
      cache->temperature_[i] = temperature;
      cache->diffusion_reference_[i] = diffusion_reference;
      cache->specific_area_[i] = specific_area;
      cache->thickness_[i] = thickness;
      cache->effective_diffusivity_[i] = effective_diffusivity;
      cache->flux_denominator_[i] = denominator;
      cache->valid_[i] = 1;
      return;
    }
  }
  using std::exp;
  const Real arrhenius = (Real{ 1 } / p.T_ref - Real{ 1 } / temperature) / p.Rg;
  effective_diffusivity = diffusion_reference * exp(p.D_T[d] * arrhenius);
  const Real denominator = specific_area * p.n * p.F * thickness;
  molar_flux = static_cast<Real>(molar_flux_sign(domain)) * ctx.i_app[lane]
               / denominator;
}

/** Compute only the transport terms required by the diffusion RHS. */
template <int NCH, class Real>
void computeSpmTransport(const SpmConcentrationParams<NCH> &p,
                         const BasicBatchView<const Real> &state,
                         const SpmStateLayout &layout,
                         const BasicStepCtx<Real> &ctx,
                         PerDomain<std::span<Real>>
                           effective_diffusivity,
                         PerDomain<std::span<Real>>
                           molar_flux,
                         SpmTransportCache *cache = nullptr)
{
  const int lanes = state.n_lanes();
  ctx.assert_valid_for(lanes);
  assert(static_cast<int>(effective_diffusivity[0].size()) == lanes
         && static_cast<int>(effective_diffusivity[1].size()) == lanes
         && static_cast<int>(molar_flux[0].size()) == lanes
         && static_cast<int>(molar_flux[1].size()) == lanes);
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane)
      computeSpmTransportLane(p, state, layout, ctx, domain, lane, effective_diffusivity[d][lane], molar_flux[d][lane], cache);
  }
}

/** Reconstruct only surface concentration plus transport for voltage/kinetics paths. */
template <int NCH, class Real>
void computeSpmSurfaceConcentrations(const SpmConcentrationParams<NCH> &p,
                                     const BasicBatchView<const Real> &state,
                                     const SpmStateLayout &layout,
                                     const BasicStepCtx<Real> &ctx,
                                     PerDomain<std::span<Real>>
                                       concentration,
                                     PerDomain<std::span<Real>>
                                       effective_diffusivity,
                                     PerDomain<std::span<Real>>
                                       molar_flux,
                                     SpmTransportCache *cache = nullptr)
{
  const int lanes = state.n_lanes();
  ctx.assert_valid_for(lanes);
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane) {
      Real diffusivity{};
      Real flux{};
      computeSpmTransportLane(p, state, layout, ctx, domain, lane, diffusivity, flux, cache);
      effective_diffusivity[d][lane] = diffusivity;
      molar_flux[d][lane] = flux;
      Real surface{};
      for (int mode = 0; mode < NCH; ++mode)
        surface += p.C[d][0][mode] * state.at(layout.z[d], mode, lane);
      concentration[d][lane] = surface + p.Dout[d][0] * flux / diffusivity;
    }
  }
}

/**
 * Reconstruct Li concentration at surface, interior nodes and centre for every lane/electrode.
 *
 * This is a free-function kernel so the same code is called by the RHS observable stage and
 * by lazy recording. The view is rebound to an integrator's current trial vector before entry.
 * The operation order inside each lane matches legacy Cell_SPM::calcSurfaceConcentration.
 */
template <int NCH, class Real>
void computeSpmConcentrations(const SpmConcentrationParams<NCH> &p,
                              const BasicBatchView<const Real> &state,
                              const SpmStateLayout &layout,
                              const BasicStepCtx<Real> &ctx,
                              PerDomain<std::span<Real>>
                                concentration,
                              PerDomain<std::span<Real>> effective_diffusivity = {},
                              PerDomain<std::span<Real>> molar_flux = {})
{
  assert(domain_value(layout.z, Domain::pos).rows == NCH
         && domain_value(layout.z, Domain::neg).rows == NCH
         && layout.temperature.rows == 1
         && domain_value(layout.diffusion_coefficient, Domain::pos).rows == 1
         && domain_value(layout.diffusion_coefficient, Domain::neg).rows == 1
         && domain_value(layout.specific_surface_area, Domain::pos).rows == 1
         && domain_value(layout.specific_surface_area, Domain::neg).rows == 1
         && domain_value(layout.electrode_thickness, Domain::pos).rows == 1
         && domain_value(layout.electrode_thickness, Domain::neg).rows == 1);
  const int L = state.n_lanes();
  ctx.assert_valid_for(L);
  [[maybe_unused]] constexpr int output_rows = NCH + 2; // surface + NCH interior + centre
  assert(static_cast<int>(concentration[0].size()) == output_rows * L
         && static_cast<int>(concentration[1].size()) == output_rows * L);
  assert((effective_diffusivity[0].empty()
          && effective_diffusivity[1].empty())
         || (static_cast<int>(effective_diffusivity[0].size()) == L
             && static_cast<int>(effective_diffusivity[1].size()) == L));
  assert((molar_flux[0].empty() && molar_flux[1].empty())
         || (static_cast<int>(molar_flux[0].size()) == L
             && static_cast<int>(molar_flux[1].size()) == L));

  const std::span<const Real> T = state.row(layout.temperature.row_begin);
  using std::exp;

  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    const Real D_Td = p.D_T[d];
    const Real sgnd = static_cast<Real>(molar_flux_sign(domain));

    for (int c = 0; c < L; ++c) {
      const Real Arr = (Real{ 1 } / p.T_ref - Real{ 1 } / T[c]) / p.Rg;
      const Real Dt = state.at(layout.diffusion_coefficient[d], 0, c) * exp(D_Td * Arr);
      const Real flux_den = state.at(layout.specific_surface_area[d], 0, c) * p.n * p.F
                            * state.at(layout.electrode_thickness[d], 0, c);
      const Real molarFlux = sgnd * ctx.i_app[c] / flux_den;
      if (!effective_diffusivity[d].empty())
        effective_diffusivity[d][c] = Dt;
      if (!molar_flux[d].empty())
        molar_flux[d][c] = molarFlux;

      for (int node = 0; node < NCH + 1; ++node) {
        Real acc{};
        for (int j = 0; j < NCH; ++j)
          acc += p.C[d][node][j] * state.at(layout.z[d], j, c);
        concentration[d][static_cast<std::size_t>(node) * L + c] = acc + p.Dout[d][node] * molarFlux / Dt;
      }

      Real centre_acc{};
      for (int node = 0; node < NCH + 1; ++node)
        centre_acc += p.Cc[node] * concentration[d][static_cast<std::size_t>(node) * L + c];
      concentration[d][static_cast<std::size_t>(NCH + 1) * L + c] = p.cc_coeff * (centre_acc + molarFlux * p.R[d] / Dt);
    }
  }
}

template <int NCH>
struct SpmElectricalParams
{
  SpmConcentrationParams<NCH> concentration{};
  PerDomain<ElectrodeParams> electrode{};
  PerDomain<IndexedPiecewiseLinear> electrode_ocv{};
  IndexedPiecewiseLinear total_entropic_coefficient{};
  IndexedPiecewiseLinear negative_entropic_coefficient{};

  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t electrolyte_concentration{ 1000.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t sei_resistivity_area{};
};

template <class Real>
struct BasicSpmObservables
{
  PerDomain<std::span<Real>> concentration{}; //!< node-major: [node*n_lanes + lane]
  PerDomain<std::span<Real>> effective_diffusivity{};
  PerDomain<std::span<Real>> molar_flux{};
  PerDomain<std::span<Real>> surface_stoichiometry{};
  PerDomain<std::span<Real>> exchange_current_density{};
  PerDomain<std::span<Real>> overpotential{};
  PerDomain<std::span<Real>> electrode_ocv{};
  std::span<Real> entropic_coefficient{};
  std::span<Real> negative_entropic_coefficient{};
  std::span<Real> open_circuit_voltage{};
  std::span<Real> resistance{};
  std::span<Real> terminal_voltage{};
  std::span<Real> reversible_heat{};
  std::span<Real> reaction_heat{};
  std::span<Real> ohmic_heat{};
  std::span<Real> total_heat{};
};

template <int NCH, class Real = real_t>
class SpmObservableScratch
{
public:
  explicit SpmObservableScratch(int n_lanes)
    : n_lanes_{ n_lanes }, storage_(required_size(n_lanes))
  {
    assert(n_lanes > 0);
  }

  BasicSpmObservables<Real> view()
  {
    return view(n_lanes_);
  }

  BasicSpmObservables<Real> view(int active_lanes)
  {
    assert(active_lanes > 0 && active_lanes <= n_lanes_);
    BasicSpmObservables<Real> output;
    std::size_t cursor = 0;
    auto take = [&](std::size_t count) {
      auto result = std::span<Real>{ storage_ }.subspan(cursor, count);
      cursor += count;
      return result;
    };
    const auto L = static_cast<std::size_t>(active_lanes);
    for (const Domain domain : domains)
      output.concentration[domain_index(domain)] = take(static_cast<std::size_t>(NCH + 2) * L);
    for (auto *field : { &output.effective_diffusivity,
                         &output.molar_flux,
                         &output.surface_stoichiometry,
                         &output.exchange_current_density,
                         &output.overpotential,
                         &output.electrode_ocv })
      for (const Domain domain : domains)
        (*field)[domain_index(domain)] = take(L);
    output.entropic_coefficient = take(L);
    output.negative_entropic_coefficient = take(L);
    output.open_circuit_voltage = take(L);
    output.resistance = take(L);
    output.terminal_voltage = take(L);
    output.reversible_heat = take(L);
    output.reaction_heat = take(L);
    output.ohmic_heat = take(L);
    output.total_heat = take(L);
    assert(cursor <= storage_.size());
    return output;
  }

  int n_lanes() const { return n_lanes_; }

private:
  static std::size_t required_size(int n_lanes)
  {
    return static_cast<std::size_t>(n_lanes) * (2 * (NCH + 2) + 21);
  }

  int n_lanes_{};
  std::vector<Real> storage_{};
};

template <int NCH, class Real>
[[nodiscard]] slide::Status computeSpmObservables(const SpmElectricalParams<NCH> &p,
                                                  const BasicBatchView<const Real> &state,
                                                  const SpmStateLayout &layout,
                                                  const BasicStepCtx<Real> &ctx,
                                                  BasicSpmObservables<Real>
                                                    output,
                                                  SpmTransportCache *cache = nullptr)
{
  const int L = state.n_lanes();
  ctx.assert_valid_for(L);
  computeSpmSurfaceConcentrations(p.concentration,
                                  state,
                                  layout,
                                  ctx,
                                  output.concentration,
                                  output.effective_diffusivity,
                                  output.molar_flux,
                                  cache);

  const std::span<const Real> temperature = state.row(layout.temperature.row_begin);
  using std::asinh;
  using std::exp;
  using std::sqrt;

  for (int lane = 0; lane < L; ++lane) {
    const Real T = temperature[lane];
    const Real arrhenius = (Real{ 1 } / p.reference_temperature - Real{ 1 } / T) / p.Rg;

    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      const auto &electrode = p.electrode[d];
      const Real cs = output.concentration[d][lane];
      const Real z_surface = cs / electrode.cs_max;
      if (!(primal_value(z_surface) > 0.0 && primal_value(z_surface) < 1.0))
        return slide::Status::Invalid_states;

      const Real reaction_rate = electrode.reaction_rate_ref
                                 * exp(electrode.reaction_activation * arrhenius);
      const Real exchange_current = reaction_rate * p.n * p.F
                                    * sqrt(p.electrolyte_concentration * cs
                                           * (electrode.cs_max - cs));
      const Real area = state.at(layout.specific_surface_area[d], 0, lane);
      const Real thickness = state.at(layout.electrode_thickness[d], 0, lane);
      const Real argument = Real{ 0.5 * molar_flux_sign(domain) } * ctx.i_app[lane]
                            / (area * thickness * exchange_current);

      output.surface_stoichiometry[d][lane] = z_surface;
      output.exchange_current_density[d][lane] = exchange_current;
      output.overpotential[d][lane] = Real{ 2 } * p.Rg * T / (p.n * p.F) * asinh(argument);
      output.electrode_ocv[d][lane] = p.electrode_ocv[d].eval(z_surface);
    }

    const auto neg = domain_index(Domain::neg);
    const auto pos = domain_index(Domain::pos);
    const Real z_pos = output.surface_stoichiometry[pos][lane];
    const Real d_ocv = p.total_entropic_coefficient.eval(z_pos);
    const Real d_ocv_neg = p.negative_entropic_coefficient.eval(z_pos);
    const Real ocv = output.electrode_ocv[pos][lane] - output.electrode_ocv[neg][lane]
                     + (T - p.reference_temperature) * d_ocv;

    const Real area_neg = state.at(layout.specific_surface_area[neg], 0, lane)
                          * p.electrode_area
                          * state.at(layout.electrode_thickness[neg], 0, lane);
    const Real area_pos = state.at(layout.specific_surface_area[pos], 0, lane)
                          * p.electrode_area
                          * state.at(layout.electrode_thickness[pos], 0, lane);
    const Real resistance = state.at(layout.sei_thickness, 0, lane) * p.sei_resistivity_area / area_neg
                            + state.at(layout.specific_resistance[neg], 0, lane) / area_neg
                            + state.at(layout.specific_resistance[pos], 0, lane) / area_pos
                            + state.at(layout.current_collector_resistance, 0, lane) / p.electrode_area;
    const Real current = ctx.i_app[lane] * p.electrode_area;
    const Real reaction_heat = current
                               * (output.overpotential[neg][lane] - output.overpotential[pos][lane]);
    const Real reversible_heat = -current * T * d_ocv;
    const Real ohmic_heat = current * current * resistance;

    output.entropic_coefficient[lane] = d_ocv;
    output.negative_entropic_coefficient[lane] = d_ocv_neg;
    output.open_circuit_voltage[lane] = ocv;
    output.resistance[lane] = resistance;
    output.terminal_voltage[lane] = ocv + output.overpotential[pos][lane]
                                    - output.overpotential[neg][lane] - resistance * current;
    output.reversible_heat[lane] = reversible_heat;
    output.reaction_heat[lane] = reaction_heat;
    output.ohmic_heat[lane] = ohmic_heat;
    output.total_heat[lane] = reversible_heat + reaction_heat + ohmic_heat;
  }
  return slide::Status::Success;
}

} // namespace slide::core
