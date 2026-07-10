/**
 * @file ThermalLumped.hpp
 * @brief Allocation-free lumped thermal RHS for the v4 core (PLAN.md Phase 1).
 */

#pragma once

#include "BatchBuilder.hpp"
#include "CellDesign.hpp"
#include "SpmObservables.hpp"

#include <cassert>
#include <cmath>

namespace slide::core {

/**
 * Cold-compiled thermal constants. `thermal_capacity` is rho*Cp*volume [J/K] and
 * `environment_conductance` is h*A [W/K], so the hot kernel needs only one multiply
 * and one division per lane beyond summing its heat inputs.
 */
struct ThermalLumpedParams
{
  real_t thermal_capacity{};
  real_t environment_conductance{};
  real_t environment_temperature{};
};

/** Compile a physical thermal description into hot SI scalars, atomically on failure. */
[[nodiscard]] inline slide::Status compileThermalLumped(const ThermalDesign &design,
                                                        ThermalLumpedParams &output)
{
  output = {};
  if (!(is_finite(design.density) && design.density > 0.0
        && is_finite(design.heat_capacity) && design.heat_capacity > 0.0
        && is_finite(design.volume) && design.volume > 0.0
        && is_finite(design.surface_area) && design.surface_area >= 0.0
        && is_finite(design.h_conv) && design.h_conv >= 0.0
        && is_finite(design.environment_temperature)
        && design.environment_temperature > 0.0))
    return slide::Status::Invalid_parameters;

  const real_t capacity = design.density * design.heat_capacity * design.volume;
  const real_t conductance = design.h_conv * design.surface_area;
  if (!(is_finite(capacity) && capacity > 0.0 && is_finite(conductance)))
    return slide::Status::Invalid_parameters;

  output = { .thermal_capacity = capacity,
             .environment_conductance = conductance,
             .environment_temperature = design.environment_temperature };
  return slide::Status::Success;
}

/**
 * State owned by the thermal component. The generated-energy and elapsed-time rows replace
 * legacy `Therm_Qgen`/`Therm_time` object members, making checkpoints self-contained (P1-G4).
 * They are monotone lifetime integrals in v4; no resettable hidden accumulator is needed because
 * temperature itself is integrated from the same instantaneous RHS.
 */
struct ThermalLumpedLayout
{
  StateSlice temperature{};
  StateSlice external_heat_flow{};    //!< q_ext, positive into the cell [W]
  StateSlice generated_heat_energy{}; //!< integral of internal heat generation [J]
  StateSlice thermal_elapsed_time{};  //!< time integrated by the thermal component [s]
};

inline ThermalLumpedLayout declareThermalLumped(BatchBuilder &builder,
                                                StateSlice temperature)
{
  assert(temperature.rows == 1);
  return { .temperature = temperature,
           .external_heat_flow = builder.reserve_thermal_flux(),
           .generated_heat_energy = builder.declare({ "generated_heat_energy", 1, Unit::J }),
           .thermal_elapsed_time = builder.declare({ "thermal_elapsed_time", 1, Unit::s }) };
}

/**
 * Add the lumped thermal contribution to an already-zeroed derivative arena.
 *
 *   C_th dT/dt = Q_internal + q_ext + hA(T_env - T)
 *
 * `q_ext` is the reserved cross-batch seam; pack thermal adjacency writes it between batch
 * steps. Internal heat comes from the one shared SPM observable stage, so recording and RHS
 * physics cannot drift apart.
 */
template <class Real>
[[nodiscard]] slide::Status addThermalLumpedRhs(const ThermalLumpedParams &p,
                                                const BasicBatchView<const Real> &state,
                                                BasicBatchView<Real>
                                                  derivative,
                                                const ThermalLumpedLayout &layout,
                                                const BasicSpmObservables<Real> &observables,
                                                const BasicStepCtx<Real> &ctx)
{
  const int lanes = state.n_lanes();
  assert(derivative.n_lanes() == lanes && derivative.n_rows() == state.n_rows());
  assert(layout.temperature.rows == 1 && layout.external_heat_flow.rows == 1
         && layout.generated_heat_energy.rows == 1 && layout.thermal_elapsed_time.rows == 1);
  assert(static_cast<int>(observables.total_heat.size()) == lanes);
  ctx.assert_valid_for(lanes);

  for (int lane = 0; lane < lanes; ++lane) {
    const Real temperature = state.at(layout.temperature, 0, lane);
    const Real internal_heat = observables.total_heat[static_cast<std::size_t>(lane)];
    const Real external_heat = state.at(layout.external_heat_flow, 0, lane);
    if (!(is_finite_primal(temperature) && primal_value(temperature) > 0.0
          && is_finite_primal(internal_heat) && is_finite_primal(external_heat)))
      return slide::Status::Invalid_states;

    const Real environment_heat = p.environment_conductance
                                  * (p.environment_temperature - temperature);
    derivative.at(layout.temperature, 0, lane) += (internal_heat + external_heat + environment_heat) / p.thermal_capacity;
    derivative.at(layout.generated_heat_energy, 0, lane) += internal_heat;
    derivative.at(layout.thermal_elapsed_time, 0, lane) += Real{ 1 };
  }
  return slide::Status::Success;
}

} // namespace slide::core
