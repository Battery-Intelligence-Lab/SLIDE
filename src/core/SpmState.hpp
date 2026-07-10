/**
 * @file SpmState.hpp
 * @brief Named arena layout for the v4 SPM composition.
 */

#pragma once

#include "BatchBuilder.hpp"
#include "CellDesign.hpp"

namespace slide::core {

struct SpmStateLayout
{
  PerDomain<StateSlice> z{};
  StateSlice temperature{};
  StateSlice sei_thickness{};
  StateSlice lost_lithium{};
  StateSlice crack_surface{};
  StateSlice plated_lithium_thickness{};
  PerDomain<StateSlice> active_fraction{};
  PerDomain<StateSlice> diffusion_coefficient{};
  PerDomain<StateSlice> electrode_thickness{};
  PerDomain<StateSlice> specific_surface_area{};
  PerDomain<StateSlice> specific_resistance{};
  StateSlice current_collector_resistance{};
};

template <int NCH>
SpmStateLayout declareSpmState(BatchBuilder &builder,
                               bool temperature_is_ode = true,
                               bool ageing_is_ode = true)
{
  const auto temperature_role = temperature_is_ode ? StateRole::ode : StateRole::algebraic;
  const auto ageing_role = ageing_is_ode ? StateRole::ode : StateRole::algebraic;
  SpmStateLayout layout;
  domain_value(layout.z, Domain::neg) = builder.declare({ "zn", NCH, Unit::none });
  domain_value(layout.z, Domain::pos) = builder.declare({ "zp", NCH, Unit::none });
  layout.temperature = builder.declare({ "T", 1, Unit::K, temperature_role });
  layout.sei_thickness = builder.declare({ "delta_sei", 1, Unit::m, ageing_role });
  layout.lost_lithium = builder.declare({ "lost_lithium", 1, Unit::C, ageing_role });
  layout.crack_surface = builder.declare({ "crack_surface", 1, Unit::m2, ageing_role });
  layout.plated_lithium_thickness = builder.declare({ "plated_lithium_thickness", 1, Unit::m, ageing_role });
  domain_value(layout.active_fraction, Domain::neg) = builder.declare({ "active_fraction_neg", 1, Unit::none, ageing_role });
  domain_value(layout.active_fraction, Domain::pos) = builder.declare({ "active_fraction_pos", 1, Unit::none, ageing_role });
  domain_value(layout.diffusion_coefficient, Domain::neg) = builder.declare({ "diffusion_neg", 1, Unit::m2_s, ageing_role });
  domain_value(layout.diffusion_coefficient, Domain::pos) = builder.declare({ "diffusion_pos", 1, Unit::m2_s, ageing_role });
  domain_value(layout.electrode_thickness, Domain::neg) = builder.declare({ "thickness_neg", 1, Unit::m, ageing_role });
  domain_value(layout.electrode_thickness, Domain::pos) = builder.declare({ "thickness_pos", 1, Unit::m, ageing_role });
  domain_value(layout.specific_surface_area, Domain::neg) = builder.declare({ "specific_area_neg", 1, Unit::inv_m, ageing_role });
  domain_value(layout.specific_surface_area, Domain::pos) = builder.declare({ "specific_area_pos", 1, Unit::inv_m, ageing_role });
  domain_value(layout.specific_resistance, Domain::neg) = builder.declare({ "specific_resistance_neg", 1, Unit::ohm_m2, StateRole::algebraic });
  domain_value(layout.specific_resistance, Domain::pos) = builder.declare({ "specific_resistance_pos", 1, Unit::ohm_m2, StateRole::algebraic });
  layout.current_collector_resistance = builder.declare({ "current_collector_resistance", 1, Unit::ohm_m2, StateRole::algebraic });
  return layout;
}

} // namespace slide::core
