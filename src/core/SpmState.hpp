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
  PerDomain<StateSlice> active_fraction{};
  PerDomain<StateSlice> diffusion_coefficient{};
  PerDomain<StateSlice> electrode_thickness{};
  PerDomain<StateSlice> specific_surface_area{};
  PerDomain<StateSlice> specific_resistance{};
  StateSlice current_collector_resistance{};
};

template <int NCH>
SpmStateLayout declareSpmState(BatchBuilder &builder)
{
  SpmStateLayout layout;
  domain_value(layout.z, Domain::neg) = builder.declare({ "zn", NCH, Unit::none });
  domain_value(layout.z, Domain::pos) = builder.declare({ "zp", NCH, Unit::none });
  layout.temperature = builder.declare({ "T", 1, Unit::K });
  layout.sei_thickness = builder.declare({ "delta_sei", 1, Unit::m });
  layout.lost_lithium = builder.declare({ "lost_lithium", 1, Unit::C });
  domain_value(layout.active_fraction, Domain::neg) = builder.declare({ "active_fraction_neg", 1, Unit::none });
  domain_value(layout.active_fraction, Domain::pos) = builder.declare({ "active_fraction_pos", 1, Unit::none });
  domain_value(layout.diffusion_coefficient, Domain::neg) = builder.declare({ "diffusion_neg", 1, Unit::m2_s });
  domain_value(layout.diffusion_coefficient, Domain::pos) = builder.declare({ "diffusion_pos", 1, Unit::m2_s });
  domain_value(layout.electrode_thickness, Domain::neg) = builder.declare({ "thickness_neg", 1, Unit::m });
  domain_value(layout.electrode_thickness, Domain::pos) = builder.declare({ "thickness_pos", 1, Unit::m });
  domain_value(layout.specific_surface_area, Domain::neg) = builder.declare({ "specific_area_neg", 1, Unit::inv_m });
  domain_value(layout.specific_surface_area, Domain::pos) = builder.declare({ "specific_area_pos", 1, Unit::inv_m });
  domain_value(layout.specific_resistance, Domain::neg) = builder.declare({ "specific_resistance_neg", 1, Unit::ohm_m2 });
  domain_value(layout.specific_resistance, Domain::pos) = builder.declare({ "specific_resistance_pos", 1, Unit::ohm_m2 });
  layout.current_collector_resistance = builder.declare({ "current_collector_resistance", 1, Unit::ohm_m2 });
  return layout;
}

} // namespace slide::core
