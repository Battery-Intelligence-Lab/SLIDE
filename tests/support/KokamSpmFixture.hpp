/**
 * @file KokamSpmFixture.hpp
 * @brief Exact legacy-Kokam inputs shared by parity and payoff verification.
 */

#pragma once

#include "../../src/slide.hpp"
#include "../../src/core/SpmFactory.hpp"

#include <vector>

namespace slide::test_support {

inline constexpr double kokam_xp_0 = 0.983999588653496;
inline constexpr double kokam_xp_100 = 0.400145394039564;
inline constexpr double kokam_xn_0 = 0.029397569380507;
inline constexpr double kokam_xn_100 = 0.932469496648387;

inline core::OCVCurve copy_curve(const XYdata_ss &curve)
{
  return { .stoichiometry = std::vector<double>(curve.x.begin(), curve.x.end()),
           .value = std::vector<double>(curve.y.begin(), curve.y.end()) };
}

inline core::SpmFactoryInput make_legacy_kokam_input(double initial_soc,
                                                     double initial_temperature,
                                                     double reference_temperature)
{
  auto curves = OCVcurves::makeOCVcurves(cellType::KokamNMC);
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 16.0;
  input.design.electrode_area = 0.1 * 0.2 * 31.0;
  input.design.electrolyte.concentration = PhyConst::C_elec;
  input.design.thermal.reference_temperature = reference_temperature;
  input.design.thermal.environment_temperature = settings::T_ENV;
  input.initial_soc = initial_soc;
  input.initial_temperature = initial_temperature;
  input.initial_sei_thickness = 1e-9;
  input.initial_lost_lithium = 0.0;
  input.initial_crack_surface_fraction = 0.01;
  input.initial_plated_lithium_thickness = 0.0;
  input.initial_current_collector_resistance = 0.2325e-3;
  input.sei_resistivity_area = 2037.4 * 50.0;
  input.total_entropic_coefficient = copy_curve(curves.dOCV_tot);
  input.negative_entropic_coefficient = copy_curve(curves.dOCV_neg);

  auto &negative = core::domain_value(input.design.electrode, core::Domain::neg);
  negative.thickness = 74.883947e-6;
  negative.porosity = 0.0;
  negative.active_fraction = 0.5;
  negative.particle_radius = 1.25e-5;
  negative.active_material.ocv = copy_curve(curves.OCV_neg);
  negative.active_material.cs_max = 30555.0;
  negative.active_material.x_0 = kokam_xn_0;
  negative.active_material.x_100 = kokam_xn_100;
  negative.active_material.D_s = { .reference_value = 7e-14,
                                   .activation_energy = 35000.0 / 5.0,
                                   .reference_temperature = reference_temperature };
  negative.active_material.k_ct = { .reference_value = 1.764e-11,
                                    .activation_energy = 20000.0,
                                    .reference_temperature = reference_temperature };

  auto &positive = core::domain_value(input.design.electrode, core::Domain::pos);
  positive.thickness = 86.87357e-6;
  positive.porosity = 0.0;
  positive.active_fraction = 0.5;
  positive.particle_radius = 8.5e-6;
  positive.active_material.ocv = copy_curve(curves.OCV_pos);
  positive.active_material.cs_max = 51385.0;
  positive.active_material.x_0 = kokam_xp_0;
  positive.active_material.x_100 = kokam_xp_100;
  positive.active_material.D_s = { .reference_value = 8e-14,
                                   .activation_energy = 29000.0,
                                   .reference_temperature = reference_temperature };
  positive.active_material.k_ct = { .reference_value = 5e-11,
                                    .activation_energy = 58000.0,
                                    .reference_temperature = reference_temperature };

  for (const auto domain : core::domains)
    input.initial_specific_resistance[core::domain_index(domain)] = 2.8e-3;
  return input;
}

inline Status initialize_legacy_kokam(Cell_SPM &cell, double initial_soc,
                                      double current_A)
{
  cell.setBlockDegAndTherm(true);
  cell.setC({ kokam_xp_0 + initial_soc * (kokam_xp_100 - kokam_xp_0),
              kokam_xn_0 + initial_soc * (kokam_xn_100 - kokam_xn_0) });
  auto status = cell.setSOC(initial_soc, false, false);
  if (status != Status::Success)
    return status;
  return cell.setCurrent(current_A, false, false);
}

} // namespace slide::test_support
