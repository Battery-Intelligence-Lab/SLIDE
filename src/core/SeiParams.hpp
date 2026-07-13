/**
 * @file SeiParams.hpp
 * @brief Cold, batch-shared SEI constants and their validity gate.
 *
 * Owns: `SeiParams` (models 1--4), its model-bit vocabulary, and `validateSeiParams`.
 * Implements PLAN.md §3.11 (cold parameter compilation). Cold: built once per batch.
 * Split from `Sei.hpp` for MC-5 — `SpmFactoryInput` stores these by value, so the type
 * must reach the public factory without dragging the SEI kernel in with it.
 * @surface support
 */

#pragma once

#include "AgeingModelMask.hpp"
#include "Numeric.hpp"

#include <array>
#include <cstdint>

namespace slide::core {

constexpr std::uint8_t sei_model_bit(unsigned model)
{
  return ageing_model_bit<4>(model);
}

/** Batch-shared, cold-built constants for legacy-compatible SEI models 1--4. */
struct SeiParams
{
  std::uint8_t model_mask{};
  bool reduce_active_fraction{};

  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t n_sei{ 1.0 };
  real_t alpha_sei{ 1.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t negative_particle_radius{};
  real_t sei_resistivity_area{};
  real_t sei_equilibrium_potential{ 0.4 };
  real_t sei_molar_volume{ 100e3 };
  real_t electrolyte_reactant_concentration{ 4.541e-3 };
  real_t main_molar_volume{ 13.0 };
  real_t side_molar_volume{ 64.39 };
  real_t porosity_coefficient{ 3.0 * 7.5e-7 };

  real_t model1_k{ 0.075e-14 };
  real_t model1_k_activation{ 130e3 };
  real_t model2_k{ 2.75e-11 };
  real_t model2_k_activation{ 130e3 };
  real_t model2_D{ 1.125e-14 };
  real_t model2_D_activation{ 20e3 };
  real_t model3_k{ 1.1458e-15 };
  real_t model3_k_activation{ 65e3 };
  real_t model3_D{ 0.25e-15 };
  real_t model3_D_activation{ 200e3 };
  real_t model4_k{ 3.75e-15 };
  real_t model4_k_activation{ 130000.0 / 1.5 };
  real_t model4_D{ 0.5e-16 / 15.0 };
  real_t model4_D_activation{ 80000.0 };
};

[[nodiscard]] inline slide::Status validateSeiParams(const SeiParams &p)
{
  if (!valid_ageing_model_mask<4>(p.model_mask))
    return slide::Status::Invalid_parameters;
  const std::array values{ p.F,
                           p.Rg,
                           p.n,
                           p.n_sei,
                           p.alpha_sei,
                           p.reference_temperature,
                           p.electrode_area,
                           p.negative_particle_radius,
                           p.sei_resistivity_area,
                           p.sei_equilibrium_potential,
                           p.sei_molar_volume,
                           p.electrolyte_reactant_concentration,
                           p.main_molar_volume,
                           p.side_molar_volume,
                           p.porosity_coefficient,
                           p.model1_k,
                           p.model1_k_activation,
                           p.model2_k,
                           p.model2_k_activation,
                           p.model2_D,
                           p.model2_D_activation,
                           p.model3_k,
                           p.model3_k_activation,
                           p.model3_D,
                           p.model3_D_activation,
                           p.model4_k,
                           p.model4_k_activation,
                           p.model4_D,
                           p.model4_D_activation };
  for (const real_t value : values)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  if (!(p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0 && p.n_sei > 0.0
        && p.reference_temperature > 0.0 && p.electrode_area > 0.0
        && p.negative_particle_radius > 0.0 && p.sei_molar_volume > 0.0
        && p.model1_k >= 0.0 && p.model2_k >= 0.0 && p.model2_D > 0.0
        && p.model3_k >= 0.0 && p.model3_D > 0.0 && p.model4_k >= 0.0
        && p.model4_D > 0.0))
    return slide::Status::Invalid_parameters;
  return slide::Status::Success;
}

} // namespace slide::core
