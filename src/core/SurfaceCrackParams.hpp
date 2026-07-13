/**
 * @file SurfaceCrackParams.hpp
 * @brief Cold, batch-shared surface-crack constants and their validity gate.
 *
 * Owns: `SurfaceCrackParams` (models 1--5), its model-bit vocabulary, and
 * `validateSurfaceCrackParams`. Implements PLAN.md §3.11. Cold: built once per batch.
 * Split from `SurfaceCrack.hpp` for MC-5 (see `SeiParams.hpp`).
 * @surface support
 */

#pragma once

#include "AgeingModelMask.hpp"
#include "Numeric.hpp"

#include <array>
#include <cstdint>

namespace slide::core {

constexpr std::uint8_t surface_crack_model_bit(unsigned model)
{
  return ageing_model_bit<5>(model);
}

struct SurfaceCrackParams
{
  std::uint8_t model_mask{};
  bool reduce_negative_diffusivity{};
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n_sei{ 1.0 };
  real_t alpha_sei{ 1.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t negative_cs_max{ 30555.0 };
  real_t sei_resistivity_area{};
  real_t sei_equilibrium_potential{ 0.4 };
  real_t model1_alpha{ 4.25e-5 };
  real_t model2_alpha{ 6.3e-7 };
  real_t model3_alpha{ 2.31e-16 };
  real_t model4_alpha{ 4.3306e-8 };
  real_t model4_max_surface{};
  real_t model5_k{ 1e-18 };
  real_t model5_k_activation{ -127040.0 };
  real_t diffusion_exponent{ 2.0 };
};

[[nodiscard]] inline slide::Status validateSurfaceCrackParams(const SurfaceCrackParams &p)
{
  if (!valid_ageing_model_mask<5>(p.model_mask))
    return slide::Status::Invalid_parameters;
  const std::array values{ p.F,
                           p.Rg,
                           p.n_sei,
                           p.alpha_sei,
                           p.reference_temperature,
                           p.electrode_area,
                           p.negative_cs_max,
                           p.sei_resistivity_area,
                           p.sei_equilibrium_potential,
                           p.model1_alpha,
                           p.model2_alpha,
                           p.model3_alpha,
                           p.model4_alpha,
                           p.model4_max_surface,
                           p.model5_k,
                           p.model5_k_activation,
                           p.diffusion_exponent };
  for (const real_t value : values)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  if (!(p.F > 0.0 && p.Rg > 0.0 && p.n_sei > 0.0
        && p.reference_temperature > 0.0 && p.electrode_area > 0.0
        && p.negative_cs_max > 0.0 && p.model4_max_surface > 0.0
        && p.model5_k >= 0.0 && p.diffusion_exponent > 0.0))
    return slide::Status::Invalid_parameters;
  return slide::Status::Success;
}

} // namespace slide::core
