/**
 * @file LithiumPlatingParams.hpp
 * @brief Cold, batch-shared lithium-plating constants, derived scales, and validity gate.
 *
 * Owns: `LithiumPlatingParams`, the `LithiumPlatingScales` derived from them,
 * `tryLithiumPlatingScales`, and `validateLithiumPlatingParams`. Implements PLAN.md §3.11.
 * Cold: built once per batch. Split from `LithiumPlating.hpp` for MC-5 (see `SeiParams.hpp`).
 * @surface support
 */

#pragma once

#include "../types/Status.hpp"
#include "Numeric.hpp"

#include <array>

namespace slide::core {

struct LithiumPlatingParams
{
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t n_plating{ 1.0 };
  real_t alpha_plating{ 1.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t sei_resistivity_area{};
  real_t equilibrium_potential{};
  real_t plated_lithium_molar_density{ 10e6 };
  real_t reaction_rate_ref{ 4.5e-10 };
  real_t reaction_rate_activation{ -2.014008e5 };
};

struct LithiumPlatingScales
{
  real_t faradaic_scale{};        //!< n_plating * F
  real_t kinetic_charge{};        //!< n * F
  real_t thickness_denominator{}; //!< n_plating * F * plated molar density
};

[[nodiscard]] inline bool tryLithiumPlatingScales(
  const LithiumPlatingParams &p,
  LithiumPlatingScales &scales) noexcept
{
  LithiumPlatingScales candidate;
  if (!try_multiply_nonnegative(
        p.n_plating, p.F, candidate.faradaic_scale)
      || !try_multiply_nonnegative(
        p.n, p.F, candidate.kinetic_charge)
      || !try_multiply_nonnegative(
        candidate.faradaic_scale,
        p.plated_lithium_molar_density,
        candidate.thickness_denominator))
    return false;
  scales = candidate;
  return true;
}

[[nodiscard]] inline slide::Status validateLithiumPlatingParams(
  const LithiumPlatingParams &p)
{
  const std::array values{ p.F,
                           p.Rg,
                           p.n,
                           p.n_plating,
                           p.alpha_plating,
                           p.reference_temperature,
                           p.electrode_area,
                           p.sei_resistivity_area,
                           p.equilibrium_potential,
                           p.plated_lithium_molar_density,
                           p.reaction_rate_ref,
                           p.reaction_rate_activation };
  for (const real_t value : values)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  LithiumPlatingScales scales;
  if (!tryLithiumPlatingScales(p, scales))
    return slide::Status::Invalid_parameters;
  return (p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0 && p.n_plating > 0.0
          && p.reference_temperature > 0.0 && p.electrode_area > 0.0
          && p.plated_lithium_molar_density > 0.0 && p.reaction_rate_ref >= 0.0)
           ? slide::Status::Success
           : slide::Status::Invalid_parameters;
}

} // namespace slide::core
