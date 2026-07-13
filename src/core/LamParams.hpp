/**
 * @file LamParams.hpp
 * @brief Cold, batch-shared loss-of-active-material constants and their validity gate.
 *
 * Owns: `LamParams` (models 1--4), its model-bit vocabulary, and `validateLamParams`.
 * Implements PLAN.md §3.11. Cold: built once per batch. The compiled positive OCV curve
 * is stored by value because LAM-3 needs it inside the lane sweep.
 * Split from `Lam.hpp` for MC-5 (see `SeiParams.hpp`).
 * @surface support
 */

#pragma once

#include "../types/Status.hpp"
#include "../types/Status.hpp"
#include "AgeingModelMask.hpp"
#include "CellDesign.hpp"
#include "CompiledCurve.hpp"
#include "Numeric.hpp"

#include <array>
#include <cstdint>

namespace slide::core {

constexpr std::uint8_t lam_model_bit(unsigned model)
{
  return ageing_model_bit<4>(model);
}

struct LamParams
{
  std::uint8_t model_mask{};
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t reference_temperature{ 298.15 };
  PerDomain<real_t> particle_radius{};

  PerDomain<real_t> model1_stress_coefficient{};
  PerDomain<real_t> model2_linear_flux{};
  PerDomain<real_t> model2_sqrt_flux{};
  real_t model2_activation{};
  real_t model3_k{};
  real_t model3_k_activation{};
  real_t model3_equilibrium_potential{ 4.1 };
  PerDomain<real_t> model4_area_coefficient{};
  IndexedPiecewiseLinear positive_ocv{};
};

[[nodiscard]] inline slide::Status validateLamParams(const LamParams &p)
{
  if (!valid_ageing_model_mask<4>(p.model_mask) || !p.positive_ocv.valid())
    return slide::Status::Invalid_parameters;
  const std::array scalars{ p.F,
                            p.Rg,
                            p.n,
                            p.reference_temperature,
                            p.model2_activation,
                            p.model3_k,
                            p.model3_k_activation,
                            p.model3_equilibrium_potential };
  for (const real_t value : scalars)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    const std::array values{ p.particle_radius[d],
                             p.model1_stress_coefficient[d],
                             p.model2_linear_flux[d],
                             p.model2_sqrt_flux[d],
                             p.model4_area_coefficient[d] };
    for (const real_t value : values)
      if (!is_finite(value))
        return slide::Status::Invalid_parameters;
    if (!(p.particle_radius[d] > 0.0))
      return slide::Status::Invalid_parameters;
  }
  return (p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0
          && p.reference_temperature > 0.0 && p.model3_k >= 0.0)
           ? slide::Status::Success
           : slide::Status::Invalid_parameters;
}

} // namespace slide::core
