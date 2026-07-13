/**
 * @file ForwardSensitivity.hpp
 * @brief True dual-number forward sensitivities for the isothermal SPM.
 * @surface api
 */

#pragma once

#include "Experiment.hpp"

#include <array>
#include <cassert>
#include <span>
#include <string_view>
#include <vector>

namespace slide::core {

enum class SensitivityParameter : unsigned char {
  negative_diffusivity,
  positive_diffusivity,
  negative_reaction_rate,
  positive_reaction_rate,
  contact_resistance,
  nominal_capacity,
  negative_minimum_stoichiometry,
  negative_maximum_stoichiometry,
  positive_minimum_stoichiometry,
  positive_maximum_stoichiometry
};

inline constexpr std::array supported_sensitivity_parameters{
  SensitivityParameter::negative_diffusivity,
  SensitivityParameter::positive_diffusivity,
  SensitivityParameter::negative_reaction_rate,
  SensitivityParameter::positive_reaction_rate,
  SensitivityParameter::contact_resistance,
  SensitivityParameter::nominal_capacity,
  SensitivityParameter::negative_minimum_stoichiometry,
  SensitivityParameter::negative_maximum_stoichiometry,
  SensitivityParameter::positive_minimum_stoichiometry,
  SensitivityParameter::positive_maximum_stoichiometry,
};

std::string_view sensitivityParameterName(SensitivityParameter parameter);
[[nodiscard]] slide::Status parseSensitivityParameter(
  std::string_view name,
  SensitivityParameter &parameter);

struct ForwardSensitivitySolution
{
  std::vector<real_t> time{};
  std::vector<real_t> terminal_voltage{};
  std::vector<SensitivityParameter> parameters{};
  /** sample-major: derivative[sample*n_parameters + parameter]. */
  std::vector<real_t> derivative{};

  std::span<const real_t> derivativeAt(std::size_t sample) const
  {
    assert(sample < time.size());
    return std::span<const real_t>{ derivative }.subspan(
      sample * parameters.size(), parameters.size());
  }
};

/**
 * Solve one fixed-duration CC segment. A C-rate control differentiates current
 * through nominal capacity; an ampere control keeps imposed current fixed.
 */
[[nodiscard]] slide::Status solveCcForwardSensitivities(
  const SpmFactoryInput &input,
  int nch,
  real_t control_magnitude,
  bool control_is_c_rate,
  Direction direction,
  real_t duration,
  real_t sample_step,
  std::span<const SensitivityParameter> parameters,
  ForwardSensitivitySolution &output);

real_t sensitivityParameterValue(const SpmFactoryInput &input,
                                 SensitivityParameter parameter);

} // namespace slide::core
