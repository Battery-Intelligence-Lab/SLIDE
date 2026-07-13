/**
 * @file AgeingModelMask.hpp
 * @brief The model-selection bit vocabulary shared by every ageing mechanism.
 *
 * Owns: `ageing_model_bit`, `ageing_model_mask`, and the two mask validity predicates.
 * Implements PLAN.md §3.11 (which models a batch runs is cold configuration).
 * Cold: read once per batch build, and by callers assembling a `*Params` mask.
 * Split from `AgeingKernel.hpp` for MC-5: the mask vocabulary is part of the public
 * parameter surface, the lane-sweep scaffolding around it is not.
 * @surface support
 */

#pragma once

#include <cstdint>

namespace slide::core {

template <unsigned ModelCount>
constexpr std::uint8_t ageing_model_bit(unsigned model) noexcept
{
  static_assert(ModelCount > 0 && ModelCount <= 8);
  return model >= 1 && model <= ModelCount
           ? static_cast<std::uint8_t>(std::uint16_t{ 1 } << (model - 1))
           : std::uint8_t{};
}

template <unsigned ModelCount>
constexpr std::uint8_t ageing_model_mask() noexcept
{
  static_assert(ModelCount > 0 && ModelCount <= 8);
  return static_cast<std::uint8_t>((std::uint16_t{ 1 } << ModelCount) - 1);
}

template <unsigned ModelCount>
constexpr bool valid_ageing_model_mask(std::uint8_t mask) noexcept
{
  return mask != 0
         && (mask & static_cast<std::uint8_t>(~ageing_model_mask<ModelCount>())) == 0;
}

template <unsigned ModelCount>
constexpr bool valid_optional_ageing_model_mask(std::uint8_t mask) noexcept
{
  return (mask & static_cast<std::uint8_t>(~ageing_model_mask<ModelCount>())) == 0;
}

} // namespace slide::core
