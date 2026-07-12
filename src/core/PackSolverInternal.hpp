/**
 * @file PackSolverInternal.hpp
 * @brief Internal numeric helpers shared only with focused solver tests.
 * @details Owns M0.3 diagnostic-bound saturation; hot, allocation-free, and
 *          deliberately excluded from the public PackSolver surface.
 */

#pragma once

#include "Numeric.hpp"

namespace slide::core::detail {

[[nodiscard]] real_t conservativePackRoundoffBound(real_t accumulation_ratio,
                                                   real_t current_scale,
                                                   real_t operation_scale) noexcept;

} // namespace slide::core::detail
