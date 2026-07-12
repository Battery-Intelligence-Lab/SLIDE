/**
 * @file BpxExpression.hpp
 * @brief Opaque cold-path seam for compiled BPX scalar expressions.
 *
 * M0.7 / 9C-3 cold-path contract: this header exposes compilation and
 * evaluation only. The AST, scanner, parser, and evaluator layout remain
 * translation-unit-local so BPX grammar internals cannot become public API.
 */

#pragma once

#include "../CellDesign.hpp"

#include <span>
#include <string>
#include <string_view>

namespace slide::core::detail {

bool evaluateBpxExpressionSamples(std::string_view source,
                                  std::span<const real_t>
                                    samples,
                                  std::span<real_t>
                                    output,
                                  std::string &diagnostic);
bool sampleBpxExpressionCurve(std::string_view source,
                              OCVCurve &output,
                              std::string &diagnostic);

} // namespace slide::core::detail
