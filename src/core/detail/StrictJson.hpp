/**
 * @file StrictJson.hpp
 * @brief Value-semantic DOM seam for the strict cold JSON parser.
 *
 * M0.7 / 9C-3 cold-path contract: this header owns only the internal DOM and
 * one bounded parse entry point. Grammar, UTF-8 scanning, and resource limits
 * remain translation-unit-local and independent of all wire-format consumers.
 */

#pragma once

#include "../StateArena.hpp"

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace slide::core::detail {

struct StrictJsonValue
{
  enum class Kind : unsigned char {
    null_value,
    boolean,
    number,
    string,
    array,
    object,
  } kind{ Kind::null_value };
  bool boolean{};
  real_t number{};
  std::string string{};
  std::vector<StrictJsonValue> array{};
  std::map<std::string, StrictJsonValue, std::less<>> object{};
};

bool parseStrictJson(std::string_view source,
                     StrictJsonValue &output,
                     std::string &diagnostic);

} // namespace slide::core::detail
