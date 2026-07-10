/**
 * @file ParameterSet.hpp
 * @brief PyBaMM/BPX-named cold parameter absorption and traceability.
 */

#pragma once

#include "SpmFactory.hpp"

#include <filesystem>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace slide::core {

using ParameterValue = std::variant<real_t, OCVCurve>;

struct ParameterDescription
{
  std::string name{};
  ParameterValue value{};
  std::string provenance{};
};

/** Cold, value-semantic parameter bag keyed by PyBaMM's public SI names. */
class ParameterSet
{
public:
  [[nodiscard]] slide::Status set(std::string name, ParameterValue value,
                                  std::string provenance = "user");
  [[nodiscard]] slide::Status update(
    std::span<const ParameterDescription> values);
  bool contains(std::string_view name) const;
  const ParameterValue *find(std::string_view name) const;
  const real_t *findScalar(std::string_view name) const;
  const OCVCurve *findCurve(std::string_view name) const;
  std::vector<ParameterDescription> describe() const;
  std::size_t size() const { return values_.size(); }

  [[nodiscard]] slide::Status toSpmInput(SpmFactoryInput &output) const;
  [[nodiscard]] static slide::Status chen2020(ParameterSet &output);
  [[nodiscard]] static slide::Status fromBpxJson(
    std::string_view json,
    ParameterSet &output,
    std::string &diagnostic);
  [[nodiscard]] static slide::Status fromBpxFile(
    const std::filesystem::path &path,
    ParameterSet &output,
    std::string &diagnostic);

  static std::string canonicalName(std::string_view name);

private:
  struct Entry
  {
    ParameterValue value{};
    std::string provenance{};
  };
  std::map<std::string, Entry, std::less<>> values_{};
};

} // namespace slide::core
