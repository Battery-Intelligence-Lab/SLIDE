/**
 * @file BpxParameterReader.cpp
 * @brief BPX 1.x wire absorption and bounded file loading.
 */

#include "ParameterSet.hpp"
#include "BoundedFileReader.hpp"
#include "detail/BpxExpression.hpp"
#include "detail/ParameterCurve.hpp"
#include "detail/StrictJson.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace slide::core {
namespace {

  constexpr std::size_t max_bpx_json_bytes = 4U * 1024U * 1024U;
  using JsonValue = detail::StrictJsonValue;

  void assignDiagnosticNoThrow(std::string &target,
                               std::string_view message) noexcept
  {
    try {
      target.assign(message);
    } catch (...) {
      target.clear();
    }
  }

  slide::Status allocationFailure(std::string &diagnostic,
                                  std::string_view message) noexcept
  {
    assignDiagnosticNoThrow(diagnostic, message);
    return slide::Status::Numerical_failure;
  }

  const JsonValue *jsonPath(const JsonValue &root,
                            std::initializer_list<std::string_view>
                              path)
  {
    const JsonValue *current = &root;
    for (const auto name : path) {
      if (current->kind != JsonValue::Kind::object)
        return nullptr;
      const auto found = current->object.find(name);
      if (found == current->object.end())
        return nullptr;
      current = &found->second;
    }
    return current;
  }

  std::optional<real_t> jsonConstant(const JsonValue &value,
                                     std::string &diagnostic)
  {
    if (value.kind == JsonValue::Kind::number)
      return value.number;
    if (value.kind != JsonValue::Kind::string)
      return std::nullopt;
    constexpr std::array sample_points{ 0.0, 0.5, 1.0 };
    std::array<real_t, 3> samples{};
    if (!detail::evaluateBpxExpressionSamples(
          value.string, sample_points, samples, diagnostic))
      return std::nullopt;
    const real_t tolerance =
      64.0 * std::numeric_limits<real_t>::epsilon()
      * std::max(std::numeric_limits<real_t>::min(), std::abs(samples[0]));
    if (std::abs(samples[1] - samples[0]) > tolerance
        || std::abs(samples[2] - samples[0]) > tolerance) {
      diagnostic = "state-dependent BPX diffusivity is not supported by the constant-D SPM composition";
      return std::nullopt;
    }
    return samples[0];
  }

  std::optional<OCVCurve> jsonCurve(const JsonValue &value,
                                    std::string &diagnostic)
  {
    if (value.kind == JsonValue::Kind::number)
      return OCVCurve{ { 0.0, 1.0 }, { value.number, value.number } };
    if (value.kind == JsonValue::Kind::string) {
      OCVCurve curve;
      if (!detail::sampleBpxExpressionCurve(
            value.string, curve, diagnostic))
        return std::nullopt;
      return curve;
    }
    if (value.kind != JsonValue::Kind::object)
      return std::nullopt;
    const auto x = value.object.find("x");
    const auto y = value.object.find("y");
    if (x == value.object.end() || y == value.object.end()
        || x->second.kind != JsonValue::Kind::array
        || y->second.kind != JsonValue::Kind::array
        || x->second.array.size() != y->second.array.size())
      return std::nullopt;
    OCVCurve curve;
    for (std::size_t i = 0; i < x->second.array.size(); ++i) {
      if (x->second.array[i].kind != JsonValue::Kind::number
          || y->second.array[i].kind != JsonValue::Kind::number)
        return std::nullopt;
      curve.stoichiometry.push_back(x->second.array[i].number);
      curve.value.push_back(y->second.array[i].number);
    }
    return detail::validParameterCurve(curve)
             ? std::optional<OCVCurve>{ std::move(curve) }
             : std::nullopt;
  }

} // namespace

slide::Status ParameterSet::fromBpxJson(std::string_view json,
                                        ParameterSet &output,
                                        std::string &diagnostic)
try {
  if (json.size() > max_bpx_json_bytes) {
    diagnostic = "BPX JSON exceeds 4194304 bytes";
    return slide::Status::Invalid_parameters;
  }
  JsonValue root;
  if (!detail::parseStrictJson(json, root, diagnostic)
      || root.kind != JsonValue::Kind::object) {
    if (diagnostic.empty()) diagnostic = "BPX root must be an object";
    return slide::Status::Invalid_parameters;
  }
  const auto *version = jsonPath(root, { "Header", "BPX" });
  const auto *model = jsonPath(root, { "Header", "Model" });
  std::string version_text;
  if (version != nullptr && version->kind == JsonValue::Kind::string
      && version->string.starts_with("1."))
    version_text = version->string;
  else if (version != nullptr && version->kind == JsonValue::Kind::number
           && version->number >= 1.0 && version->number < 2.0)
    version_text = "1.0 (legacy numeric header)";
  if (version_text.empty() || model == nullptr
      || model->kind != JsonValue::Kind::string
      || (model->string != "SPM" && model->string != "SPMe"
          && model->string != "DFN" && model->string != "Partial")) {
    diagnostic =
      "expected BPX 1.x Header with Model SPM, SPMe, DFN, or Partial";
    return slide::Status::Invalid_parameters;
  }
  ParameterSet candidate;
  auto addScalar = [&](std::initializer_list<std::string_view> path,
                       std::string_view name,
                       bool required = true) -> slide::Status {
    const auto *value = jsonPath(root, path);
    if (value == nullptr)
      return required ? slide::Status::Invalid_parameters
                      : slide::Status::Success;
    if (value->kind != JsonValue::Kind::number)
      return slide::Status::Invalid_parameters;
    return candidate.set(
      std::string{ name }, value->number, "BPX " + version_text);
  };
  const auto requireCellScalar = [&](std::string_view bpx_name,
                                     std::string_view parameter_name) {
    const auto status = addScalar(
      { "Parameterisation", "Cell", bpx_name }, parameter_name);
    if (status != slide::Status::Success)
      diagnostic = "missing, non-numeric, or invalid BPX Cell parameter: "
                   + std::string{ parameter_name };
    return status;
  };
  for (const auto &[bpx_name, parameter_name] : {
         std::pair{ std::string_view{ "Electrode area [m2]" },
                    std::string_view{ "Electrode area [m2]" } },
         std::pair{ std::string_view{ "Nominal cell capacity [A.h]" },
                    std::string_view{ "Nominal cell capacity [A.h]" } },
         std::pair{ std::string_view{ "Reference temperature [K]" },
                    std::string_view{ "Reference temperature [K]" } },
       }) {
    const auto status = requireCellScalar(bpx_name, parameter_name);
    if (status != slide::Status::Success)
      return status;
  }
  const auto optionalScalar = [&](std::initializer_list<std::string_view> path,
                                  std::string_view name) {
    const auto status = addScalar(path, name, false);
    if (status != slide::Status::Success)
      diagnostic = "invalid optional BPX scalar: " + std::string{ name };
    return status;
  };
  for (const auto &[bpx_name, parameter_name] : {
         std::pair{ std::string_view{ "External surface area [m2]" },
                    std::string_view{ "Cell cooling surface area [m2]" } },
         std::pair{ std::string_view{ "Volume [m3]" },
                    std::string_view{ "Cell volume [m3]" } },
         std::pair{ std::string_view{ "Density [kg.m-3]" },
                    std::string_view{ "Cell density [kg.m-3]" } },
         std::pair{ std::string_view{ "Specific heat capacity [J.K-1.kg-1]" },
                    std::string_view{ "Cell specific heat capacity [J.kg-1.K-1]" } },
       }) {
    const auto status = optionalScalar(
      { "Parameterisation", "Cell", bpx_name }, parameter_name);
    if (status != slide::Status::Success)
      return status;
  }

  struct BpxElectrode
  {
    std::string_view section;
    std::string_view prefix;
    std::string_view concentration_name;
    std::string_view diffusivity_name;
    std::string_view ocp_name;
  };
  constexpr std::array electrodes{
    BpxElectrode{ "Negative electrode", "Negative", "Maximum concentration in negative electrode [mol.m-3]", "Negative particle diffusivity [m2.s-1]", "Negative electrode OCP [V]" },
    BpxElectrode{ "Positive electrode", "Positive", "Maximum concentration in positive electrode [mol.m-3]", "Positive particle diffusivity [m2.s-1]", "Positive electrode OCP [V]" },
  };
  for (const auto &electrode : electrodes) {
    const std::string prefix{ electrode.prefix };
    const auto path = [&](std::string_view field) {
      return jsonPath(
        root, { "Parameterisation", electrode.section, field });
    };
    auto requireNumber = [&](std::string_view field, std::string name) {
      const auto *value = path(field);
      if (value == nullptr || value->kind != JsonValue::Kind::number)
        return slide::Status::Invalid_parameters;
      return candidate.set(
        std::move(name), value->number, "BPX " + version_text);
    };
    auto requireConstant = [&](std::string_view field, std::string name) {
      const auto *value = path(field);
      if (value == nullptr)
        return slide::Status::Invalid_parameters;
      auto constant = jsonConstant(*value, diagnostic);
      if (!constant.has_value())
        return slide::Status::Invalid_parameters;
      return candidate.set(
        std::move(name), *constant, "BPX " + version_text);
    };
    const auto requireElectrode = [&](slide::Status status) {
      if (status != slide::Status::Success && diagnostic.empty())
        diagnostic = "missing or unsupported BPX electrode scalar";
      return status;
    };
    auto status = requireElectrode(
      requireNumber("Thickness [m]", prefix + " electrode thickness [m]"));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireNumber(
      "Minimum stoichiometry", prefix + " electrode minimum stoichiometry"));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireNumber(
      "Maximum stoichiometry", prefix + " electrode maximum stoichiometry"));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireNumber(
      "Maximum concentration [mol.m-3]",
      std::string{ electrode.concentration_name }));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireNumber(
      "Particle radius [m]", prefix + " particle radius [m]"));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireConstant(
      "Diffusivity [m2.s-1]", std::string{ electrode.diffusivity_name }));
    if (status != slide::Status::Success)
      return status;
    status = requireElectrode(requireNumber(
      "Reaction rate constant [mol.m-2.s-1]",
      prefix + " electrode reaction rate constant [mol.m-2.s-1]"));
    if (status != slide::Status::Success)
      return status;
    const auto *area = path("Surface area per unit volume [m-1]");
    const auto *radius = path("Particle radius [m]");
    if (area == nullptr || radius == nullptr
        || area->kind != JsonValue::Kind::number
        || radius->kind != JsonValue::Kind::number) {
      diagnostic = "BPX electrode needs numeric surface area and particle radius";
      return slide::Status::Invalid_parameters;
    }
    const real_t fraction = area->number * radius->number / 3.0;
    const auto *porosity = path("Porosity");
    if (porosity != nullptr && porosity->kind != JsonValue::Kind::number) {
      diagnostic = "invalid optional BPX electrode porosity";
      return slide::Status::Invalid_parameters;
    }
    const real_t porosity_value = porosity != nullptr
                                    ? porosity->number
                                    : 1.0 - fraction;
    const std::string porosity_provenance = porosity != nullptr
                                              ? "BPX " + version_text
                                              : "SPM complement of BPX a*R/3";
    status = candidate.set(prefix + " electrode active material volume fraction",
                           fraction,
                           "derived exactly from BPX a*R/3");
    if (status != slide::Status::Success) {
      diagnostic = "invalid BPX derived active fraction";
      return status;
    }
    status = candidate.set(
      prefix + " electrode porosity", porosity_value, porosity_provenance);
    if (status != slide::Status::Success) {
      diagnostic = "invalid BPX electrode porosity";
      return status;
    }
    const auto *ocp = path("OCP [V]");
    const auto curve =
      ocp == nullptr ? std::nullopt : jsonCurve(*ocp, diagnostic);
    if (!curve.has_value()) {
      if (diagnostic.empty())
        diagnostic = "BPX OCP must be a numeric constant, function, or exact {x,y} table";
      return slide::Status::Invalid_parameters;
    }
    status = candidate.set(std::string{ electrode.ocp_name },
                           *curve,
                           "BPX " + version_text + " canonical curve");
    if (status != slide::Status::Success) {
      diagnostic = "failed to store BPX OCP canonical curve";
      return status;
    }
    const auto *activation =
      path("Reaction rate constant activation energy [J.mol-1]");
    if (activation != nullptr) {
      if (activation->kind != JsonValue::Kind::number) {
        diagnostic = "invalid BPX reaction-rate activation energy";
        return slide::Status::Invalid_parameters;
      }
      status = candidate.set(
        prefix + " electrode reaction rate activation energy [J.mol-1]",
        activation->number,
        "BPX " + version_text);
      if (status != slide::Status::Success) {
        diagnostic = "failed to store BPX reaction-rate activation energy";
        return status;
      }
    }
    const auto *diffusion_activation =
      path("Diffusivity activation energy [J.mol-1]");
    if (diffusion_activation != nullptr) {
      if (diffusion_activation->kind != JsonValue::Kind::number) {
        diagnostic = "invalid BPX diffusivity activation energy";
        return slide::Status::Invalid_parameters;
      }
      status = candidate.set(
        prefix + " particle diffusivity activation energy [J.mol-1]",
        diffusion_activation->number,
        "BPX " + version_text);
      if (status != slide::Status::Success) {
        diagnostic = "failed to store BPX diffusivity activation energy";
        return status;
      }
    }
  }

  auto setDefault = [&](std::string name, real_t value) -> slide::Status {
    return candidate.contains(name)
             ? slide::Status::Success
             : candidate.set(std::move(name), value, "BPX SPM default");
  };
  for (const auto &[section, bpx_name, parameter_name] : {
         std::tuple{ std::string_view{ "Initial conditions" },
                     std::string_view{ "Initial state-of-charge" },
                     std::string_view{ "Initial state-of-charge" } },
         std::tuple{ std::string_view{ "Initial conditions" },
                     std::string_view{ "Initial temperature [K]" },
                     std::string_view{ "Initial temperature [K]" } },
         std::tuple{ std::string_view{ "Initial conditions" },
                     std::string_view{ "Initial electrolyte concentration [mol.m-3]" },
                     std::string_view{ "Initial concentration in electrolyte [mol.m-3]" } },
         std::tuple{ std::string_view{ "Thermal environment" },
                     std::string_view{ "Ambient temperature [K]" },
                     std::string_view{ "Ambient temperature [K]" } },
         std::tuple{ std::string_view{ "Thermal environment" },
                     std::string_view{ "Heat transfer coefficient [W.m-2.K-1]" },
                     std::string_view{ "Total heat transfer coefficient [W.m-2.K-1]" } },
       }) {
    const auto status = optionalScalar(
      { "State", section, bpx_name }, parameter_name);
    if (status != slide::Status::Success)
      return status;
  }
  const real_t reference_temperature =
    *candidate.findScalar("Reference temperature [K]");
  for (const auto &[name, value] : {
         std::pair{ std::string_view{ "Initial state-of-charge" }, 0.5 },
         std::pair{ std::string_view{ "Initial temperature [K]" },
                    reference_temperature },
         std::pair{ std::string_view{ "Ambient temperature [K]" },
                    reference_temperature },
         std::pair{ std::string_view{ "Initial concentration in electrolyte [mol.m-3]" },
                    1000.0 },
         std::pair{ std::string_view{ "Initial SEI thickness [m]" }, 1e-9 },
         std::pair{ std::string_view{ "Contact resistance [Ohm]" }, 0.0 },
       }) {
    const auto status = setDefault(std::string{ name }, value);
    if (status != slide::Status::Success) {
      diagnostic = "failed to install BPX SPM default: "
                   + std::string{ name };
      return status;
    }
  }
  diagnostic.clear();
  output = std::move(candidate);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return allocationFailure(diagnostic, "BPX JSON allocation failed");
} catch (const std::length_error &) {
  return allocationFailure(diagnostic, "BPX JSON size is not representable");
}

slide::Status ParameterSet::fromBpxFile(const std::filesystem::path &path,
                                        ParameterSet &output,
                                        std::string &diagnostic)
try {
  std::string contents;
  const auto read =
    detail::readBoundedFile(path, max_bpx_json_bytes, contents);
  if (read != detail::BoundedFileRead::success) {
    if (read == detail::BoundedFileRead::open_failed)
      diagnostic = "could not open BPX file";
    else if (read == detail::BoundedFileRead::too_large)
      diagnostic = "BPX JSON exceeds 4194304 bytes";
    else
      diagnostic = "could not read complete BPX file";
    return detail::boundedFileStatus(read);
  }
  return fromBpxJson(contents, output, diagnostic);
} catch (const std::bad_alloc &) {
  return allocationFailure(diagnostic, "BPX file allocation failed");
} catch (const std::length_error &) {
  return allocationFailure(diagnostic, "BPX file size is not representable");
}

} // namespace slide::core
