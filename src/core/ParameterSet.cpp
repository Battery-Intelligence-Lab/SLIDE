/**
 * @file ParameterSet.cpp
 * @brief PyBaMM 26.6.2.0 Chen2020 and BPX 1.x absorption.
 */

#include "ParameterSet.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <new>
#include <optional>
#include <queue>
#include <utility>

namespace slide::core {
namespace {

  bool validCurve(const OCVCurve &curve)
  {
    if (curve.stoichiometry.size() != curve.value.size()
        || curve.stoichiometry.size() < 2)
      return false;
    for (std::size_t i = 0; i < curve.stoichiometry.size(); ++i)
      if (!is_finite(curve.stoichiometry[i]) || !is_finite(curve.value[i])
          || (i > 0 && curve.stoichiometry[i] <= curve.stoichiometry[i - 1]))
        return false;
    return true;
  }

  real_t graphiteOcp(real_t x)
  {
    return 1.9793 * std::exp(-39.3631 * x) + 0.2482
           - 0.0909 * std::tanh(29.8538 * (x - 0.1234))
           - 0.04478 * std::tanh(14.9159 * (x - 0.2769))
           - 0.0205 * std::tanh(30.4444 * (x - 0.6103));
  }

  real_t nmcOcp(real_t x)
  {
    return -0.8090 * x + 4.4875
           - 0.0428 * std::tanh(18.5138 * (x - 0.5542))
           - 17.7326 * std::tanh(15.7890 * (x - 0.3117))
           + 17.5842 * std::tanh(15.9308 * (x - 0.3120));
  }

  template <class Function>
  OCVCurve sampleCurve(Function function)
  {
    struct Segment
    {
      real_t left;
      real_t right;
      real_t left_value;
      real_t right_value;
      real_t relative_error;
    };
    const auto make_segment = [&function](real_t left, real_t right, real_t left_value, real_t right_value) {
      real_t error{};
      for (const real_t fraction : { real_t{ 0.25 }, real_t{ 0.5 }, real_t{ 0.75 } }) {
        const real_t x = left + fraction * (right - left);
        const real_t exact = function(x);
        const real_t linear = left_value + fraction * (right_value - left_value);
        error = std::max(error, std::abs(linear - exact) / std::max(std::abs(exact), real_t{ 1e-12 }));
      }
      return Segment{ left, right, left_value, right_value, error };
    };
    const auto lower_error = [](const Segment &left, const Segment &right) {
      return left.relative_error < right.relative_error;
    };

    constexpr std::size_t maximum_points = 4096;
    // Quarter-point checks plus this safety factor bound interpolation between
    // the probes without spending the table budget in traversal order.
    constexpr real_t relative_tolerance = 1e-7;

    std::priority_queue<Segment, std::vector<Segment>, decltype(lower_error)>
      work{ lower_error };
    work.push(make_segment(0.0, 1.0, function(0.0), function(1.0)));
    while (work.size() + 1 < maximum_points
           && work.top().relative_error > relative_tolerance) {
      const Segment segment = work.top();
      work.pop();
      const real_t middle = 0.5 * (segment.left + segment.right);
      const real_t middle_value = function(middle);
      work.push(make_segment(segment.left, middle, segment.left_value, middle_value));
      work.push(make_segment(middle, segment.right, middle_value, segment.right_value));
    }

    std::vector<Segment> segments;
    segments.reserve(work.size());
    while (!work.empty()) {
      segments.push_back(work.top());
      work.pop();
    }
    std::sort(segments.begin(), segments.end(), [](const Segment &left, const Segment &right) {
      return left.left < right.left;
    });
    OCVCurve curve;
    curve.stoichiometry.reserve(segments.size() + 1);
    curve.value.reserve(segments.size() + 1);
    for (const auto &segment : segments) {
      curve.stoichiometry.push_back(segment.left);
      curve.value.push_back(segment.left_value);
    }
    curve.stoichiometry.push_back(segments.back().right);
    curve.value.push_back(segments.back().right_value);
    return curve;
  }

  const real_t *requiredScalar(const ParameterSet &parameters,
                               std::string_view name)
  {
    return parameters.findScalar(name);
  }

} // namespace

std::string ParameterSet::canonicalName(std::string_view name)
{
  static constexpr std::array aliases{
    std::pair{ std::string_view{ "Negative electrode diffusivity [m2.s-1]" },
               std::string_view{ "Negative particle diffusivity [m2.s-1]" } },
    std::pair{ std::string_view{ "Positive electrode diffusivity [m2.s-1]" },
               std::string_view{ "Positive particle diffusivity [m2.s-1]" } },
    std::pair{ std::string_view{ "Exchange-current density for lithium plating [A.m-2]" },
               std::string_view{ "Exchange-current density for lithium metal electrode [A.m-2]" } },
    std::pair{ std::string_view{ "1 + dlnf/dlnc" },
               std::string_view{ "Thermodynamic factor" } },
  };
  for (const auto &[old_name, current_name] : aliases)
    if (name == old_name)
      return std::string{ current_name };
  return std::string{ name };
}

slide::Status ParameterSet::set(std::string name, ParameterValue value,
                                std::string provenance)
{
  name = canonicalName(name);
  if (name.empty() || provenance.empty())
    return slide::Status::Invalid_parameters;
  if (const auto *scalar = std::get_if<real_t>(&value)) {
    if (!is_finite(*scalar))
      return slide::Status::Invalid_parameters;
  } else if (!validCurve(std::get<OCVCurve>(value))) {
    return slide::Status::Invalid_parameters;
  }
  try {
    values_.insert_or_assign(std::move(name),
                             Entry{ std::move(value), std::move(provenance) });
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  }
  return slide::Status::Success;
}

slide::Status ParameterSet::update(
  std::span<const ParameterDescription> values)
{
  try {
    ParameterSet candidate = *this;
    for (const auto &description : values) {
      const auto status = candidate.set(description.name, description.value, description.provenance);
      if (status != slide::Status::Success)
        return status;
    }
    *this = std::move(candidate);
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  }
  return slide::Status::Success;
}

bool ParameterSet::contains(std::string_view name) const
{
  return find(name) != nullptr;
}

const ParameterValue *ParameterSet::find(std::string_view name) const
{
  const auto canonical = canonicalName(name);
  const auto found = values_.find(canonical);
  return found == values_.end() ? nullptr : &found->second.value;
}

const real_t *ParameterSet::findScalar(std::string_view name) const
{
  const auto *value = find(name);
  return value == nullptr ? nullptr : std::get_if<real_t>(value);
}

const OCVCurve *ParameterSet::findCurve(std::string_view name) const
{
  const auto *value = find(name);
  return value == nullptr ? nullptr : std::get_if<OCVCurve>(value);
}

std::vector<ParameterDescription> ParameterSet::describe() const
{
  std::vector<ParameterDescription> result;
  result.reserve(values_.size());
  for (const auto &[name, entry] : values_)
    result.push_back({ name, entry.value, entry.provenance });
  return result;
}

slide::Status ParameterSet::chen2020(ParameterSet &output)
{
  // Exact scalar entries from PyBaMM 26.6.2.0's Chen2020.py. Function-valued
  // OCPs are canonicalised below to the D-16 4096-knot table form.
  static constexpr std::array scalars{
    std::pair{ "Ratio of lithium moles to SEI moles", 2.0 },
    std::pair{ "SEI partial molar volume [m3.mol-1]", 9.585e-5 },
    std::pair{ "SEI reaction exchange current density [A.m-2]", 1.5e-7 },
    std::pair{ "SEI resistivity [Ohm.m]", 200000.0 },
    std::pair{ "SEI solvent diffusivity [m2.s-1]", 2.5e-22 },
    std::pair{ "Initial SEI thickness [m]", 5e-9 },
    std::pair{ "Negative current collector thickness [m]", 1.2e-5 },
    std::pair{ "Negative electrode thickness [m]", 8.52e-5 },
    std::pair{ "Separator thickness [m]", 1.2e-5 },
    std::pair{ "Positive electrode thickness [m]", 7.56e-5 },
    std::pair{ "Positive current collector thickness [m]", 1.6e-5 },
    std::pair{ "Electrode height [m]", 0.065 },
    std::pair{ "Electrode width [m]", 1.58 },
    std::pair{ "Electrode area [m2]", 0.065 * 1.58 },
    std::pair{ "Cell cooling surface area [m2]", 0.00531 },
    std::pair{ "Cell volume [m3]", 2.42e-5 },
    std::pair{ "Nominal cell capacity [A.h]", 5.0 },
    std::pair{ "Current function [A]", 5.0 },
    std::pair{ "Contact resistance [Ohm]", 0.0 },
    std::pair{ "Maximum concentration in negative electrode [mol.m-3]", 33133.0 },
    std::pair{ "Negative particle diffusivity [m2.s-1]", 3.3e-14 },
    std::pair{ "Negative electrode porosity", 0.25 },
    std::pair{ "Negative electrode active material volume fraction", 0.75 },
    std::pair{ "Negative particle radius [m]", 5.86e-6 },
    std::pair{ "Negative electrode OCP entropic change [V.K-1]", 0.0 },
    std::pair{ "Negative electrode minimum stoichiometry", 0.02634579027064577 },
    std::pair{ "Negative electrode maximum stoichiometry", 0.910618046652409 },
    std::pair{ "Negative electrode reaction rate constant [mol.m-2.s-1]", 6.48e-7 / 96487.0 },
    std::pair{ "Negative electrode reaction rate activation energy [J.mol-1]", 35000.0 },
    std::pair{ "Maximum concentration in positive electrode [mol.m-3]", 63104.0 },
    std::pair{ "Positive particle diffusivity [m2.s-1]", 4e-15 },
    std::pair{ "Positive electrode porosity", 0.335 },
    std::pair{ "Positive electrode active material volume fraction", 0.665 },
    std::pair{ "Positive particle radius [m]", 5.22e-6 },
    std::pair{ "Positive electrode OCP entropic change [V.K-1]", 0.0 },
    std::pair{ "Positive electrode minimum stoichiometry", 0.2638452245913301 },
    std::pair{ "Positive electrode maximum stoichiometry", 0.853974674630047 },
    std::pair{ "Positive electrode reaction rate constant [mol.m-2.s-1]", 3.42e-6 / 96487.0 },
    std::pair{ "Positive electrode reaction rate activation energy [J.mol-1]", 17800.0 },
    std::pair{ "Separator porosity", 0.47 },
    std::pair{ "Initial concentration in electrolyte [mol.m-3]", 1000.0 },
    std::pair{ "Cation transference number", 0.2594 },
    std::pair{ "Thermodynamic factor", 1.0 },
    std::pair{ "Reference temperature [K]", 298.15 },
    std::pair{ "Total heat transfer coefficient [W.m-2.K-1]", 10.0 },
    std::pair{ "Ambient temperature [K]", 298.15 },
    std::pair{ "Lower voltage cut-off [V]", 2.5 },
    std::pair{ "Upper voltage cut-off [V]", 4.2 },
    std::pair{ "Initial concentration in negative electrode [mol.m-3]", 29866.0 },
    std::pair{ "Initial concentration in positive electrode [mol.m-3]", 17038.0 },
    std::pair{ "Initial temperature [K]", 298.15 },
  };

  ParameterSet candidate;
  for (const auto &[name, value] : scalars) {
    const auto status = candidate.set(name, value, "PyBaMM 26.6.2.0 Chen2020");
    if (status != slide::Status::Success)
      return status;
  }
  auto status = candidate.set("Negative electrode OCP [V]",
                              sampleCurve(graphiteOcp),
                              "PyBaMM 26.6.2.0 Chen2020 function; adaptive D-16 table");
  if (status != slide::Status::Success)
    return status;
  status = candidate.set("Positive electrode OCP [V]",
                         sampleCurve(nmcOcp),
                         "PyBaMM 26.6.2.0 Chen2020 function; adaptive D-16 table");
  if (status != slide::Status::Success)
    return status;
  output = std::move(candidate);
  return slide::Status::Success;
}

slide::Status ParameterSet::toSpmInput(SpmFactoryInput &output) const
{
  constexpr std::array required_names{
    "Nominal cell capacity [A.h]",
    "Electrode area [m2]",
    "Initial concentration in electrolyte [mol.m-3]",
    "Reference temperature [K]",
    "Initial temperature [K]",
    "Initial SEI thickness [m]",
  };
  for (const auto name : required_names)
    if (requiredScalar(*this, name) == nullptr)
      return slide::Status::Invalid_parameters;
  const auto *negative_ocp = findCurve("Negative electrode OCP [V]");
  const auto *positive_ocp = findCurve("Positive electrode OCP [V]");
  if (negative_ocp == nullptr || positive_ocp == nullptr)
    return slide::Status::Invalid_parameters;

  SpmFactoryInput candidate;
  candidate.design.capacity_Ah = *findScalar("Nominal cell capacity [A.h]");
  candidate.design.electrode_area = *findScalar("Electrode area [m2]");
  candidate.design.electrolyte.concentration = *findScalar(
    "Initial concentration in electrolyte [mol.m-3]");
  candidate.design.thermal.reference_temperature = *findScalar("Reference temperature [K]");
  candidate.design.thermal.environment_temperature = findScalar("Ambient temperature [K]") != nullptr
                                                       ? *findScalar("Ambient temperature [K]")
                                                       : candidate.design.thermal.reference_temperature;
  candidate.design.thermal.volume = findScalar("Cell volume [m3]") != nullptr
                                      ? *findScalar("Cell volume [m3]")
                                      : 1.0;
  candidate.design.thermal.surface_area = findScalar("Cell cooling surface area [m2]") != nullptr
                                            ? *findScalar("Cell cooling surface area [m2]")
                                            : 1.0;
  candidate.design.thermal.h_conv = findScalar("Total heat transfer coefficient [W.m-2.K-1]") != nullptr
                                      ? *findScalar("Total heat transfer coefficient [W.m-2.K-1]")
                                      : 0.0;
  candidate.design.thermal.density = 1626.0;
  if (const auto *density = findScalar("Cell density [kg.m-3]"))
    candidate.design.thermal.density = *density;
  candidate.design.thermal.heat_capacity = 750.0;
  if (const auto *heat_capacity = findScalar("Cell specific heat capacity [J.kg-1.K-1]"))
    candidate.design.thermal.heat_capacity = *heat_capacity;
  candidate.initial_temperature = *findScalar("Initial temperature [K]");
  candidate.initial_sei_thickness = *findScalar("Initial SEI thickness [m]");
  candidate.sei_resistivity_area = findScalar("SEI resistivity [Ohm.m]") != nullptr
                                     ? *findScalar("SEI resistivity [Ohm.m]")
                                     : 0.0;
  const real_t contact = findScalar("Contact resistance [Ohm]") != nullptr
                           ? *findScalar("Contact resistance [Ohm]")
                           : 0.0;
  candidate.initial_current_collector_resistance = contact * candidate.design.electrode_area;

  struct ElectrodeNames
  {
    Domain domain;
    std::string_view prefix;
    std::string_view concentration;
    std::string_view diffusivity;
    std::string_view ocp;
  };
  constexpr std::array electrode_names{
    ElectrodeNames{ Domain::neg, "Negative", "Maximum concentration in negative electrode [mol.m-3]", "Negative particle diffusivity [m2.s-1]", "Negative electrode OCP [V]" },
    ElectrodeNames{ Domain::pos, "Positive", "Maximum concentration in positive electrode [mol.m-3]", "Positive particle diffusivity [m2.s-1]", "Positive electrode OCP [V]" },
  };
  for (const auto &names : electrode_names) {
    auto &electrode = domain_value(candidate.design.electrode, names.domain);
    const std::string prefix{ names.prefix };
    const auto scalar = [&](std::string suffix) {
      return findScalar(prefix + std::move(suffix));
    };
    const auto *thickness = scalar(" electrode thickness [m]");
    const auto *porosity = scalar(" electrode porosity");
    const auto *fraction = scalar(" electrode active material volume fraction");
    const auto *radius = scalar(" particle radius [m]");
    const auto *cs_max = findScalar(names.concentration);
    const auto *diffusivity = findScalar(names.diffusivity);
    const auto *diffusivity_activation = scalar(
      " particle diffusivity activation energy [J.mol-1]");
    const auto *minimum = scalar(" electrode minimum stoichiometry");
    const auto *maximum = scalar(" electrode maximum stoichiometry");
    const auto *reaction = scalar(" electrode reaction rate constant [mol.m-2.s-1]");
    const auto *reaction_activation = scalar(" electrode reaction rate activation energy [J.mol-1]");
    if (thickness == nullptr || porosity == nullptr || fraction == nullptr
        || radius == nullptr || cs_max == nullptr || diffusivity == nullptr
        || minimum == nullptr || maximum == nullptr || reaction == nullptr)
      return slide::Status::Invalid_parameters;
    electrode.thickness = *thickness;
    electrode.porosity = *porosity;
    electrode.active_fraction = *fraction;
    electrode.particle_radius = *radius;
    electrode.active_material.cs_max = *cs_max;
    electrode.active_material.x_0 = names.domain == Domain::neg ? *minimum : *maximum;
    electrode.active_material.x_100 = names.domain == Domain::neg ? *maximum : *minimum;
    electrode.active_material.D_s = { .reference_value = *diffusivity,
                                      .activation_energy = diffusivity_activation == nullptr
                                                             ? 0.0
                                                             : *diffusivity_activation,
                                      .reference_temperature = candidate.design.thermal.reference_temperature };
    electrode.active_material.k_ct = { .reference_value = *reaction,
                                       .activation_energy = reaction_activation == nullptr ? 0.0 : *reaction_activation,
                                       .reference_temperature = candidate.design.thermal.reference_temperature };
    electrode.active_material.ocv = *findCurve(names.ocp);
  }
  const auto &negative = domain_value(candidate.design.electrode, Domain::neg).active_material;
  if (const auto *initial_soc = findScalar("Initial state-of-charge")) {
    candidate.initial_soc = *initial_soc;
  } else {
    const auto *initial_concentration = findScalar(
      "Initial concentration in negative electrode [mol.m-3]");
    if (initial_concentration == nullptr)
      return slide::Status::Invalid_parameters;
    const real_t negative_initial = *initial_concentration
                                    / *findScalar("Maximum concentration in negative electrode [mol.m-3]");
    candidate.initial_soc = (negative_initial - negative.x_0)
                            / (negative.x_100 - negative.x_0);
  }
  constexpr std::array zero_x{ 0.0, 1.0 };
  constexpr std::array zero_y{ 0.0, 0.0 };
  candidate.total_entropic_coefficient = { std::vector<real_t>{ zero_x.begin(), zero_x.end() },
                                           std::vector<real_t>{ zero_y.begin(), zero_y.end() } };
  candidate.negative_entropic_coefficient = candidate.total_entropic_coefficient;
  output = std::move(candidate);
  return slide::Status::Success;
}

namespace {

  struct JsonValue
  {
    enum class Kind : unsigned char { null_value,
                                      boolean,
                                      number,
                                      string,
                                      array,
                                      object } kind{ Kind::null_value };
    bool boolean{};
    real_t number{};
    std::string string{};
    std::vector<JsonValue> array{};
    std::map<std::string, JsonValue, std::less<>> object{};
  };

  /** Strict evaluator for BPX's documented one-variable expression language. */
  class BpxExpression
  {
  public:
    bool compile(std::string_view source, std::string &diagnostic)
    {
      source_ = source;
      cursor_ = 0;
      nodes_.clear();
      diagnostic.clear();
      if (source.empty() || source.size() > 65'536)
        return fail(diagnostic, "BPX expression is empty or too long");
      root_ = parseExpression(diagnostic, 0);
      skipSpace();
      if (root_ < 0 || cursor_ != source_.size()) {
        if (diagnostic.empty())
          fail(diagnostic, "unexpected BPX expression token");
        return false;
      }
      return true;
    }

    bool evaluate(real_t x, real_t &output) const
    {
      return root_ >= 0 && evaluateNode(root_, x, output, 0)
             && is_finite(output);
    }

  private:
    enum class Kind : unsigned char {
      literal,
      variable,
      add,
      subtract,
      multiply,
      divide,
      power,
      negate,
      exponential,
      hyperbolic_tangent,
      hyperbolic_cosine
    };

    struct Node
    {
      Kind kind{ Kind::literal };
      real_t value{};
      int left{ -1 };
      int right{ -1 };
    };

    int append(Node node, std::string &diagnostic)
    {
      if (nodes_.size() >= 1024) {
        fail(diagnostic, "BPX expression exceeds 1024 operations");
        return -1;
      }
      nodes_.push_back(node);
      return static_cast<int>(nodes_.size() - 1);
    }

    int parseExpression(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      int left = parseTerm(diagnostic, depth + 1);
      while (left >= 0) {
        skipSpace();
        Kind kind;
        if (take('+'))
          kind = Kind::add;
        else if (take('-'))
          kind = Kind::subtract;
        else
          break;
        const int right = parseTerm(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = kind, .left = left, .right = right }, diagnostic);
      }
      return left;
    }

    int parseTerm(std::string &diagnostic, int depth)
    {
      int left = parseUnary(diagnostic, depth + 1);
      while (left >= 0) {
        skipSpace();
        Kind kind;
        if (source_.substr(cursor_).starts_with("**"))
          break;
        if (take('*'))
          kind = Kind::multiply;
        else if (take('/'))
          kind = Kind::divide;
        else
          break;
        const int right = parseUnary(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = kind, .left = left, .right = right }, diagnostic);
      }
      return left;
    }

    int parseUnary(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      skipSpace();
      if (take('+'))
        return parseUnary(diagnostic, depth + 1);
      if (take('-')) {
        const int child = parseUnary(diagnostic, depth + 1);
        return child < 0 ? -1
                         : append({ .kind = Kind::negate, .left = child }, diagnostic);
      }
      return parsePower(diagnostic, depth + 1);
    }

    int parsePower(std::string &diagnostic, int depth)
    {
      int left = parsePrimary(diagnostic, depth + 1);
      skipSpace();
      if (left >= 0 && consume("**")) {
        const int right = parseUnary(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = Kind::power, .left = left, .right = right }, diagnostic);
      }
      return left;
    }

    int parsePrimary(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      skipSpace();
      if (take('(')) {
        const int result = parseExpression(diagnostic, depth + 1);
        skipSpace();
        if (result < 0 || !take(')')) {
          fail(diagnostic, "expected ')' in BPX expression");
          return -1;
        }
        return result;
      }
      if (cursor_ < source_.size()
          && (std::isdigit(static_cast<unsigned char>(source_[cursor_]))
              || source_[cursor_] == '.')) {
        const char *first = source_.data() + cursor_;
        const char *last = source_.data() + source_.size();
        real_t value{};
        const auto parsed = std::from_chars(first, last, value, std::chars_format::general);
        if (parsed.ec != std::errc{} || parsed.ptr == first || !is_finite(value)) {
          fail(diagnostic, "invalid number in BPX expression");
          return -1;
        }
        cursor_ = static_cast<std::size_t>(parsed.ptr - source_.data());
        return append({ .kind = Kind::literal, .value = value }, diagnostic);
      }
      if (cursor_ < source_.size()
          && std::isalpha(static_cast<unsigned char>(source_[cursor_]))) {
        const std::size_t begin = cursor_++;
        while (cursor_ < source_.size()
               && std::isalnum(static_cast<unsigned char>(source_[cursor_])))
          ++cursor_;
        const auto identifier = source_.substr(begin, cursor_ - begin);
        if (identifier == "x")
          return append({ .kind = Kind::variable }, diagnostic);
        Kind kind;
        if (identifier == "exp")
          kind = Kind::exponential;
        else if (identifier == "tanh")
          kind = Kind::hyperbolic_tangent;
        else if (identifier == "cosh")
          kind = Kind::hyperbolic_cosine;
        else {
          fail(diagnostic, "unsupported BPX expression identifier");
          return -1;
        }
        skipSpace();
        if (!take('(')) {
          fail(diagnostic, "expected '(' after BPX function");
          return -1;
        }
        const int child = parseExpression(diagnostic, depth + 1);
        skipSpace();
        if (child < 0 || !take(')')) {
          fail(diagnostic, "expected ')' after BPX function argument");
          return -1;
        }
        return append({ .kind = kind, .left = child }, diagnostic);
      }
      fail(diagnostic, "expected value in BPX expression");
      return -1;
    }

    bool evaluateNode(int index, real_t x, real_t &output, int depth) const
    {
      if (depth > 128 || index < 0
          || static_cast<std::size_t>(index) >= nodes_.size())
        return false;
      const auto &node = nodes_[static_cast<std::size_t>(index)];
      if (node.kind == Kind::literal) {
        output = node.value;
        return true;
      }
      if (node.kind == Kind::variable) {
        output = x;
        return is_finite(output);
      }
      real_t left{};
      if (!evaluateNode(node.left, x, left, depth + 1))
        return false;
      if (node.kind == Kind::negate)
        output = -left;
      else if (node.kind == Kind::exponential)
        output = std::exp(left);
      else if (node.kind == Kind::hyperbolic_tangent)
        output = std::tanh(left);
      else if (node.kind == Kind::hyperbolic_cosine)
        output = std::cosh(left);
      else {
        real_t right{};
        if (!evaluateNode(node.right, x, right, depth + 1))
          return false;
        switch (node.kind) {
        case Kind::add:
          output = left + right;
          break;
        case Kind::subtract:
          output = left - right;
          break;
        case Kind::multiply:
          output = left * right;
          break;
        case Kind::divide:
          if (right == 0.0)
            return false;
          output = left / right;
          break;
        case Kind::power:
          output = std::pow(left, right);
          break;
        default:
          return false;
        }
      }
      return is_finite(output);
    }

    void skipSpace()
    {
      while (cursor_ < source_.size()
             && std::isspace(static_cast<unsigned char>(source_[cursor_])))
        ++cursor_;
    }

    bool take(char token)
    {
      if (cursor_ < source_.size() && source_[cursor_] == token) {
        ++cursor_;
        return true;
      }
      return false;
    }

    bool consume(std::string_view token)
    {
      if (!source_.substr(cursor_).starts_with(token))
        return false;
      cursor_ += token.size();
      return true;
    }

    bool fail(std::string &diagnostic, std::string_view message)
    {
      if (diagnostic.empty())
        diagnostic = std::string{ message } + " at expression byte "
                     + std::to_string(cursor_);
      return false;
    }

    std::string_view source_{};
    std::size_t cursor_{};
    std::vector<Node> nodes_{};
    int root_{ -1 };
  };

  class JsonParser
  {
  public:
    JsonParser(std::string_view source, std::string &diagnostic)
      : source_{ source }, diagnostic_{ diagnostic }
    {}

    bool parse(JsonValue &output)
    {
      diagnostic_.clear();
      skipSpace();
      if (!parseValue(output, 0))
        return false;
      skipSpace();
      if (cursor_ != source_.size())
        return fail("unexpected trailing JSON data");
      return true;
    }

  private:
    bool parseValue(JsonValue &output, int depth)
    {
      if (depth > 64)
        return fail("JSON nesting exceeds 64 levels");
      skipSpace();
      if (cursor_ == source_.size())
        return fail("unexpected end of JSON");
      const char token = source_[cursor_];
      if (token == '{') return parseObject(output, depth + 1);
      if (token == '[') return parseArray(output, depth + 1);
      if (token == '"') {
        output.kind = JsonValue::Kind::string;
        return parseString(output.string);
      }
      if (token == '-' || (token >= '0' && token <= '9'))
        return parseNumber(output);
      if (consume("true")) {
        output.kind = JsonValue::Kind::boolean;
        output.boolean = true;
        return true;
      }
      if (consume("false")) {
        output.kind = JsonValue::Kind::boolean;
        output.boolean = false;
        return true;
      }
      if (consume("null")) {
        output.kind = JsonValue::Kind::null_value;
        return true;
      }
      return fail("invalid JSON value");
    }

    bool parseObject(JsonValue &output, int depth)
    {
      ++cursor_;
      output.kind = JsonValue::Kind::object;
      skipSpace();
      if (take('}')) return true;
      while (true) {
        std::string key;
        if (!parseString(key))
          return false;
        skipSpace();
        if (!take(':'))
          return fail("expected ':' after object key");
        JsonValue value;
        if (!parseValue(value, depth))
          return false;
        if (!output.object.emplace(std::move(key), std::move(value)).second)
          return fail("duplicate JSON object key");
        skipSpace();
        if (take('}')) return true;
        if (!take(','))
          return fail("expected ',' or '}' in object");
        skipSpace();
      }
    }

    bool parseArray(JsonValue &output, int depth)
    {
      ++cursor_;
      output.kind = JsonValue::Kind::array;
      skipSpace();
      if (take(']')) return true;
      while (true) {
        JsonValue value;
        if (!parseValue(value, depth))
          return false;
        output.array.push_back(std::move(value));
        skipSpace();
        if (take(']')) return true;
        if (!take(','))
          return fail("expected ',' or ']' in array");
        skipSpace();
      }
    }

    bool parseString(std::string &output)
    {
      skipSpace();
      if (!take('"'))
        return fail("expected JSON string");
      output.clear();
      while (cursor_ < source_.size()) {
        const unsigned char c = static_cast<unsigned char>(source_[cursor_++]);
        if (c == '"') return true;
        if (c < 0x20)
          return fail("control byte in JSON string");
        if (c != '\\') {
          output.push_back(static_cast<char>(c));
          continue;
        }
        if (cursor_ == source_.size())
          return fail("unfinished JSON escape");
        const char escape = source_[cursor_++];
        switch (escape) {
        case '"':
          output.push_back('"');
          break;
        case '\\':
          output.push_back('\\');
          break;
        case '/':
          output.push_back('/');
          break;
        case 'b':
          output.push_back('\b');
          break;
        case 'f':
          output.push_back('\f');
          break;
        case 'n':
          output.push_back('\n');
          break;
        case 'r':
          output.push_back('\r');
          break;
        case 't':
          output.push_back('\t');
          break;
        case 'u': {
          std::uint32_t codepoint{};
          if (!parseHex4(codepoint))
            return false;
          if (codepoint >= 0xd800U && codepoint <= 0xdbffU) {
            if (cursor_ + 2 > source_.size() || source_[cursor_] != '\\'
                || source_[cursor_ + 1] != 'u')
              return fail("high surrogate without low surrogate");
            cursor_ += 2;
            std::uint32_t low{};
            if (!parseHex4(low))
              return false;
            if (low < 0xdc00U || low > 0xdfffU)
              return fail("invalid low surrogate");
            codepoint = 0x10000U + ((codepoint - 0xd800U) << 10U)
                        + (low - 0xdc00U);
          } else if (codepoint >= 0xdc00U && codepoint <= 0xdfffU) {
            return fail("unpaired low surrogate");
          }
          appendUtf8(output, codepoint);
          break;
        }
        default:
          return fail("unsupported JSON escape (use UTF-8 directly)");
        }
      }
      return fail("unterminated JSON string");
    }

    bool parseNumber(JsonValue &output)
    {
      const char *first = source_.data() + cursor_;
      const char *last = source_.data() + source_.size();
      real_t value{};
      const auto parsed = std::from_chars(first, last, value, std::chars_format::general);
      if (parsed.ec != std::errc{} || parsed.ptr == first || !is_finite(value))
        return fail("invalid JSON number");
      cursor_ = static_cast<std::size_t>(parsed.ptr - source_.data());
      output.kind = JsonValue::Kind::number;
      output.number = value;
      return true;
    }

    bool parseHex4(std::uint32_t &value)
    {
      if (cursor_ + 4 > source_.size())
        return fail("unfinished Unicode escape");
      value = 0;
      for (int i = 0; i < 4; ++i) {
        const char c = source_[cursor_++];
        unsigned digit{};
        if (c >= '0' && c <= '9')
          digit = static_cast<unsigned>(c - '0');
        else if (c >= 'a' && c <= 'f')
          digit = 10U + static_cast<unsigned>(c - 'a');
        else if (c >= 'A' && c <= 'F')
          digit = 10U + static_cast<unsigned>(c - 'A');
        else
          return fail("invalid hexadecimal Unicode escape");
        value = (value << 4U) | digit;
      }
      return true;
    }

    static void appendUtf8(std::string &output, std::uint32_t codepoint)
    {
      if (codepoint <= 0x7fU) {
        output.push_back(static_cast<char>(codepoint));
      } else if (codepoint <= 0x7ffU) {
        output.push_back(static_cast<char>(0xc0U | (codepoint >> 6U)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      } else if (codepoint <= 0xffffU) {
        output.push_back(static_cast<char>(0xe0U | (codepoint >> 12U)));
        output.push_back(static_cast<char>(0x80U | ((codepoint >> 6U) & 0x3fU)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      } else {
        output.push_back(static_cast<char>(0xf0U | (codepoint >> 18U)));
        output.push_back(static_cast<char>(0x80U | ((codepoint >> 12U) & 0x3fU)));
        output.push_back(static_cast<char>(0x80U | ((codepoint >> 6U) & 0x3fU)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      }
    }

    void skipSpace()
    {
      while (cursor_ < source_.size()
             && std::isspace(static_cast<unsigned char>(source_[cursor_])))
        ++cursor_;
    }

    bool take(char expected)
    {
      if (cursor_ < source_.size() && source_[cursor_] == expected) {
        ++cursor_;
        return true;
      }
      return false;
    }

    bool consume(std::string_view text)
    {
      if (!source_.substr(cursor_).starts_with(text))
        return false;
      cursor_ += text.size();
      return true;
    }

    bool fail(std::string_view message)
    {
      diagnostic_ = std::string{ message } + " at byte " + std::to_string(cursor_);
      return false;
    }

    std::string_view source_{};
    std::string &diagnostic_;
    std::size_t cursor_{};
  };

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
    BpxExpression expression;
    if (!expression.compile(value.string, diagnostic))
      return std::nullopt;
    std::array<real_t, 3> samples{};
    if (!expression.evaluate(0.0, samples[0])
        || !expression.evaluate(0.5, samples[1])
        || !expression.evaluate(1.0, samples[2])) {
      diagnostic = "BPX scalar expression is non-finite";
      return std::nullopt;
    }
    const real_t tolerance = 64.0 * std::numeric_limits<real_t>::epsilon()
                             * std::max(std::numeric_limits<real_t>::min(),
                                        std::abs(samples[0]));
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
      BpxExpression expression;
      if (!expression.compile(value.string, diagnostic))
        return std::nullopt;
      const auto curve = sampleCurve([&expression](real_t x) {
        real_t value_at_x{};
        return expression.evaluate(x, value_at_x)
                 ? value_at_x
                 : std::numeric_limits<real_t>::quiet_NaN();
      });
      if (!validCurve(curve)) {
        diagnostic = "BPX function is non-finite on stoichiometry [0,1]";
        return std::nullopt;
      }
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
    return validCurve(curve) ? std::optional<OCVCurve>{ std::move(curve) }
                             : std::nullopt;
  }

} // namespace

slide::Status ParameterSet::fromBpxJson(std::string_view json,
                                        ParameterSet &output,
                                        std::string &diagnostic)
{
  JsonValue root;
  JsonParser parser{ json, diagnostic };
  if (!parser.parse(root) || root.kind != JsonValue::Kind::object) {
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
    diagnostic = "expected BPX 1.x Header with Model SPM, SPMe, DFN, or Partial";
    return slide::Status::Invalid_parameters;
  }
  ParameterSet candidate;
  auto addScalar = [&](std::initializer_list<std::string_view> path,
                       std::string_view name,
                       bool required = true) {
    const auto *value = jsonPath(root, path);
    if (value == nullptr)
      return !required;
    return value->kind == JsonValue::Kind::number
           && candidate.set(std::string{ name }, value->number, "BPX " + version_text)
                == slide::Status::Success;
  };
  const auto requireCellScalar = [&](std::string_view bpx_name,
                                     std::string_view parameter_name) {
    if (addScalar({ "Parameterisation", "Cell", bpx_name }, parameter_name))
      return true;
    diagnostic = "missing, non-numeric, or invalid BPX Cell parameter: "
                 + std::string{ parameter_name };
    return false;
  };
  if (!requireCellScalar("Electrode area [m2]", "Electrode area [m2]")
      || !requireCellScalar("Nominal cell capacity [A.h]",
                            "Nominal cell capacity [A.h]")
      || !requireCellScalar("Reference temperature [K]",
                            "Reference temperature [K]"))
    return slide::Status::Invalid_parameters;
  addScalar({ "Parameterisation", "Cell", "External surface area [m2]" },
            "Cell cooling surface area [m2]",
            false);
  addScalar({ "Parameterisation", "Cell", "Volume [m3]" },
            "Cell volume [m3]",
            false);
  addScalar({ "Parameterisation", "Cell", "Density [kg.m-3]" },
            "Cell density [kg.m-3]",
            false);
  addScalar({ "Parameterisation", "Cell", "Specific heat capacity [J.K-1.kg-1]" },
            "Cell specific heat capacity [J.kg-1.K-1]",
            false);

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
      return jsonPath(root, { "Parameterisation", electrode.section, field });
    };
    auto requireNumber = [&](std::string_view field, std::string name) {
      const auto *value = path(field);
      return value != nullptr && value->kind == JsonValue::Kind::number
             && candidate.set(std::move(name), value->number, "BPX " + version_text)
                  == slide::Status::Success;
    };
    auto requireConstant = [&](std::string_view field, std::string name) {
      const auto *value = path(field);
      if (value == nullptr)
        return false;
      auto constant = jsonConstant(*value, diagnostic);
      return constant.has_value()
             && candidate.set(std::move(name), *constant, "BPX " + version_text)
                  == slide::Status::Success;
    };
    if (!requireNumber("Thickness [m]", prefix + " electrode thickness [m]")
        || !requireNumber("Minimum stoichiometry", prefix + " electrode minimum stoichiometry")
        || !requireNumber("Maximum stoichiometry", prefix + " electrode maximum stoichiometry")
        || !requireNumber("Maximum concentration [mol.m-3]", std::string{ electrode.concentration_name })
        || !requireNumber("Particle radius [m]", prefix + " particle radius [m]")
        || !requireConstant("Diffusivity [m2.s-1]", std::string{ electrode.diffusivity_name })
        || !requireNumber("Reaction rate constant [mol.m-2.s-1]",
                          prefix + " electrode reaction rate constant [mol.m-2.s-1]")) {
      if (diagnostic.empty())
        diagnostic = "missing or unsupported BPX electrode scalar";
      return slide::Status::Invalid_parameters;
    }
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
    const real_t porosity_value = porosity != nullptr
                                      && porosity->kind == JsonValue::Kind::number
                                    ? porosity->number
                                    : 1.0 - fraction;
    const std::string porosity_provenance = porosity != nullptr
                                                && porosity->kind == JsonValue::Kind::number
                                              ? "BPX " + version_text
                                              : "SPM complement of BPX a*R/3";
    if (candidate.set(prefix + " electrode active material volume fraction", fraction, "derived exactly from BPX a*R/3") != slide::Status::Success
        || candidate.set(prefix + " electrode porosity", porosity_value, porosity_provenance) != slide::Status::Success) {
      diagnostic = "invalid BPX derived active fraction";
      return slide::Status::Invalid_parameters;
    }
    const auto *ocp = path("OCP [V]");
    const auto curve = ocp == nullptr ? std::nullopt : jsonCurve(*ocp, diagnostic);
    if (!curve.has_value()
        || candidate.set(std::string{ electrode.ocp_name }, *curve, "BPX " + version_text + " canonical curve")
             != slide::Status::Success) {
      if (diagnostic.empty())
        diagnostic = "BPX OCP must be a numeric constant, function, or exact {x,y} table";
      return slide::Status::Invalid_parameters;
    }
    const auto *activation = path("Reaction rate constant activation energy [J.mol-1]");
    if (activation != nullptr
        && (activation->kind != JsonValue::Kind::number
            || candidate.set(prefix + " electrode reaction rate activation energy [J.mol-1]",
                             activation->number,
                             "BPX " + version_text)
                 != slide::Status::Success)) {
      diagnostic = "invalid BPX reaction-rate activation energy";
      return slide::Status::Invalid_parameters;
    }
    const auto *diffusion_activation = path("Diffusivity activation energy [J.mol-1]");
    if (diffusion_activation != nullptr
        && (diffusion_activation->kind != JsonValue::Kind::number
            || candidate.set(prefix + " particle diffusivity activation energy [J.mol-1]",
                             diffusion_activation->number,
                             "BPX " + version_text)
                 != slide::Status::Success)) {
      diagnostic = "invalid BPX diffusivity activation energy";
      return slide::Status::Invalid_parameters;
    }
  }

  auto setDefault = [&](std::string name, real_t value) {
    return candidate.contains(name)
           || candidate.set(std::move(name), value, "BPX SPM default")
                == slide::Status::Success;
  };
  addScalar({ "State", "Initial conditions", "Initial state-of-charge" },
            "Initial state-of-charge",
            false);
  addScalar({ "State", "Initial conditions", "Initial temperature [K]" },
            "Initial temperature [K]",
            false);
  addScalar({ "State", "Initial conditions", "Initial electrolyte concentration [mol.m-3]" },
            "Initial concentration in electrolyte [mol.m-3]",
            false);
  addScalar({ "State", "Thermal environment", "Ambient temperature [K]" },
            "Ambient temperature [K]",
            false);
  addScalar({ "State", "Thermal environment", "Heat transfer coefficient [W.m-2.K-1]" },
            "Total heat transfer coefficient [W.m-2.K-1]",
            false);
  if (!setDefault("Initial state-of-charge", 0.5)
      || !setDefault("Initial temperature [K]",
                     *candidate.findScalar("Reference temperature [K]"))
      || !setDefault("Ambient temperature [K]",
                     *candidate.findScalar("Reference temperature [K]"))
      || !setDefault("Initial concentration in electrolyte [mol.m-3]", 1000.0)
      || !setDefault("Initial SEI thickness [m]", 1e-9)
      || !setDefault("Contact resistance [Ohm]", 0.0)) {
    diagnostic = "failed to install BPX SPM defaults";
    return slide::Status::Numerical_failure;
  }
  diagnostic.clear();
  output = std::move(candidate);
  return slide::Status::Success;
}

slide::Status ParameterSet::fromBpxFile(const std::filesystem::path &path,
                                        ParameterSet &output,
                                        std::string &diagnostic)
{
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    diagnostic = "could not open BPX file";
    return slide::Status::Invalid_parameters;
  }
  try {
    const std::string contents{ std::istreambuf_iterator<char>{ input },
                                std::istreambuf_iterator<char>{} };
    return fromBpxJson(contents, output, diagnostic);
  } catch (const std::bad_alloc &) {
    diagnostic = "BPX file is too large";
    return slide::Status::Numerical_failure;
  }
}

} // namespace slide::core
