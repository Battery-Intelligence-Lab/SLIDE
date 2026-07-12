/**
 * @file core_ParameterSet_test.cpp
 * @brief Phase-7 Chen2020/BPX parameter fidelity and M0.7 split-file gates.
 */

#include "../support/RecordedBits.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../../src/core/detail/StrictJson.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

using namespace slide;

namespace {

double graphiteOcp(double x)
{
  return 1.9793 * std::exp(-39.3631 * x) + 0.2482
         - 0.0909 * std::tanh(29.8538 * (x - 0.1234))
         - 0.04478 * std::tanh(14.9159 * (x - 0.2769))
         - 0.0205 * std::tanh(30.4444 * (x - 0.6103));
}

double interpolate(const core::OCVCurve &curve, double x)
{
  const auto upper = std::upper_bound(curve.stoichiometry.begin(),
                                      curve.stoichiometry.end(),
                                      x);
  const auto i = static_cast<std::size_t>(upper - curve.stoichiometry.begin() - 1);
  const double fraction = (x - curve.stoichiometry[i])
                          / (curve.stoichiometry[i + 1] - curve.stoichiometry[i]);
  return curve.value[i] + fraction * (curve.value[i + 1] - curve.value[i]);
}

void replaceRequired(std::string &text,
                     std::string_view needle,
                     std::string_view replacement)
{
  const auto position = text.find(needle);
  REQUIRE(position != std::string::npos);
  text.replace(position, needle.size(), replacement);
}

core::ParameterSet copyWithout(
  const core::ParameterSet &source,
  std::initializer_list<std::string_view>
    omitted)
{
  core::ParameterSet result;
  for (const auto &description : source.describe()) {
    if (std::find(omitted.begin(), omitted.end(), description.name)
        != omitted.end())
      continue;
    REQUIRE(result.set(description.name,
                       description.value,
                       description.provenance)
            == Status::Success);
  }
  return result;
}

constexpr std::string_view bpx_fixture = R"json({
  "Header": {"BPX": "1.0.0", "Title": "LGM50 \u03bc fixture", "Model": "SPM"},
  "Parameterisation": {
    "Cell": {
      "Electrode area [m2]": 0.1027,
      "Number of electrode pairs connected in parallel to make a cell": 1,
      "Lower voltage cut-off [V]": 2.5,
      "Upper voltage cut-off [V]": 4.2,
      "Nominal cell capacity [A.h]": 5.0,
      "Reference temperature [K]": 298.15,
      "External surface area [m2]": 0.00531,
      "Volume [m3]": 2.42e-5
    },
    "Negative electrode": {
      "Thickness [m]": 8.52e-5,
      "Minimum stoichiometry": 0.02634579027064577,
      "Maximum stoichiometry": 0.910618046652409,
      "Maximum concentration [mol.m-3]": 33133.0,
      "Particle radius [m]": 5.86e-6,
      "Surface area per unit volume [m-1]": 383959.0443686007,
      "Diffusivity [m2.s-1]": 3.3e-14,
      "OCP [V]": {"x": [0.0, 0.5, 1.0], "y": [1.9793, 0.25, 0.1]},
      "Reaction rate constant [mol.m-2.s-1]": 6.716e-12,
      "Reaction rate constant activation energy [J.mol-1]": 35000.0
    },
    "Positive electrode": {
      "Thickness [m]": 7.56e-5,
      "Minimum stoichiometry": 0.2638452245913301,
      "Maximum stoichiometry": 0.853974674630047,
      "Maximum concentration [mol.m-3]": 63104.0,
      "Particle radius [m]": 5.22e-6,
      "Surface area per unit volume [m-1]": 382183.908045977,
      "Diffusivity [m2.s-1]": 4e-15,
      "OCP [V]": {"x": [0.0, 0.5, 1.0], "y": [4.5, 3.8, 3.0]},
      "Reaction rate constant [mol.m-2.s-1]": 3.544e-11,
      "Reaction rate constant activation energy [J.mol-1]": 17800.0
    }
  },
  "State": {
    "Initial conditions": {
      "Initial state-of-charge": 0.75,
      "Initial temperature [K]": 299.15,
      "Initial electrolyte concentration [mol.m-3]": 1000.0
    },
    "Thermal environment": {
      "Ambient temperature [K]": 298.15,
      "Heat transfer coefficient [W.m-2.K-1]": 10.0
    }
  }
})json";

/** Independent byte recorder for ordered names, kinds, and provenance. */
struct RecordedText
{
  std::uint64_t fnv1a{ UINT64_C(14695981039346656037) };
  std::uint64_t mixed{ UINT64_C(0x6a09e667f3bcc909) };
  std::size_t strings{};
  std::size_t bytes{};

  void append(std::string_view value) noexcept
  {
    absorb(static_cast<std::uint64_t>(value.size()));
    for (const unsigned char byte : value)
      absorb(byte);
    ++strings;
    bytes += value.size();
  }

private:
  void absorb(std::uint64_t word) noexcept
  {
    constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);
    for (unsigned shift = 0; shift < 64; shift += 8) {
      fnv1a ^= (word >> shift) & UINT64_C(0xff);
      fnv1a *= fnv_prime;
    }
    mixed ^= word + UINT64_C(0x9e3779b97f4a7c15) + std::rotl(mixed, 17);
    mixed *= UINT64_C(0xbf58476d1ce4e5b9);
    mixed ^= mixed >> 29;
  }
};

struct RecordedParameters
{
  test_support::RecordedBits values{};
  RecordedText metadata{};
  std::size_t entries{};
};

class ScopedTemporaryDirectory
{
public:
  ScopedTemporaryDirectory()
  {
    const auto nonce = static_cast<std::uint64_t>(
      std::chrono::steady_clock::now().time_since_epoch().count());
    const auto root = std::filesystem::temp_directory_path();
    for (std::uint64_t attempt = 0; attempt < 1024; ++attempt) {
      const auto candidate = root
                             / ("slide_bpx_parser_"
                                + std::to_string(nonce) + "_"
                                + std::to_string(attempt));
      std::error_code error;
      if (std::filesystem::create_directory(candidate, error)) {
        path_ = candidate;
        return;
      }
      if (error)
        throw std::runtime_error{ "cannot create BPX test directory: "
                                  + error.message() };
    }
    throw std::runtime_error{ "cannot reserve a unique BPX test directory" };
  }

  ~ScopedTemporaryDirectory()
  {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
  }

  ScopedTemporaryDirectory(const ScopedTemporaryDirectory &) = delete;
  ScopedTemporaryDirectory &operator=(const ScopedTemporaryDirectory &) = delete;

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  std::filesystem::path path_{};
};

RecordedParameters recordParameters(const core::ParameterSet &parameters)
{
  RecordedParameters recorded;
  const auto descriptions = parameters.describe();
  recorded.entries = descriptions.size();
  for (const auto &description : descriptions) {
    recorded.metadata.append(description.name);
    recorded.metadata.append(
      std::holds_alternative<core::real_t>(description.value) ? "scalar" : "curve");
    recorded.metadata.append(description.provenance);
    if (const auto *scalar = std::get_if<core::real_t>(&description.value)) {
      const std::array singleton{ *scalar };
      recorded.values.append(singleton);
    } else {
      const auto &curve = std::get<core::OCVCurve>(description.value);
      recorded.values.append(curve.stoichiometry);
      recorded.values.append(curve.value);
    }
  }
  return recorded;
}

void recordCurve(test_support::RecordedBits &recorded,
                 const core::OCVCurve &curve)
{
  recorded.append(curve.stoichiometry);
  recorded.append(curve.value);
}

void recordScalar(test_support::RecordedBits &recorded, double value)
{
  const std::array singleton{ value };
  recorded.append(singleton);
}

void recordArrhenius(test_support::RecordedBits &recorded,
                     const core::Arrhenius &value)
{
  recordScalar(recorded, value.reference_value);
  recordScalar(recorded, value.activation_energy);
  recordScalar(recorded, value.reference_temperature);
}

/**
 * Explicit declaration-order traversal of every public SpmFactoryInput numeric
 * field and curve. Empty vectors are appended too, so shape is part of the
 * trace rather than being inferred from the floating-point value count.
 */
test_support::RecordedBits recordSpmInput(const core::SpmFactoryInput &input)
{
  test_support::RecordedBits recorded;
  for (const core::Domain domain : core::domains) {
    const auto &electrode = input.design.electrode[core::domain_index(domain)];
    const auto &material = electrode.active_material;
    recordCurve(recorded, material.ocv);
    recordScalar(recorded, material.cs_max);
    recordScalar(recorded, material.x_0);
    recordScalar(recorded, material.x_100);
    recordArrhenius(recorded, material.D_s);
    recordArrhenius(recorded, material.k_ct);
    recordScalar(recorded, electrode.thickness);
    recordScalar(recorded, electrode.porosity);
    recordScalar(recorded, electrode.active_fraction);
    recordScalar(recorded, electrode.particle_radius);
    recordScalar(recorded, electrode.stress.youngs_modulus);
    recordScalar(recorded, electrode.stress.poisson_ratio);
    recordScalar(recorded, electrode.stress.partial_molar_volume);
    recordScalar(recorded, static_cast<double>(electrode.aging.size()));
    for (const auto &mechanism : electrode.aging) {
      recordScalar(recorded, static_cast<double>(mechanism.kind));
      recorded.append(mechanism.coefficients);
    }
  }
  recordScalar(recorded, input.design.separator.thickness);
  recordScalar(recorded, input.design.separator.porosity);
  recordScalar(recorded, input.design.electrolyte.concentration);
  recordScalar(recorded, input.design.electrolyte.diffusivity);
  recordScalar(recorded, input.design.electrolyte.transference_number);
  recordScalar(recorded, input.design.thermal.density);
  recordScalar(recorded, input.design.thermal.heat_capacity);
  recordScalar(recorded, input.design.thermal.volume);
  recordScalar(recorded, input.design.thermal.surface_area);
  recordScalar(recorded, input.design.thermal.h_conv);
  recordScalar(recorded, input.design.thermal.reference_temperature);
  recordScalar(recorded, input.design.thermal.environment_temperature);
  recordScalar(recorded, input.design.capacity_Ah);
  recordScalar(recorded, input.design.electrode_area);

  recordCurve(recorded, input.total_entropic_coefficient);
  recordCurve(recorded, input.negative_entropic_coefficient);
  recordCurve(recorded, input.negative_laresgoiti_stress);
  recordScalar(recorded, input.initial_soc);
  recordScalar(recorded, input.initial_temperature);
  recordScalar(recorded, input.initial_sei_thickness);
  recordScalar(recorded, input.initial_lost_lithium);
  recordScalar(recorded, input.initial_crack_surface_fraction);
  recordScalar(recorded, input.initial_plated_lithium_thickness);
  for (const auto value : input.initial_specific_resistance)
    recordScalar(recorded, value);
  recordScalar(recorded, input.initial_current_collector_resistance);
  recordScalar(recorded, input.initial_stress_interval);
  recordScalar(recorded, input.sei_resistivity_area);

  const auto &sei = input.sei;
  recordScalar(recorded, static_cast<double>(sei.model_mask));
  recordScalar(recorded, sei.reduce_active_fraction ? 1.0 : 0.0);
  for (const auto value : {
         sei.F,
         sei.Rg,
         sei.n,
         sei.n_sei,
         sei.alpha_sei,
         sei.reference_temperature,
         sei.electrode_area,
         sei.negative_particle_radius,
         sei.sei_resistivity_area,
         sei.sei_equilibrium_potential,
         sei.sei_molar_volume,
         sei.electrolyte_reactant_concentration,
         sei.main_molar_volume,
         sei.side_molar_volume,
         sei.porosity_coefficient,
         sei.model1_k,
         sei.model1_k_activation,
         sei.model2_k,
         sei.model2_k_activation,
         sei.model2_D,
         sei.model2_D_activation,
         sei.model3_k,
         sei.model3_k_activation,
         sei.model3_D,
         sei.model3_D_activation,
         sei.model4_k,
         sei.model4_k_activation,
         sei.model4_D,
         sei.model4_D_activation,
       })
    recordScalar(recorded, value);

  const auto &crack = input.surface_crack;
  recordScalar(recorded, static_cast<double>(crack.model_mask));
  recordScalar(recorded, crack.reduce_negative_diffusivity ? 1.0 : 0.0);
  for (const auto value : {
         crack.F,
         crack.Rg,
         crack.n_sei,
         crack.alpha_sei,
         crack.reference_temperature,
         crack.electrode_area,
         crack.negative_cs_max,
         crack.sei_resistivity_area,
         crack.sei_equilibrium_potential,
         crack.model1_alpha,
         crack.model2_alpha,
         crack.model3_alpha,
         crack.model4_alpha,
         crack.model4_max_surface,
         crack.model5_k,
         crack.model5_k_activation,
         crack.diffusion_exponent,
       })
    recordScalar(recorded, value);

  const auto &lam = input.lam;
  recordScalar(recorded, static_cast<double>(lam.model_mask));
  recordScalar(recorded, lam.F);
  recordScalar(recorded, lam.Rg);
  recordScalar(recorded, lam.n);
  recordScalar(recorded, lam.reference_temperature);
  for (const auto value : lam.particle_radius)
    recordScalar(recorded, value);
  for (const auto value : lam.model1_stress_coefficient)
    recordScalar(recorded, value);
  for (const auto value : lam.model2_linear_flux)
    recordScalar(recorded, value);
  for (const auto value : lam.model2_sqrt_flux)
    recordScalar(recorded, value);
  recordScalar(recorded, lam.model2_activation);
  recordScalar(recorded, lam.model3_k);
  recordScalar(recorded, lam.model3_k_activation);
  recordScalar(recorded, lam.model3_equilibrium_potential);
  for (const auto value : lam.model4_area_coefficient)
    recordScalar(recorded, value);
  recordScalar(recorded, lam.positive_ocv.valid() ? 1.0 : 0.0);
  std::vector<double> compiled_queries;
  std::vector<double> compiled_values;
  std::vector<double> compiled_derivatives;
  if (lam.positive_ocv.valid()) {
    recordScalar(recorded, lam.positive_ocv.x_min());
    recordScalar(recorded, lam.positive_ocv.x_max());
    recordScalar(recorded, static_cast<double>(lam.positive_ocv.knots()));
    const auto &source = input.design.electrode[core::domain_index(core::Domain::pos)]
                           .active_material.ocv.stoichiometry;
    compiled_queries.reserve(source.size() * 2U - 1U);
    for (std::size_t i = 0; i < source.size(); ++i) {
      compiled_queries.push_back(source[i]);
      if (i + 1 < source.size())
        compiled_queries.push_back(source[i] + 0.5 * (source[i + 1] - source[i]));
    }
    compiled_values.reserve(compiled_queries.size());
    compiled_derivatives.reserve(compiled_queries.size());
    for (const auto query : compiled_queries) {
      compiled_values.push_back(lam.positive_ocv.eval(query));
      compiled_derivatives.push_back(lam.positive_ocv.derivative(query));
    }
  }
  recorded.append(compiled_queries);
  recorded.append(compiled_values);
  recorded.append(compiled_derivatives);

  const auto &plating = input.lithium_plating;
  for (const auto value : {
         plating.F,
         plating.Rg,
         plating.n,
         plating.n_plating,
         plating.alpha_plating,
         plating.reference_temperature,
         plating.electrode_area,
         plating.sei_resistivity_area,
         plating.equilibrium_potential,
         plating.plated_lithium_molar_density,
         plating.reaction_rate_ref,
         plating.reaction_rate_activation,
       })
    recordScalar(recorded, value);
  return recorded;
}

constexpr std::string_view bpx_all_operators_expression =
  "0.25 + (0.4*x - 0.03)/(1.0 + x) + exp(-2*x) + 0.01*cosh(x) "
  "- 0.02*tanh(3*(x-0.5)) + 0.001*x**2**3 - 0.002*(-x**2) "
  "+ (+0.0001*x)";

std::string allOperatorsBpxFixture()
{
  std::string functional{ bpx_fixture };
  replaceRequired(functional,
                  "\"Diffusivity [m2.s-1]\": 3.3e-14",
                  "\"Diffusivity [m2.s-1]\": \"3.3e-14*(2+2)/4\",\n"
                  "      \"Diffusivity activation energy [J.mol-1]\": 30000.0,\n"
                  "      \"Porosity\": 0.25398416831491194");
  replaceRequired(
    functional,
    "\"OCP [V]\": {\"x\": [0.0, 0.5, 1.0], \"y\": [1.9793, 0.25, 0.1]}",
    std::string{ "\"OCP [V]\": \"" } + std::string{ bpx_all_operators_expression }
      + "\"");
  replaceRequired(
    functional,
    "\"OCP [V]\": {\"x\": [0.0, 0.5, 1.0], \"y\": [4.5, 3.8, 3.0]}",
    "\"OCP [V]\": \"4.5 - 1.5*x\"");
  return functional;
}

} // namespace

TEST_CASE("strict JSON replaces successful outputs and preserves rejected outputs",
          "[core][parameters][JSON][atomic][9C-3]")
{
  core::detail::StrictJsonValue output;
  output.kind = core::detail::StrictJsonValue::Kind::object;
  output.object["stale"].kind =
    core::detail::StrictJsonValue::Kind::boolean;
  output.object["stale"].boolean = true;

  std::string diagnostic{ "stale diagnostic" };
  REQUIRE(core::detail::parseStrictJson(
    R"json({"fresh":[1,2]})json", output, diagnostic));
  CHECK(diagnostic.empty());
  REQUIRE(output.kind == core::detail::StrictJsonValue::Kind::object);
  REQUIRE(output.object.size() == 1);
  REQUIRE(output.object.contains("fresh"));
  REQUIRE(output.object.at("fresh").kind
          == core::detail::StrictJsonValue::Kind::array);
  REQUIRE(output.object.at("fresh").array.size() == 2);
  CHECK(output.object.at("fresh").array[0].number == 1.0);
  CHECK(output.object.at("fresh").array[1].number == 2.0);

  CHECK_FALSE(core::detail::parseStrictJson(
    R"json({"corrupt":[3,]})json", output, diagnostic));
  CHECK_FALSE(diagnostic.empty());
  REQUIRE(output.kind == core::detail::StrictJsonValue::Kind::object);
  REQUIRE(output.object.size() == 1);
  REQUIRE(output.object.contains("fresh"));
  REQUIRE(output.object.at("fresh").array.size() == 2);
  CHECK(output.object.at("fresh").array[0].number == 1.0);
  CHECK(output.object.at("fresh").array[1].number == 2.0);
}

TEST_CASE("M0.7 pre-split ParameterSet behavior is bit-exact and fully framed",
          "[core][parameters][BPX][9C-3][recorded]")
{
  core::ParameterSet chen;
  REQUIRE(core::ParameterSet::chen2020(chen) == Status::Success);
  const auto chen_recorded = recordParameters(chen);
  core::SpmFactoryInput chen_input;
  REQUIRE(chen.toSpmInput(chen_input) == Status::Success);
  const auto input_recorded = recordSpmInput(chen_input);

  const std::string functional = allOperatorsBpxFixture();
  core::ParameterSet bpx;
  std::string diagnostic;
  REQUIRE(core::ParameterSet::fromBpxJson(functional, bpx, diagnostic)
          == Status::Success);
  REQUIRE(diagnostic.empty());
  const auto bpx_recorded = recordParameters(bpx);

  CAPTURE(chen_recorded.entries,
          chen_recorded.values.values,
          chen_recorded.values.fnv1a,
          chen_recorded.values.mixed,
          chen_recorded.metadata.strings,
          chen_recorded.metadata.bytes,
          chen_recorded.metadata.fnv1a,
          chen_recorded.metadata.mixed);
  CAPTURE(input_recorded.values,
          input_recorded.fnv1a,
          input_recorded.mixed);
  CAPTURE(bpx_recorded.entries,
          bpx_recorded.values.values,
          bpx_recorded.values.fnv1a,
          bpx_recorded.values.mixed,
          bpx_recorded.metadata.strings,
          bpx_recorded.metadata.bytes,
          bpx_recorded.metadata.fnv1a,
          bpx_recorded.metadata.mixed);

  CHECK(chen_recorded.entries == 53);
  CHECK(chen_recorded.values.values == 11'537);
  CHECK(chen_recorded.metadata.strings == 159);
  CHECK(chen_recorded.metadata.bytes == 3'491);
  CHECK(input_recorded.values == 11'635);
  CHECK(bpx_recorded.entries == 35);
  CHECK(bpx_recorded.values.values == 4'135);
  CHECK(bpx_recorded.metadata.strings == 105);
  CHECK(bpx_recorded.metadata.bytes == 1'906);
  CHECK(chen_recorded.metadata.fnv1a == UINT64_C(0x8ea717ef0459a7b7));
  CHECK(chen_recorded.metadata.mixed == UINT64_C(0x745379e447872b1f));
  CHECK(bpx_recorded.metadata.fnv1a == UINT64_C(0xf4a1c008bf74be9a));
  CHECK(bpx_recorded.metadata.mixed == UINT64_C(0xfa654b7fe7b9ab89));

#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_RELEASE) && defined(SLIDE_TEST_IPO)
  constexpr auto expected_chen_fnv = UINT64_C(0xe5c3dc16e932ad03);
  constexpr auto expected_chen_mixed = UINT64_C(0xb4e4cfcce5167740);
  constexpr auto expected_input_fnv = UINT64_C(0x705555174f1b2e23);
  constexpr auto expected_input_mixed = UINT64_C(0xd9357b2ee6af55e4);
  constexpr auto expected_bpx_fnv = UINT64_C(0xd3f393daa406ad80);
  constexpr auto expected_bpx_mixed = UINT64_C(0xf7d10dea65f3c8f5);
#elif defined(SLIDE_TEST_RELEASE)
  constexpr auto expected_chen_fnv = UINT64_C(0xe5c3dc16e932ad03);
  constexpr auto expected_chen_mixed = UINT64_C(0xb4e4cfcce5167740);
  constexpr auto expected_input_fnv = UINT64_C(0x705555174f1b2e23);
  constexpr auto expected_input_mixed = UINT64_C(0xd9357b2ee6af55e4);
  constexpr auto expected_bpx_fnv = UINT64_C(0xd3f393daa406ad80);
  constexpr auto expected_bpx_mixed = UINT64_C(0xf7d10dea65f3c8f5);
#else
  constexpr auto expected_chen_fnv = UINT64_C(0xa5d211a63dd77522);
  constexpr auto expected_chen_mixed = UINT64_C(0x8ccfefe72e03094e);
  constexpr auto expected_input_fnv = UINT64_C(0x134b5b4650709d1a);
  constexpr auto expected_input_mixed = UINT64_C(0x1e63d406b3d98297);
  constexpr auto expected_bpx_fnv = UINT64_C(0x04bdc30e3e6ca305);
  constexpr auto expected_bpx_mixed = UINT64_C(0xe75964cae8725037);
#endif
  CHECK(chen_recorded.values.fnv1a == expected_chen_fnv);
  CHECK(chen_recorded.values.mixed == expected_chen_mixed);
  CHECK(input_recorded.fnv1a == expected_input_fnv);
  CHECK(input_recorded.mixed == expected_input_mixed);
  CHECK(bpx_recorded.values.fnv1a == expected_bpx_fnv);
  CHECK(bpx_recorded.values.mixed == expected_bpx_mixed);
#endif

  core::ParameterSet rejected = chen;
  const auto before_rejection = recordParameters(rejected);
  const auto check_atomic_rejection = [&] {
    const auto retained = recordParameters(rejected);
    CHECK(retained.entries == before_rejection.entries);
    CHECK(retained.values.values == before_rejection.values.values);
    CHECK(retained.values.fnv1a == before_rejection.values.fnv1a);
    CHECK(retained.values.mixed == before_rejection.values.mixed);
    CHECK(retained.metadata.strings == before_rejection.metadata.strings);
    CHECK(retained.metadata.bytes == before_rejection.metadata.bytes);
    CHECK(retained.metadata.fnv1a == before_rejection.metadata.fnv1a);
    CHECK(retained.metadata.mixed == before_rejection.metadata.mixed);
  };

  diagnostic = "poison";
  CHECK(core::ParameterSet::fromBpxJson("{\"Header\":", rejected, diagnostic)
        == Status::Invalid_parameters);
  INFO(diagnostic);
  CHECK(diagnostic == "unexpected end of JSON at byte 10");
  check_atomic_rejection();

  std::string invalid_expression = functional;
  replaceRequired(invalid_expression, bpx_all_operators_expression, "sin(x)");
  diagnostic = "poison";
  CHECK(core::ParameterSet::fromBpxJson(
          invalid_expression, rejected, diagnostic)
        == Status::Invalid_parameters);
  INFO(diagnostic);
  CHECK(diagnostic
        == "unsupported BPX expression identifier at expression byte 3");
  check_atomic_rejection();

  std::string state_dependent_diffusivity = functional;
  replaceRequired(state_dependent_diffusivity,
                  "3.3e-14*(2+2)/4",
                  "3.3e-14*(1+x)");
  diagnostic = "poison";
  CHECK(core::ParameterSet::fromBpxJson(
          state_dependent_diffusivity, rejected, diagnostic)
        == Status::Invalid_parameters);
  INFO(diagnostic);
  CHECK(diagnostic
        == "state-dependent BPX diffusivity is not supported by the constant-D SPM composition");
  check_atomic_rejection();
}

TEST_CASE("P7-G3 Chen2020 scalar and curve absorption is traceable",
          "[core][parameters][Chen2020][P7-G3]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  REQUIRE(parameters.size() >= 50);
  CHECK(*parameters.findScalar("Nominal cell capacity [A.h]") == 5.0);
  CHECK(*parameters.findScalar("Negative electrode thickness [m]") == 8.52e-5);
  CHECK(*parameters.findScalar("Negative electrode diffusivity [m2.s-1]") == 3.3e-14);
  CHECK(*parameters.findScalar("1 + dlnf/dlnc") == 1.0);
  const auto descriptions = parameters.describe();
  REQUIRE(descriptions.size() == parameters.size());
  for (const auto &description : descriptions) {
    const auto *round_trip = parameters.find(description.name);
    REQUIRE(round_trip != nullptr);
    CHECK(*round_trip == description.value);
    CHECK_FALSE(description.provenance.empty());
  }

  const auto *curve = parameters.findCurve("Negative electrode OCP [V]");
  REQUIRE(curve != nullptr);
  REQUIRE(curve->stoichiometry.size() <= 4096);
  REQUIRE(curve->stoichiometry.size() > 100);
  double maximum_error{};
  for (int i = 1; i < 1000; ++i) {
    const double x = static_cast<double>(i) / 1000.0;
    maximum_error = std::max(maximum_error,
                             std::abs(interpolate(*curve, x) - graphiteOcp(x)));
  }
  CAPTURE(maximum_error);
  CHECK(maximum_error <= 1e-6);

  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  CHECK(input.design.capacity_Ah == 5.0);
  CHECK(input.initial_soc > 0.98);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
}

TEST_CASE("P7-G3 BPX 1.x tables absorb exactly and malformed input is atomic",
          "[core][parameters][BPX][P7-G3]")
{
  core::ParameterSet parameters;
  std::string diagnostic;
  const auto status =
    core::ParameterSet::fromBpxJson(bpx_fixture, parameters, diagnostic);
  INFO(diagnostic);
  REQUIRE(status == Status::Success);
  CHECK(diagnostic.empty());
  CHECK(*parameters.findScalar("Electrode area [m2]") == 0.1027);
  CHECK(*parameters.findScalar("Initial state-of-charge") == 0.75);
  CHECK(*parameters.findScalar("Negative electrode active material volume fraction")
        == 0.75);
  const auto *negative_ocp = parameters.findCurve("Negative electrode OCP [V]");
  REQUIRE(negative_ocp != nullptr);
  CHECK(negative_ocp->stoichiometry == std::vector<double>{ 0.0, 0.5, 1.0 });
  CHECK(negative_ocp->value == std::vector<double>{ 1.9793, 0.25, 0.1 });
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  CHECK(input.initial_soc == 0.75);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);

  const auto previous_size = parameters.size();
  CHECK(core::ParameterSet::fromBpxJson("{\"Header\":", parameters, diagnostic)
        == Status::Invalid_parameters);
  CHECK(parameters.size() == previous_size);
  CHECK_FALSE(diagnostic.empty());
  std::string old_version{ bpx_fixture };
  const auto version = old_version.find("1.0.0");
  REQUIRE(version != std::string::npos);
  old_version.replace(version, 5, "0.4.0");
  CHECK(core::ParameterSet::fromBpxJson(old_version, parameters, diagnostic)
        == Status::Invalid_parameters);
  CHECK(parameters.size() == previous_size);
}

TEST_CASE("ParameterSet public guards reject incomplete scalar and curve tables",
          "[core][parameters][validation][coverage]")
{
  core::ParameterSet parameters;
  CHECK(parameters.set("", 1.0, "test") == Status::Invalid_parameters);
  CHECK(parameters.set("x", 1.0, "") == Status::Invalid_parameters);
  CHECK(parameters.set(
          "curve",
          core::OCVCurve{ .stoichiometry = { 0.0, 1.0 },
                          .value = { 1.0 } },
          "test")
        == Status::Invalid_parameters);

  core::SpmFactoryInput sentinel;
  sentinel.design.capacity_Ah = 123.0;
  CHECK(parameters.toSpmInput(sentinel) == Status::Invalid_parameters);
  CHECK(sentinel.design.capacity_Ah == 123.0);

  core::ParameterSet complete;
  REQUIRE(core::ParameterSet::chen2020(complete) == Status::Success);

  auto missing_curves = copyWithout(
    complete,
    { "Negative electrode OCP [V]", "Positive electrode OCP [V]" });
  CHECK(missing_curves.toSpmInput(sentinel) == Status::Invalid_parameters);
  CHECK(sentinel.design.capacity_Ah == 123.0);

  auto missing_electrode = copyWithout(
    complete, { "Negative electrode thickness [m]" });
  CHECK(missing_electrode.toSpmInput(sentinel)
        == Status::Invalid_parameters);
  CHECK(sentinel.design.capacity_Ah == 123.0);

  auto missing_initial_state = copyWithout(
    complete,
    { "Initial state-of-charge",
      "Initial concentration in negative electrode [mol.m-3]" });
  CHECK(missing_initial_state.toSpmInput(sentinel)
        == Status::Invalid_parameters);
  CHECK(sentinel.design.capacity_Ah == 123.0);
}

TEST_CASE("BPX 1.x semantic functions and legacy headers are absorbed safely",
          "[core][parameters][BPX][functions]")
{
  std::string functional{ bpx_fixture };
  replaceRequired(functional, "\"BPX\": \"1.0.0\"", "\"BPX\": 1.0");
  replaceRequired(functional, "\"Model\": \"SPM\"", "\"Model\": \"DFN\"");
  replaceRequired(functional,
                  "\"Diffusivity [m2.s-1]\": 3.3e-14",
                  "\"Diffusivity [m2.s-1]\": \"3.3e-14\",\n"
                  "      \"Diffusivity activation energy [J.mol-1]\": 30000.0,\n"
                  "      \"Porosity\": 0.25398416831491194");
  replaceRequired(
    functional,
    "\"OCP [V]\": {\"x\": [0.0, 0.5, 1.0], \"y\": [1.9793, 0.25, 0.1]}",
    "\"OCP [V]\": \"0.1 + 0.2*x + exp(-2*x) + 0.01*cosh(x) - "
    "0.02*tanh(3*(x-0.5)) + 0.001*x**2\"");
  replaceRequired(
    functional,
    "\"OCP [V]\": {\"x\": [0.0, 0.5, 1.0], \"y\": [4.5, 3.8, 3.0]}",
    "\"OCP [V]\": \"4.5 - 1.5*x\"");

  core::ParameterSet parameters;
  std::string diagnostic;
  INFO(diagnostic);
  REQUIRE(core::ParameterSet::fromBpxJson(functional, parameters, diagnostic)
          == Status::Success);
  CHECK(diagnostic.empty());
  CHECK(*parameters.findScalar("Negative electrode porosity")
        == 0.25398416831491194);
  CHECK(*parameters.findScalar(
          "Negative particle diffusivity activation energy [J.mol-1]")
        == 30000.0);
  const auto *curve = parameters.findCurve("Negative electrode OCP [V]");
  REQUIRE(curve != nullptr);
  const double x = 0.37;
  const double expected = 0.1 + 0.2 * x + std::exp(-2.0 * x)
                          + 0.01 * std::cosh(x)
                          - 0.02 * std::tanh(3.0 * (x - 0.5))
                          + 0.001 * x * x;
  CHECK(std::abs(interpolate(*curve, x) - expected) <= 1e-6);

  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  CHECK(input.design.electrode[core::domain_index(core::Domain::neg)]
          .active_material.D_s.activation_energy
        == 30000.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);

  const auto original_size = parameters.size();
  std::string malicious{ functional };
  replaceRequired(malicious,
                  "\"OCP [V]\": \"4.5 - 1.5*x\"",
                  "\"OCP [V]\": \"__import__(x)\"");
  CHECK(core::ParameterSet::fromBpxJson(malicious, parameters, diagnostic)
        == Status::Invalid_parameters);
  CHECK(parameters.size() == original_size);
  CHECK(diagnostic.find("BPX expression") != std::string::npos);

  std::string variable_diffusivity{ functional };
  replaceRequired(variable_diffusivity,
                  "\"Diffusivity [m2.s-1]\": \"3.3e-14\"",
                  "\"Diffusivity [m2.s-1]\": \"3.3e-14*(1+x)\"");
  CHECK(core::ParameterSet::fromBpxJson(variable_diffusivity,
                                        parameters,
                                        diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.find("state-dependent BPX diffusivity")
        != std::string::npos);
  CHECK(parameters.size() == original_size);
}

TEST_CASE("BPX JSON rejects hostile grammar and resource amplification atomically",
          "[core][parameters][BPX][parser][P9]")
{
  core::ParameterSet parameters;
  std::string diagnostic;
  REQUIRE(core::ParameterSet::fromBpxJson(bpx_fixture, parameters, diagnostic)
          == Status::Success);
  const auto original = parameters.describe();

  auto require_atomic_rejection = [&](const std::string &source) {
    core::ParameterSet target = parameters;
    diagnostic.clear();
    CHECK(core::ParameterSet::fromBpxJson(source, target, diagnostic)
          == Status::Invalid_parameters);
    CHECK_FALSE(diagnostic.empty());
    const auto retained = target.describe();
    CHECK(retained.size() == original.size());
    const auto common = std::min(retained.size(), original.size());
    for (std::size_t i = 0; i < common; ++i) {
      CHECK(retained[i].name == original[i].name);
      CHECK(retained[i].value == original[i].value);
      CHECK(retained[i].provenance == original[i].provenance);
    }
  };

  std::string leading_zero{ bpx_fixture };
  replaceRequired(leading_zero,
                  "\"Nominal cell capacity [A.h]\": 5.0",
                  "\"Nominal cell capacity [A.h]\": 05.0");
  require_atomic_rejection(leading_zero);

  std::string invalid_utf8{ bpx_fixture };
  replaceRequired(invalid_utf8,
                  "LGM50 ",
                  std::string{ "LGM50 \xC0\xAF", 8 });
  require_atomic_rejection(invalid_utf8);

  for (const char invalid_space : { '\v', '\f' }) {
    std::string invalid_whitespace{ bpx_fixture };
    const auto newline = invalid_whitespace.find('\n');
    REQUIRE(newline != std::string::npos);
    invalid_whitespace[newline] = invalid_space;
    require_atomic_rejection(invalid_whitespace);
  }

  std::string invalid_optional{ bpx_fixture };
  replaceRequired(invalid_optional,
                  "\"Initial state-of-charge\": 0.75",
                  "\"Initial state-of-charge\": \"unknown\"");
  require_atomic_rejection(invalid_optional);

  std::string invalid_optional_cell{ bpx_fixture };
  replaceRequired(invalid_optional_cell,
                  "\"External surface area [m2]\": 0.00531",
                  "\"External surface area [m2]\": \"unknown\"");
  require_atomic_rejection(invalid_optional_cell);

  std::string invalid_optional_porosity{ bpx_fixture };
  replaceRequired(
    invalid_optional_porosity,
    "\"Surface area per unit volume [m-1]\": 383959.0443686007,",
    "\"Surface area per unit volume [m-1]\": 383959.0443686007,\n"
    "      \"Porosity\": \"unknown\",");
  require_atomic_rejection(invalid_optional_porosity);

  std::string missing_cell_scalar{ bpx_fixture };
  replaceRequired(missing_cell_scalar,
                  "\"Electrode area [m2]\"",
                  "\"Missing electrode area [m2]\"");
  require_atomic_rejection(missing_cell_scalar);

  std::string missing_electrode_number{ bpx_fixture };
  replaceRequired(missing_electrode_number,
                  "\"Thickness [m]\"",
                  "\"Missing thickness [m]\"");
  require_atomic_rejection(missing_electrode_number);

  std::string missing_electrode_constant{ bpx_fixture };
  replaceRequired(missing_electrode_constant,
                  "\"Diffusivity [m2.s-1]\"",
                  "\"Missing diffusivity [m2.s-1]\"");
  require_atomic_rejection(missing_electrode_constant);

  std::string missing_surface_area{ bpx_fixture };
  replaceRequired(missing_surface_area,
                  "\"Surface area per unit volume [m-1]\"",
                  "\"Missing surface area per unit volume [m-1]\"");
  require_atomic_rejection(missing_surface_area);

  std::string invalid_reaction_activation{ bpx_fixture };
  replaceRequired(invalid_reaction_activation,
                  "\"Reaction rate constant activation energy [J.mol-1]\": 35000.0",
                  "\"Reaction rate constant activation energy [J.mol-1]\": \"bad\"");
  require_atomic_rejection(invalid_reaction_activation);

  std::string invalid_diffusion_activation{ bpx_fixture };
  replaceRequired(invalid_diffusion_activation,
                  "\"Diffusivity [m2.s-1]\": 3.3e-14,",
                  "\"Diffusivity [m2.s-1]\": 3.3e-14,\n"
                  "      \"Diffusivity activation energy [J.mol-1]\": \"bad\",");
  require_atomic_rejection(invalid_diffusion_activation);

  // Bound the wire representation before constructing an amplified JSON tree.
  std::string oversized{ bpx_fixture };
  oversized.append(4U * 1024U * 1024U, ' ');
  require_atomic_rejection(oversized);

  // A compact ignored array must still consume the global parsed-value budget.
  std::string amplified{ bpx_fixture };
  const auto insertion = amplified.find("\"Parameterisation\"");
  REQUIRE(insertion != std::string::npos);
  std::string ignored{ "\"Ignored\":[null" };
  for (std::size_t i = 1; i < 65'537; ++i)
    ignored += ",null";
  ignored += "],\n  ";
  amplified.insert(insertion, ignored);
  require_atomic_rejection(amplified);

  std::string overflow{ bpx_fixture };
  replaceRequired(overflow,
                  "\"Nominal cell capacity [A.h]\": 5.0",
                  "\"Nominal cell capacity [A.h]\": 1e309");
  require_atomic_rejection(overflow);

  // Each operand is finite JSON, but the derived a*R/3 active fraction is not.
  std::string derived_overflow{ bpx_fixture };
  replaceRequired(
    derived_overflow,
    "\"Surface area per unit volume [m-1]\": 383959.0443686007",
    "\"Surface area per unit volume [m-1]\": 1e308");
  replaceRequired(derived_overflow,
                  "\"Particle radius [m]\": 5.86e-6",
                  "\"Particle radius [m]\": 1e308");
  require_atomic_rejection(derived_overflow);

  std::string raw_utf8{ bpx_fixture };
  replaceRequired(raw_utf8,
                  "\\u03bc",
                  std::string{ "\xCE\xBC", 2 });
  core::ParameterSet raw_parameters;
  diagnostic.clear();
  REQUIRE(core::ParameterSet::fromBpxJson(
            raw_utf8, raw_parameters, diagnostic)
          == Status::Success);
  CHECK(diagnostic.empty());
  CHECK(raw_parameters.size() == parameters.size());
}

TEST_CASE("BPX file reads are bounded, complete, and atomic",
          "[core][parameters][BPX][file][P9]")
{
  const ScopedTemporaryDirectory temporary;
  const auto path = temporary.path() / "fixture.json";

  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.write(bpx_fixture.data(),
                 static_cast<std::streamsize>(bpx_fixture.size()));
    REQUIRE(output.good());
  }
  core::ParameterSet parameters;
  std::string diagnostic;
  REQUIRE(core::ParameterSet::fromBpxFile(path, parameters, diagnostic)
          == Status::Success);
  CHECK(diagnostic.empty());
  const auto original = parameters.describe();

  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.seekp(static_cast<std::streamoff>(4U * 1024U * 1024U));
    output.put('x');
    REQUIRE(output.good());
  }
  CHECK(core::ParameterSet::fromBpxFile(path, parameters, diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.find("4194304") != std::string::npos);
  const auto retained = parameters.describe();
  REQUIRE(retained.size() == original.size());
  for (std::size_t i = 0; i < original.size(); ++i) {
    CHECK(retained[i].name == original[i].name);
    CHECK(retained[i].value == original[i].value);
    CHECK(retained[i].provenance == original[i].provenance);
  }

  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  CHECK(core::ParameterSet::fromBpxFile(path, parameters, diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.find("open") != std::string::npos);
}
