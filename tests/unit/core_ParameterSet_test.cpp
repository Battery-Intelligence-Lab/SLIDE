/**
 * @file core_ParameterSet_test.cpp
 * @brief Phase-7 Chen2020/BPX parameter fidelity gates.
 */

#include "../../src/core/ParameterSet.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
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

} // namespace

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
