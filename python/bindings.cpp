/**
 * @file bindings.cpp
 * @brief Small nanobind boundary for the PyBaMM-compatible Python facade.
 */

#include "core/Experiment.hpp"
#include "core/ExponentialModal.hpp"
#include "core/ParameterSet.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

using slide::core::ParameterSet;

[[noreturn]] void fail(std::string message)
{
  throw std::invalid_argument(std::move(message));
}

ParameterSet loadParameters(const std::string &source)
{
  ParameterSet parameters;
  if (source == "Chen2020") {
    if (ParameterSet::chen2020(parameters) != slide::Status::Success)
      fail("failed to construct Chen2020 parameters");
    return parameters;
  }
  std::string diagnostic;
  if (ParameterSet::fromBpxFile(std::filesystem::path{ source }, parameters, diagnostic)
      != slide::Status::Success)
    fail(diagnostic.empty() ? "failed to read BPX parameter file" : diagnostic);
  return parameters;
}

slide::core::SpmModelOptions modelOptions(
  const std::map<std::string, int> &values)
{
  const auto get = [&values](std::string_view name, int fallback = 0) {
    const auto found = values.find(std::string{ name });
    return found == values.end() ? fallback : found->second;
  };
  slide::core::SpmModelOptions options;
  options.nch = get("nch", 8);
  options.thermal = get("thermal") != 0;
  options.sei_model_mask = static_cast<std::uint8_t>(get("sei_mask"));
  options.sei_porosity = get("sei_porosity") != 0;
  options.surface_crack_model_mask =
    static_cast<std::uint8_t>(get("crack_mask"));
  options.surface_crack_diffusivity = get("crack_diffusivity") != 0;
  options.lam_model_mask = static_cast<std::uint8_t>(get("lam_mask"));
  options.lithium_plating = get("plating") != 0;
  return options;
}

nb::dict describeParameters(const std::string &source)
{
  const auto parameters = loadParameters(source);
  nb::dict result;
  for (const auto &description : parameters.describe()) {
    if (const auto *scalar = std::get_if<double>(&description.value)) {
      result[nb::str(description.name.c_str())] = *scalar;
    } else {
      const auto &curve = std::get<slide::core::OCVCurve>(description.value);
      nb::dict encoded;
      encoded["x"] = nb::cast(curve.stoichiometry);
      encoded["y"] = nb::cast(curve.value);
      result[nb::str(description.name.c_str())] = std::move(encoded);
    }
  }
  return result;
}

nb::dict solveExperiment(
  const std::string &source,
  const std::map<std::string, double> &overrides,
  const std::map<std::string, int> &option_values,
  const std::vector<std::string> &steps,
  double sample_step)
{
  auto parameters = loadParameters(source);
  for (const auto &[name, value] : overrides)
    if (parameters.set(name, value, "Python override")
        != slide::Status::Success)
      fail("invalid parameter override: " + name);

  slide::core::SpmFactoryInput input;
  if (parameters.toSpmInput(input) != slide::Status::Success)
    fail("the parameter set is incomplete for the SPM registry");
  slide::core::SpmBatch batch;
  if (slide::core::buildSpmBatch(input, modelOptions(option_values), 1, batch)
      != slide::Status::Success)
    fail("the selected parameters/options could not build an SPM batch");

  slide::core::Experiment experiment;
  slide::core::ParseDiagnostic diagnostic;
  if (slide::core::Experiment::parse(steps, experiment, diagnostic)
      != slide::Status::Success)
    fail("experiment step " + std::to_string(diagnostic.step)
         + " at byte " + std::to_string(diagnostic.offset) + ": "
         + diagnostic.message);

  slide::core::CyclerV2 cycler;
  if (cycler.configure(batch, slide::core::CyclerIntegrator::exponential)
      != slide::Status::Success)
    fail("failed to configure the SPM experiment runner");
  slide::core::ExperimentSolution solution;
  const auto status = cycler.run(experiment, sample_step, solution);
  if (status != slide::Status::Success && solution.time.empty())
    fail("experiment failed before producing a sample");

  static constexpr const char *reasons[]{ "event", "limit", "error", "final time" };
  nb::dict result;
  result["time"] = nb::cast(std::move(solution.time));
  result["voltage"] = nb::cast(std::move(solution.voltage));
  result["current"] = nb::cast(std::move(solution.current));
  result["termination"] =
    reasons[static_cast<unsigned>(solution.reason)];
  result["segment"] = solution.segment;
  result["status"] = static_cast<int>(solution.status);
  return result;
}

enum class VariationTarget : unsigned char {
  initial_soc,
  temperature,
  negative_diffusivity,
  positive_diffusivity,
  negative_fraction,
  positive_fraction,
  negative_thickness,
  positive_thickness,
  contact_resistance
};

struct ValidatedVariation
{
  VariationTarget target{};
  const std::vector<double> *values{};
};

std::vector<ValidatedVariation> validateVariations(
  const std::map<std::string, std::vector<double>> &variations,
  int lanes)
{
  std::vector<ValidatedVariation> validated;
  validated.reserve(variations.size());
  for (const auto &[raw_name, values] : variations) {
    const std::string name = ParameterSet::canonicalName(raw_name);
    VariationTarget target;
    bool strictly_positive{};
    bool unit_interval{};
    if (name == "Initial state-of-charge") {
      target = VariationTarget::initial_soc;
      unit_interval = true;
    } else if (name == "Initial temperature [K]") {
      target = VariationTarget::temperature;
      strictly_positive = true;
    } else if (name == "Negative particle diffusivity [m2.s-1]") {
      target = VariationTarget::negative_diffusivity;
      strictly_positive = true;
    } else if (name == "Positive particle diffusivity [m2.s-1]") {
      target = VariationTarget::positive_diffusivity;
      strictly_positive = true;
    } else if (name == "Negative electrode active material volume fraction") {
      target = VariationTarget::negative_fraction;
      unit_interval = true;
    } else if (name == "Positive electrode active material volume fraction") {
      target = VariationTarget::positive_fraction;
      unit_interval = true;
    } else if (name == "Negative electrode thickness [m]") {
      target = VariationTarget::negative_thickness;
      strictly_positive = true;
    } else if (name == "Positive electrode thickness [m]") {
      target = VariationTarget::positive_thickness;
      strictly_positive = true;
    } else if (name == "Contact resistance [Ohm]") {
      target = VariationTarget::contact_resistance;
    } else {
      fail("varied() is not state-backed for parameter: " + raw_name);
    }
    if (static_cast<int>(values.size()) != lanes)
      fail("all varied() parameters must have the same lane count");
    for (const double value : values) {
      if (!slide::core::is_finite(value) || (strictly_positive && !(value > 0.0))
          || (unit_interval && !(value >= 0.0 && value <= 1.0))
          || ((target == VariationTarget::negative_fraction
               || target == VariationTarget::positive_fraction)
              && !(value > 0.0))
          || (target == VariationTarget::contact_resistance && value < 0.0))
        fail("invalid lane value for parameter: " + raw_name);
    }
    validated.push_back({ target, &values });
  }
  return validated;
}

void applyVariations(slide::core::SpmBatch &batch,
                     const slide::core::SpmFactoryInput &input,
                     const std::vector<ValidatedVariation> &variations)
{
  auto &state = batch.state();
  const auto &layout = batch.layout().spm;
  using slide::core::Domain;
  using slide::core::domain_index;
  for (const auto &variation : variations) {
    const auto &values = *variation.values;
    for (int lane = 0; lane < batch.n_lanes(); ++lane) {
      const double value = values[static_cast<std::size_t>(lane)];
      switch (variation.target) {
      case VariationTarget::initial_soc:
        for (const Domain domain : slide::core::domains) {
          const auto d = domain_index(domain);
          const auto &material = input.design.electrode[d].active_material;
          const double initial_stoichiometry = material.x_0
                                               + input.initial_soc
                                                   * (material.x_100 - material.x_0);
          const double varied_stoichiometry = material.x_0
                                              + value
                                                  * (material.x_100 - material.x_0);
          const double scale = varied_stoichiometry / initial_stoichiometry;
          for (int mode = 0; mode < layout.z[d].rows; ++mode)
            state.at(layout.z[d], mode, lane) *= scale;
        }
        break;
      case VariationTarget::temperature:
        state.at(layout.temperature, 0, lane) = value;
        break;
      case VariationTarget::negative_diffusivity:
        state.at(layout.diffusion_coefficient[domain_index(Domain::neg)], 0, lane) = value;
        break;
      case VariationTarget::positive_diffusivity:
        state.at(layout.diffusion_coefficient[domain_index(Domain::pos)], 0, lane) = value;
        break;
      case VariationTarget::negative_fraction:
      case VariationTarget::positive_fraction: {
        const Domain domain = variation.target == VariationTarget::negative_fraction
                                ? Domain::neg
                                : Domain::pos;
        const auto d = domain_index(domain);
        state.at(layout.active_fraction[d], 0, lane) = value;
        state.at(layout.specific_surface_area[d], 0, lane) =
          3.0 * value / input.design.electrode[d].particle_radius;
        break;
      }
      case VariationTarget::negative_thickness:
        state.at(layout.electrode_thickness[domain_index(Domain::neg)], 0, lane) = value;
        break;
      case VariationTarget::positive_thickness:
        state.at(layout.electrode_thickness[domain_index(Domain::pos)], 0, lane) = value;
        break;
      case VariationTarget::contact_resistance:
        state.at(layout.current_collector_resistance, 0, lane) =
          value * input.design.electrode_area;
        break;
      }
    }
  }
}

nb::dict solveEnsemble(
  const std::string &source,
  const std::map<std::string, double> &overrides,
  const std::map<std::string, std::vector<double>> &variations,
  const std::map<std::string, int> &option_values,
  const std::vector<std::string> &steps,
  double sample_step)
{
  if (variations.empty())
    fail("an ensemble requires at least one varied() parameter");
  const int lanes = static_cast<int>(variations.begin()->second.size());
  const auto validated = validateVariations(variations, lanes);

  auto parameters = loadParameters(source);
  for (const auto &[name, value] : overrides)
    if (parameters.set(name, value, "Python override")
        != slide::Status::Success)
      fail("invalid parameter override: " + name);
  slide::core::SpmFactoryInput input;
  if (parameters.toSpmInput(input) != slide::Status::Success)
    fail("the parameter set is incomplete for the SPM registry");

  slide::core::Experiment experiment;
  slide::core::ParseDiagnostic diagnostic;
  if (slide::core::Experiment::parse(steps, experiment, diagnostic)
        != slide::Status::Success
      || experiment.segments.size() != 1)
    fail("ensemble solve requires exactly one valid experiment segment");
  const auto &segment = experiment.segments.front();
  if (segment.mode != slide::core::ControlMode::current
      || !(segment.duration > 0.0) || segment.voltage_limit != 0.0
      || segment.current_cutoff != 0.0)
    fail("ensemble solve currently requires a fixed-duration CC segment without an event");
  if (!slide::core::is_finite(sample_step) || !(sample_step > 0.0))
    fail("sample_step must be finite and positive");

  slide::core::SpmBatch batch;
  if (slide::core::buildSpmBatch(input, modelOptions(option_values), lanes, batch)
      != slide::Status::Success)
    fail("the selected parameters/options could not build an ensemble batch");
  applyVariations(batch, input, validated);

  const double magnitude = segment.value_is_c_rate
                             ? segment.value * batch.capacity_Ah()
                             : segment.value;
  const double current = static_cast<int>(segment.direction) * magnitude;
  std::vector<double> current_density(static_cast<std::size_t>(lanes),
                                      current / batch.electrode_area());
  const std::size_t number_of_steps = static_cast<std::size_t>(
    std::ceil(segment.duration / sample_step));
  if (number_of_steps > std::numeric_limits<std::size_t>::max()
                          / static_cast<std::size_t>(lanes))
    fail("ensemble result shape is too large");
  std::vector<double> time(number_of_steps + 1);
  std::vector<double> voltage((number_of_steps + 1)
                              * static_cast<std::size_t>(lanes));
  std::vector<double> currents(voltage.size(), current);

  slide::core::StepCtx initial_ctx{ .time = 0.0,
                                    .dt = 0.0,
                                    .i_app = current_density };
  if (batch.terminalVoltage(initial_ctx,
                            std::span<double>{ voltage }.first(lanes))
      != slide::Status::Success)
    fail("ensemble initial observation failed");
  slide::core::ExponentialModal stepper;
  if (stepper.configure(batch) != slide::Status::Success)
    fail("ensemble stepper configuration failed");
  double simulation_time{};
  for (std::size_t step = 0; step < number_of_steps; ++step) {
    const double dt = std::min(sample_step, segment.duration - simulation_time);
    if (stepper.step(batch, current_density, simulation_time, dt)
        != slide::Status::Success)
      fail("ensemble integration failed at sample " + std::to_string(step));
    simulation_time += dt;
    time[step + 1] = simulation_time;
    std::copy(stepper.terminalVoltage().begin(),
              stepper.terminalVoltage().end(),
              voltage.begin()
                + static_cast<std::ptrdiff_t>((step + 1)
                                              * static_cast<std::size_t>(lanes)));
  }

  nb::dict result;
  result["time"] = nb::cast(std::move(time));
  result["voltage"] = nb::cast(std::move(voltage));
  result["current"] = nb::cast(std::move(currents));
  result["n_lanes"] = lanes;
  result["termination"] = "final time";
  result["segment"] = 0;
  result["status"] = 0;
  return result;
}

} // namespace

NB_MODULE(_slide_core, module)
{
  module.doc() = "Compiled SLIDE v4 SPM boundary";
  module.def("parameter_values", &describeParameters, "source"_a = "Chen2020");
  module.def("solve_experiment", &solveExperiment, "source"_a, "overrides"_a, "options"_a, "steps"_a, "sample_step"_a);
  module.def("solve_ensemble", &solveEnsemble, "source"_a, "overrides"_a, "variations"_a, "options"_a, "steps"_a, "sample_step"_a);
}
