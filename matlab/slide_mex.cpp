/**
 * @file slide_mex.cpp
 * @brief One stateless dispatcher for the MATLAB +slide facade.
 */

#include "core/CudaSpmBatch.hpp"
#include "core/Experiment.hpp"
#include "core/ExponentialModal.hpp"
#include "core/ParameterSet.hpp"

#include "mex.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace {

namespace core = slide::core;

struct MexFailure : std::runtime_error
{
  MexFailure(std::string identifier, std::string message)
    : std::runtime_error(std::move(message)), identifier(std::move(identifier))
  {}
  std::string identifier;
};

[[noreturn]] void fail(std::string identifier, std::string message)
{
  throw MexFailure(std::move(identifier), std::move(message));
}

std::string text(const mxArray *value,
                 std::string_view identifier = "slide:InvalidParameter")
{
  if (value == nullptr || !mxIsChar(value))
    fail(std::string{ identifier }, "expected a character vector");
  std::vector<char> buffer(mxGetNumberOfElements(value) + 1);
  if (mxGetString(value, buffer.data(), buffer.size()) != 0)
    fail(std::string{ identifier }, "character conversion failed");
  return buffer.data();
}

double scalar(const mxArray *value,
              std::string_view identifier = "slide:InvalidParameter")
{
  if (value == nullptr || !mxIsDouble(value) || mxIsComplex(value)
      || mxGetNumberOfElements(value) != 1)
    fail(std::string{ identifier }, "expected one real double scalar");
  const double result = mxGetScalar(value);
  if (!core::is_finite(result))
    fail(std::string{ identifier }, "scalar must be finite");
  return result;
}

std::vector<std::string> strings(
  const mxArray *value,
  std::string_view identifier = "slide:InvalidParameter")
{
  if (value == nullptr || !mxIsCell(value))
    fail(std::string{ identifier }, "expected a cell array of character vectors");
  std::vector<std::string> result;
  result.reserve(mxGetNumberOfElements(value));
  for (std::size_t i = 0; i < mxGetNumberOfElements(value); ++i)
    result.push_back(text(mxGetCell(value, i), identifier));
  return result;
}

std::vector<double> doubles(
  const mxArray *value,
  std::string_view identifier = "slide:InvalidParameter")
{
  if (value == nullptr || !mxIsDouble(value) || mxIsComplex(value))
    fail(std::string{ identifier }, "expected a real double array");
  const auto count = mxGetNumberOfElements(value);
  const double *data = mxGetDoubles(value);
  std::vector<double> result(data, data + count);
  for (const double item : result)
    if (!core::is_finite(item))
      fail(std::string{ identifier }, "array values must be finite");
  return result;
}

mxArray *makeVector(std::span<const double> values)
{
  mxArray *result = mxCreateDoubleMatrix(values.size(), 1, mxREAL);
  std::copy(values.begin(), values.end(), mxGetDoubles(result));
  return result;
}

mxArray *makeIndexVector(std::span<const std::size_t> values)
{
  mxArray *result = mxCreateDoubleMatrix(values.size(), 1, mxREAL);
  auto *data = mxGetDoubles(result);
  for (std::size_t i = 0; i < values.size(); ++i)
    data[i] = static_cast<double>(values[i] + 1); // MATLAB segment indexing
  return result;
}

mxArray *makeSamples(std::span<const double> sample_major,
                     std::size_t samples,
                     std::size_t lanes)
{
  if (sample_major.size() != samples * lanes)
    fail("slide:Internal", "invalid sample-major result shape");
  mxArray *result = mxCreateDoubleMatrix(samples, lanes, mxREAL);
  auto *output = mxGetDoubles(result);
  for (std::size_t sample = 0; sample < samples; ++sample)
    for (std::size_t lane = 0; lane < lanes; ++lane)
      output[lane * samples + sample] = sample_major[sample * lanes + lane];
  return result;
}

mxArray *makeStringCell(const std::vector<std::string> &values)
{
  mxArray *result = mxCreateCellMatrix(values.size(), 1);
  for (std::size_t i = 0; i < values.size(); ++i)
    mxSetCell(result, i, mxCreateString(values[i].c_str()));
  return result;
}

core::ParameterSet loadParameters(const std::string &source)
{
  core::ParameterSet parameters;
  if (source == "Chen2020") {
    if (core::ParameterSet::chen2020(parameters) != slide::Status::Success)
      fail("slide:ParameterValues", "failed to construct Chen2020 parameters");
    return parameters;
  }
  std::string diagnostic;
  if (core::ParameterSet::fromBpxFile(source, parameters, diagnostic)
      != slide::Status::Success)
    fail("slide:BPX",
         diagnostic.empty() ? "failed to load BPX parameter file" : diagnostic);
  return parameters;
}

std::map<std::string, double> namedScalars(const mxArray *names,
                                           const mxArray *values,
                                           std::string_view identifier)
{
  const auto decoded_names = strings(names, identifier);
  const auto decoded_values = doubles(values, identifier);
  if (decoded_names.size() != decoded_values.size())
    fail(std::string{ identifier }, "name/value lengths differ");
  std::map<std::string, double> result;
  for (std::size_t i = 0; i < decoded_names.size(); ++i)
    if (!result.emplace(decoded_names[i], decoded_values[i]).second)
      fail(std::string{ identifier }, "duplicate name: " + decoded_names[i]);
  return result;
}

core::SpmModelOptions modelOptions(const std::map<std::string, double> &values)
{
  const auto integer = [&](std::string_view name, int fallback = 0) {
    const auto found = values.find(std::string{ name });
    if (found == values.end())
      return fallback;
    if (std::floor(found->second) != found->second
        || found->second < 0.0
        || found->second > static_cast<double>(std::numeric_limits<int>::max()))
      fail("slide:Options", "option must be a non-negative integer: "
                              + std::string{ name });
    return static_cast<int>(found->second);
  };
  core::SpmModelOptions options;
  options.nch = integer("nch", 8);
  options.thermal = integer("thermal") != 0;
  options.sei_model_mask = static_cast<std::uint8_t>(integer("sei_mask"));
  options.sei_porosity = integer("sei_porosity") != 0;
  options.surface_crack_model_mask = static_cast<std::uint8_t>(
    integer("crack_mask"));
  options.surface_crack_diffusivity = integer("crack_diffusivity") != 0;
  options.lam_model_mask = static_cast<std::uint8_t>(integer("lam_mask"));
  options.lithium_plating = integer("plating") != 0;
  return options;
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

struct Variation
{
  VariationTarget target{};
  std::vector<double> values{};
};

std::vector<Variation> decodeVariations(const mxArray *names,
                                        const mxArray *value_cells,
                                        int &lanes)
{
  const auto decoded_names = strings(names, "slide:Variation");
  if (!mxIsCell(value_cells)
      || mxGetNumberOfElements(value_cells) != decoded_names.size())
    fail("slide:Variation", "variation names/values differ");
  std::vector<Variation> result;
  lanes = decoded_names.empty() ? 1 : -1;
  for (std::size_t i = 0; i < decoded_names.size(); ++i) {
    const std::string name = core::ParameterSet::canonicalName(decoded_names[i]);
    VariationTarget target;
    bool positive{};
    bool unit_interval{};
    if (name == "Initial state-of-charge") {
      target = VariationTarget::initial_soc;
      unit_interval = true;
    } else if (name == "Initial temperature [K]") {
      target = VariationTarget::temperature;
      positive = true;
    } else if (name == "Negative particle diffusivity [m2.s-1]") {
      target = VariationTarget::negative_diffusivity;
      positive = true;
    } else if (name == "Positive particle diffusivity [m2.s-1]") {
      target = VariationTarget::positive_diffusivity;
      positive = true;
    } else if (name == "Negative electrode active material volume fraction") {
      target = VariationTarget::negative_fraction;
      positive = true;
      unit_interval = true;
    } else if (name == "Positive electrode active material volume fraction") {
      target = VariationTarget::positive_fraction;
      positive = true;
      unit_interval = true;
    } else if (name == "Negative electrode thickness [m]") {
      target = VariationTarget::negative_thickness;
      positive = true;
    } else if (name == "Positive electrode thickness [m]") {
      target = VariationTarget::positive_thickness;
      positive = true;
    } else if (name == "Contact resistance [Ohm]") {
      target = VariationTarget::contact_resistance;
    } else {
      fail("slide:Variation", "parameter is not state-backed: "
                                + decoded_names[i]);
    }
    auto values = doubles(mxGetCell(value_cells, i), "slide:Variation");
    if (values.empty()
        || values.size() > static_cast<std::size_t>(
             std::numeric_limits<int>::max()))
      fail("slide:Variation", "variation vectors must be non-empty");
    if (lanes < 0)
      lanes = static_cast<int>(values.size());
    if (static_cast<int>(values.size()) != lanes)
      fail("slide:Variation", "all variation vectors need the same lane count");
    for (const double value : values)
      if ((positive && !(value > 0.0))
          || (unit_interval && !(value >= 0.0 && value <= 1.0))
          || (target == VariationTarget::contact_resistance && value < 0.0))
        fail("slide:Variation", "invalid varied value for " + decoded_names[i]);
    result.push_back({ target, std::move(values) });
  }
  return result;
}

void applyVariations(core::SpmBatch &batch,
                     const core::SpmFactoryInput &input,
                     const std::vector<Variation> &variations)
{
  auto &state = batch.state();
  const auto &layout = batch.layout().spm;
  for (const auto &variation : variations)
    for (int lane = 0; lane < batch.n_lanes(); ++lane) {
      const double value = variation.values[static_cast<std::size_t>(lane)];
      switch (variation.target) {
      case VariationTarget::initial_soc:
        for (const core::Domain domain : core::domains) {
          const auto d = core::domain_index(domain);
          const auto &material = input.design.electrode[d].active_material;
          const double base = material.x_0
                              + input.initial_soc
                                  * (material.x_100 - material.x_0);
          const double varied = material.x_0
                                + value * (material.x_100 - material.x_0);
          for (int mode = 0; mode < layout.z[d].rows; ++mode)
            state.at(layout.z[d], mode, lane) *= varied / base;
        }
        break;
      case VariationTarget::temperature:
        state.at(layout.temperature, 0, lane) = value;
        break;
      case VariationTarget::negative_diffusivity:
        state.at(layout.diffusion_coefficient[core::domain_index(
                   core::Domain::neg)], 0, lane) = value;
        break;
      case VariationTarget::positive_diffusivity:
        state.at(layout.diffusion_coefficient[core::domain_index(
                   core::Domain::pos)], 0, lane) = value;
        break;
      case VariationTarget::negative_fraction:
      case VariationTarget::positive_fraction: {
        const auto domain = variation.target == VariationTarget::negative_fraction
                              ? core::Domain::neg
                              : core::Domain::pos;
        const auto d = core::domain_index(domain);
        state.at(layout.active_fraction[d], 0, lane) = value;
        state.at(layout.specific_surface_area[d], 0, lane) =
          3.0 * value / input.design.electrode[d].particle_radius;
        break;
      }
      case VariationTarget::negative_thickness:
        state.at(layout.electrode_thickness[core::domain_index(
                   core::Domain::neg)], 0, lane) = value;
        break;
      case VariationTarget::positive_thickness:
        state.at(layout.electrode_thickness[core::domain_index(
                   core::Domain::pos)], 0, lane) = value;
        break;
      case VariationTarget::contact_resistance:
        state.at(layout.current_collector_resistance, 0, lane) =
          value * input.design.electrode_area;
        break;
      }
    }
}

struct Result
{
  std::vector<double> time{};
  std::vector<double> voltage{};
  std::vector<double> current{};
  std::vector<std::size_t> segment{};
  std::string termination{ "final time" };
  std::string termination_detail{ "final time" };
  std::size_t terminal_segment{};
  int status{};
  int lanes{ 1 };
};

Result solveSingle(core::SpmFactoryInput input,
                   const core::SpmModelOptions &options,
                   const std::vector<std::string> &steps,
                   double sample_step)
{
  core::SpmBatch batch;
  if (core::buildSpmBatch(input, options, 1, batch) != slide::Status::Success)
    fail("slide:Options", "parameters/options could not build an SPM batch");
  core::Experiment experiment;
  core::ParseDiagnostic diagnostic;
  if (core::Experiment::parse(steps, experiment, diagnostic)
      != slide::Status::Success)
    fail("slide:Experiment",
         "step " + std::to_string(diagnostic.step + 1)
           + " at byte " + std::to_string(diagnostic.offset + 1)
           + ": " + diagnostic.message);
  core::CyclerV2 cycler;
  if (cycler.configure(batch, core::CyclerIntegrator::exponential)
      != slide::Status::Success)
    fail("slide:Solve", "failed to configure the experiment runner");
  core::ExperimentSolution native;
  const auto status = cycler.run(experiment, sample_step, native);
  if (status != slide::Status::Success && native.time.empty())
    fail("slide:Solve", "experiment failed before producing a sample");
  static constexpr const char *reasons[]{ "event", "limit", "error", "final time" };
  Result result;
  result.time = std::move(native.time);
  result.voltage = std::move(native.voltage);
  result.current = std::move(native.current);
  result.segment = std::move(native.sample_segment);
  result.termination = reasons[static_cast<unsigned>(native.reason)];
  result.termination_detail = native.termination_name.empty()
                                ? result.termination
                                : native.termination_name;
  result.terminal_segment = native.segment;
  result.status = static_cast<int>(native.status);
  return result;
}

Result solveBatch(core::SpmFactoryInput input,
                  const core::SpmModelOptions &options,
                  const std::vector<std::string> &steps,
                  double sample_step,
                  const std::string &device,
                  const std::vector<Variation> &variations,
                  int lanes)
{
  core::Experiment experiment;
  core::ParseDiagnostic diagnostic;
  if (core::Experiment::parse(steps, experiment, diagnostic)
        != slide::Status::Success
      || experiment.segments.size() != 1)
    fail("slide:Experiment",
         "varied/CUDA solves require exactly one valid experiment segment");
  const auto &segment = experiment.segments.front();
  if (segment.mode != core::ControlMode::current || !(segment.duration > 0.0)
      || segment.voltage_limit != 0.0 || segment.current_cutoff != 0.0)
    fail("slide:Experiment",
         "varied/CUDA solves require fixed-duration CC without an event");

  core::SpmBatch cpu;
  core::CudaSpmBatch cuda;
  const auto build_status = device == "cuda"
                              ? cuda.build(input, options, lanes)
                              : core::buildSpmBatch(input, options, lanes, cpu);
  if (build_status != slide::Status::Success)
    fail(device == "cuda" ? "slide:Device" : "slide:Options",
         "parameters/options could not build the selected batch");
  auto &batch = device == "cuda" ? cuda.hostBatch() : cpu;
  applyVariations(batch, input, variations);
  const double magnitude = segment.value_is_c_rate
                             ? segment.value * batch.capacity_Ah()
                             : segment.value;
  const double current = static_cast<int>(segment.direction) * magnitude;
  std::vector<double> density(static_cast<std::size_t>(lanes),
                              current / batch.electrode_area());
  const auto number_of_steps = static_cast<std::size_t>(
    std::ceil(segment.duration / sample_step));
  Result result;
  result.lanes = lanes;
  result.time.resize(number_of_steps + 1);
  result.voltage.resize((number_of_steps + 1)
                        * static_cast<std::size_t>(lanes));
  result.current.assign(result.voltage.size(), current);
  result.segment.assign(number_of_steps + 1, 0);
  const core::StepCtx initial{ .time = 0.0, .dt = 0.0, .i_app = density };
  if (batch.terminalVoltage(initial,
                            std::span{ result.voltage }.first(
                              static_cast<std::size_t>(lanes)))
      != slide::Status::Success)
    fail("slide:Solve", "initial observation failed");
  core::ExponentialModal stepper;
  if (device == "cpu" && stepper.configure(batch) != slide::Status::Success)
    fail("slide:Solve", "CPU stepper configuration failed");
  if (device == "cuda" && cuda.uploadState() != slide::Status::Success)
    fail("slide:Device", "CUDA state upload failed");
  double time{};
  for (std::size_t step = 0; step < number_of_steps; ++step) {
    const double dt = std::min(sample_step, segment.duration - time);
    std::span<const double> sample_voltage;
    if (device == "cuda") {
      if (cuda.step(density, time, dt) != slide::Status::Success
          || cuda.synchronize() != slide::Status::Success)
        fail("slide:Solve", "CUDA integration failed at sample "
                              + std::to_string(step + 1));
      sample_voltage = cuda.terminalVoltage();
    } else {
      if (stepper.step(batch, density, time, dt) != slide::Status::Success)
        fail("slide:Solve", "CPU integration failed at sample "
                              + std::to_string(step + 1));
      sample_voltage = stepper.terminalVoltage();
    }
    time += dt;
    result.time[step + 1] = time;
    std::copy(sample_voltage.begin(), sample_voltage.end(),
              result.voltage.begin()
                + static_cast<std::ptrdiff_t>((step + 1)
                                              * static_cast<std::size_t>(lanes)));
  }
  return result;
}

mxArray *encodeResult(const Result &result)
{
  const char *fields[]{ "time", "voltage", "current",
                        "sample_segment", "termination",
                        "termination_detail", "segment",
                        "status", "n_lanes" };
  mxArray *output = mxCreateStructMatrix(1, 1, 9, fields);
  mxSetField(output, 0, "time", makeVector(result.time));
  mxSetField(output, 0, "voltage",
             makeSamples(result.voltage, result.time.size(), result.lanes));
  mxSetField(output, 0, "current",
             makeSamples(result.current, result.time.size(), result.lanes));
  mxSetField(output, 0, "sample_segment", makeIndexVector(result.segment));
  mxSetField(output, 0, "termination",
             mxCreateString(result.termination.c_str()));
  mxSetField(output, 0, "termination_detail",
             mxCreateString(result.termination_detail.c_str()));
  mxSetField(output, 0, "segment",
             mxCreateDoubleScalar(static_cast<double>(result.terminal_segment + 1)));
  mxSetField(output, 0, "status", mxCreateDoubleScalar(result.status));
  mxSetField(output, 0, "n_lanes", mxCreateDoubleScalar(result.lanes));
  return output;
}

void commandDevices(int nlhs, mxArray **plhs)
{
  if (nlhs > 1)
    fail("slide:InvalidParameter", "devices returns at most one output");
  std::vector<std::string> devices{ "cpu" };
  if (core::CudaSpmBatch::available())
    devices.emplace_back("cuda");
  if (nlhs == 1)
    plhs[0] = makeStringCell(devices);
}

void commandParameters(int nlhs, mxArray **plhs,
                       int nrhs, const mxArray **prhs)
{
  if (nrhs != 2 || nlhs != 2)
    fail("slide:InvalidParameter", "parameters needs a source and two outputs");
  const auto parameters = loadParameters(text(prhs[1], "slide:ParameterValues"));
  const auto descriptions = parameters.describe();
  mxArray *names = mxCreateCellMatrix(descriptions.size(), 1);
  mxArray *values = mxCreateCellMatrix(descriptions.size(), 1);
  for (std::size_t i = 0; i < descriptions.size(); ++i) {
    mxSetCell(names, i, mxCreateString(descriptions[i].name.c_str()));
    if (const auto *number = std::get_if<double>(&descriptions[i].value)) {
      mxSetCell(values, i, mxCreateDoubleScalar(*number));
    } else {
      const auto &curve = std::get<core::OCVCurve>(descriptions[i].value);
      mxArray *matrix = mxCreateDoubleMatrix(2, curve.value.size(), mxREAL);
      auto *data = mxGetDoubles(matrix);
      for (std::size_t knot = 0; knot < curve.value.size(); ++knot) {
        data[2 * knot] = curve.stoichiometry[knot];
        data[2 * knot + 1] = curve.value[knot];
      }
      mxSetCell(values, i, matrix);
    }
  }
  plhs[0] = names;
  plhs[1] = values;
}

void commandSolve(int nlhs, mxArray **plhs,
                  int nrhs, const mxArray **prhs)
{
  if (nrhs != 11 || nlhs != 1)
    fail("slide:InvalidParameter",
         "solve expects source, overrides, options, steps, period, device, and variations");
  auto parameters = loadParameters(text(prhs[1], "slide:ParameterValues"));
  const auto overrides = namedScalars(prhs[2], prhs[3],
                                      "slide:ParameterValues");
  for (const auto &[name, value] : overrides)
    if (parameters.set(name, value, "MATLAB override")
        != slide::Status::Success)
      fail("slide:ParameterValues", "invalid parameter override: " + name);
  core::SpmFactoryInput input;
  if (parameters.toSpmInput(input) != slide::Status::Success)
    fail("slide:ParameterValues", "parameter set is incomplete for SPM");
  const auto options = modelOptions(namedScalars(prhs[4], prhs[5],
                                                  "slide:Options"));
  const auto steps = strings(prhs[6], "slide:Experiment");
  const double sample_step = scalar(prhs[7], "slide:Experiment");
  if (!(sample_step > 0.0))
    fail("slide:Experiment", "sample period must be positive");
  const std::string device = text(prhs[8], "slide:Device");
  if (device != "cpu" && device != "cuda")
    fail("slide:Device", "device must be 'cpu' or 'cuda'");
  if (device == "cuda" && !core::CudaSpmBatch::available())
    fail("slide:Device", "CUDA is unavailable in this build/runtime");
  int lanes{};
  const auto variations = decodeVariations(prhs[9], prhs[10], lanes);
  const bool batch_path = device == "cuda" || !variations.empty();
  const Result result = batch_path
                          ? solveBatch(input, options, steps, sample_step,
                                       device, variations, lanes)
                          : solveSingle(input, options, steps, sample_step);
  plhs[0] = encodeResult(result);
}

} // namespace

#if defined(_WIN32)
#define SLIDE_MEX_EXPORT __declspec(dllexport)
#else
#define SLIDE_MEX_EXPORT __attribute__((visibility("default")))
#endif

extern "C" SLIDE_MEX_EXPORT void mexFunction(
  int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  try {
    if (nrhs < 1)
      fail("slide:InvalidParameter", "first argument must be a command");
    const std::string command = text(prhs[0]);
    if (command == "devices") {
      if (nrhs != 1)
        fail("slide:InvalidParameter", "devices takes no arguments");
      commandDevices(nlhs, plhs);
    } else if (command == "parameters") {
      commandParameters(nlhs, plhs, nrhs, prhs);
    } else if (command == "solve") {
      commandSolve(nlhs, plhs, nrhs, prhs);
    } else {
      fail("slide:InvalidCommand", "unknown dispatcher command: " + command);
    }
  } catch (const MexFailure &error) {
    mexErrMsgIdAndTxt(error.identifier.c_str(), "%s", error.what());
  } catch (const std::exception &error) {
    mexErrMsgIdAndTxt("slide:Internal", "%s", error.what());
  } catch (...) {
    mexErrMsgIdAndTxt("slide:Internal", "unknown native exception");
  }
}
