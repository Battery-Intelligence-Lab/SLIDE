/**
 * @file ForwardSensitivity.cpp
 * @brief Dual-number modal propagation and voltage observation.
 */

#include "ForwardSensitivity.hpp"

#include "Dual.hpp"
#include "CompiledCurve.hpp"
#include "ParameterSet.hpp"
#include "SpmScalarKernels.hpp"
#include "SpectralModel.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <new>

namespace slide::core {
namespace {

  static_assert(registered_spm_nch == std::array{ 5, 8, 12 },
                "forward-sensitivity dispatch must match the SPM registry");

  Dual seed(real_t value, SensitivityParameter active,
            SensitivityParameter candidate)
  {
    return { value, active == candidate ? 1.0 : 0.0 };
  }

  template <int NCH>
  class DualSpm
  {
  public:
    DualSpm(const SpmFactoryInput &input,
            SensitivityParameter active,
            real_t control_magnitude,
            bool control_is_c_rate,
            Direction direction)
      : input_{ input }, active_{ active }
    {
      const PerDomain<real_t> radius{
        input.design.electrode[domain_index(Domain::neg)].particle_radius,
        input.design.electrode[domain_index(Domain::pos)].particle_radius,
      };
      status_ = compileSpectralModel<NCH>(radius, spectral_);
      if (status_ != slide::Status::Success)
        return;
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const auto &curve = input.design.electrode[d].active_material.ocv;
        status_ = ocv_[d].build(curve.stoichiometry, curve.value);
        if (status_ != slide::Status::Success)
          return;
      }
      const Dual capacity = seed(input.design.capacity_Ah, active, SensitivityParameter::nominal_capacity);
      current_ = static_cast<real_t>(direction)
                 * (control_is_c_rate ? Dual{ control_magnitude } * capacity
                                      : Dual{ control_magnitude });
      current_density_ = current_ / input.design.electrode_area;
      initialiseState();
    }

    slide::Status status() const { return status_; }

    Dual observe() const
    {
      PerDomain<Dual> stoichiometry{};
      PerDomain<Dual> overpotential{};
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const auto &electrode = input_.design.electrode[d];
        const auto &material = electrode.active_material;
        const Dual diffusivity = effectiveDiffusivity(domain);
        const Dual flux = molarFlux(domain);
        Dual surface{};
        for (int mode = 0; mode < NCH; ++mode)
          surface += spectral_.C[d][0][static_cast<std::size_t>(mode)]
                     * z_[d][static_cast<std::size_t>(mode)];
        surface = spm_scalar::concentrationOutput(
          surface, spectral_.D[d][0], flux, diffusivity);
        stoichiometry[d] = spm_scalar::surfaceStoichiometry(
          surface, material.cs_max);
        if (!(stoichiometry[d].value > 0.0 && stoichiometry[d].value < 1.0))
          return { std::numeric_limits<real_t>::quiet_NaN(), 0.0 };

        const SensitivityParameter reaction_parameter =
          domain == Domain::neg
            ? SensitivityParameter::negative_reaction_rate
            : SensitivityParameter::positive_reaction_rate;
        const Dual reaction_ref = seed(material.k_ct.reference_value, active_, reaction_parameter);
        const real_t temperature = temperature_;
        const real_t arrhenius = spm_scalar::arrheniusFactor(
          material.k_ct.reference_temperature, temperature, 8.314);
        const Dual reaction_rate = spm_scalar::activatedValue(
          reaction_ref, material.k_ct.activation_energy, arrhenius);
        const Dual exchange_current = spm_scalar::exchangeCurrent(
          reaction_rate,
          1.0,
          96487.0,
          input_.design.electrolyte.concentration,
          surface,
          material.cs_max);
        const real_t specific_area = 3.0 * electrode.active_fraction
                                     / electrode.particle_radius;
        const Dual argument = spm_scalar::activationArgument(
          molar_flux_sign(domain),
          current_density_,
          specific_area,
          electrode.thickness,
          exchange_current);
        overpotential[d] = spm_scalar::activationOverpotential(
          temperature, 8.314, 1.0, 96487.0, argument);
      }

      const auto neg = domain_index(Domain::neg);
      const auto pos = domain_index(Domain::pos);
      const Dual open_circuit = spm_scalar::cellOpenCircuitVoltage(
        ocv_[neg].eval(stoichiometry[neg]),
        ocv_[pos].eval(stoichiometry[pos]),
        temperature_,
        temperature_,
        0.0);
      const auto &negative = input_.design.electrode[neg];
      const auto &positive = input_.design.electrode[pos];
      const real_t area_neg = spm_scalar::activeArea(
        3.0 * negative.active_fraction / negative.particle_radius,
        input_.design.electrode_area,
        negative.thickness);
      const real_t area_pos = spm_scalar::activeArea(
        3.0 * positive.active_fraction / positive.particle_radius,
        input_.design.electrode_area,
        positive.thickness);
      const Dual collector_resistance_area{
        input_.initial_current_collector_resistance,
        active_ == SensitivityParameter::contact_resistance
          ? input_.design.electrode_area
          : 0.0
      };
      const Dual resistance = spm_scalar::seriesResistance(
        input_.initial_sei_thickness,
        input_.sei_resistivity_area,
        input_.initial_specific_resistance[neg],
        input_.initial_specific_resistance[pos],
        collector_resistance_area,
        area_neg,
        area_pos,
        input_.design.electrode_area);
      return spm_scalar::terminalVoltage(open_circuit,
                                         overpotential[neg],
                                         overpotential[pos],
                                         resistance,
                                         current_);
    }

    slide::Status step(real_t dt)
    {
      // solveOne is the sole caller. Its validated schedule uses either a
      // positive sample_step or the positive final remainder; near-integral
      // ratios are snapped before the sample count is constructed so rounding
      // cannot manufacture a zero-length final step.
      assert(is_finite(dt) && dt > 0.0);
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const Dual diffusivity = effectiveDiffusivity(domain);
        const Dual flux = molarFlux(domain);
        for (int mode = 0; mode < NCH; ++mode) {
          auto &value = z_[d][static_cast<std::size_t>(mode)];
          SLIDE_SPM_ADVANCE_MODAL_ADL(
            value,
            diffusivity,
            spectral_.A[d][static_cast<std::size_t>(mode)],
            dt,
            spectral_.B[d][static_cast<std::size_t>(mode)],
            flux);
        }
      }
      const Dual voltage = observe();
      return is_finite(voltage.value) && is_finite(voltage.derivative)
               ? slide::Status::Success
               : slide::Status::Invalid_states;
    }

  private:
    Dual materialX0(Domain domain) const
    {
      const auto &material = input_.design.electrode[domain_index(domain)].active_material;
      const auto parameter = domain == Domain::neg
                               ? SensitivityParameter::negative_minimum_stoichiometry
                               : SensitivityParameter::positive_maximum_stoichiometry;
      return seed(material.x_0, active_, parameter);
    }

    Dual materialX100(Domain domain) const
    {
      const auto &material = input_.design.electrode[domain_index(domain)].active_material;
      const auto parameter = domain == Domain::neg
                               ? SensitivityParameter::negative_maximum_stoichiometry
                               : SensitivityParameter::positive_minimum_stoichiometry;
      return seed(material.x_100, active_, parameter);
    }

    void initialiseState()
    {
      temperature_ = input_.initial_temperature > 0.0
                       ? input_.initial_temperature
                       : input_.design.thermal.reference_temperature;
      for (const Domain domain : domains) {
        const auto d = domain_index(domain);
        const auto &electrode = input_.design.electrode[d];
        const auto &material = electrode.active_material;
        const Dual stoichiometry = materialX0(domain)
                                   + input_.initial_soc
                                       * (materialX100(domain) - materialX0(domain));
        const Dual concentration = stoichiometry * material.cs_max;
        const int zero = spectral_.zero_mode[d];
        Dual uniform_mode{};
        for (int node = 0; node < NCH; ++node) {
          const Dual radial_state = electrode.particle_radius * concentration
                                    * spectral_.x_inner[static_cast<std::size_t>(node)];
          uniform_mode += spectral_.state_transform[d][static_cast<std::size_t>(zero)]
                                                   [static_cast<std::size_t>(node)]
                          * radial_state;
        }
        z_[d][static_cast<std::size_t>(zero)] = uniform_mode;
      }
    }

    Dual effectiveDiffusivity(Domain domain) const
    {
      const auto &material = input_.design.electrode[domain_index(domain)].active_material;
      const auto parameter = domain == Domain::neg
                               ? SensitivityParameter::negative_diffusivity
                               : SensitivityParameter::positive_diffusivity;
      const Dual reference = seed(material.D_s.reference_value, active_, parameter);
      const real_t arrhenius = spm_scalar::arrheniusFactor(
        material.D_s.reference_temperature, temperature_, 8.314);
      return spm_scalar::activatedValue(
        reference, material.D_s.activation_energy, arrhenius);
    }

    Dual molarFlux(Domain domain) const
    {
      const auto &electrode = input_.design.electrode[domain_index(domain)];
      const real_t specific_area = 3.0 * electrode.active_fraction
                                   / electrode.particle_radius;
      const real_t denominator = spm_scalar::fluxDenominator(
        specific_area, 1.0, 96487.0, electrode.thickness);
      return spm_scalar::molarFlux(
        molar_flux_sign(domain), current_density_, denominator);
    }

    const SpmFactoryInput &input_;
    SensitivityParameter active_{};
    CompiledSpectralModel<NCH> spectral_{};
    PerDomain<IndexedPiecewiseLinear> ocv_{};
    PerDomain<std::array<Dual, NCH>> z_{};
    Dual current_{};
    Dual current_density_{};
    real_t temperature_{};
    slide::Status status_{ slide::Status::Success };
  };

  template <int NCH>
  slide::Status solveOne(const SpmFactoryInput &input,
                         real_t control_magnitude,
                         bool control_is_c_rate,
                         Direction direction,
                         real_t duration,
                         real_t sample_step,
                         SensitivityParameter parameter,
                         std::span<real_t>
                           time,
                         std::span<real_t>
                           voltage,
                         std::span<real_t>
                           derivative,
                         bool write_primal)
  {
    DualSpm<NCH> model{ input, parameter, control_magnitude, control_is_c_rate, direction };
    if (model.status() != slide::Status::Success)
      return model.status();
    Dual observation = model.observe();
    if (!(is_finite(observation.value) && is_finite(observation.derivative)))
      return slide::Status::Invalid_states;
    if (write_primal) {
      time[0] = 0.0;
      voltage[0] = observation.value;
    }
    derivative[0] = observation.derivative;
    real_t simulation_time{};
    for (std::size_t sample = 1; sample < time.size(); ++sample) {
      const bool final_sample = sample + 1 == time.size();
      const real_t dt = final_sample ? duration - simulation_time
                                     : sample_step;
      const auto status = model.step(dt);
      if (status != slide::Status::Success)
        return status;
      simulation_time = final_sample ? duration : simulation_time + dt;
      observation = model.observe();
      if (write_primal) {
        time[sample] = simulation_time;
        voltage[sample] = observation.value;
      }
      derivative[sample] = observation.derivative;
    }
    return slide::Status::Success;
  }

} // namespace

std::string_view sensitivityParameterName(SensitivityParameter parameter)
{
  switch (parameter) {
  case SensitivityParameter::negative_diffusivity:
    return "Negative particle diffusivity [m2.s-1]";
  case SensitivityParameter::positive_diffusivity:
    return "Positive particle diffusivity [m2.s-1]";
  case SensitivityParameter::negative_reaction_rate:
    return "Negative electrode reaction rate constant [mol.m-2.s-1]";
  case SensitivityParameter::positive_reaction_rate:
    return "Positive electrode reaction rate constant [mol.m-2.s-1]";
  case SensitivityParameter::contact_resistance:
    return "Contact resistance [Ohm]";
  case SensitivityParameter::nominal_capacity:
    return "Nominal cell capacity [A.h]";
  case SensitivityParameter::negative_minimum_stoichiometry:
    return "Negative electrode minimum stoichiometry";
  case SensitivityParameter::negative_maximum_stoichiometry:
    return "Negative electrode maximum stoichiometry";
  case SensitivityParameter::positive_minimum_stoichiometry:
    return "Positive electrode minimum stoichiometry";
  case SensitivityParameter::positive_maximum_stoichiometry:
    return "Positive electrode maximum stoichiometry";
  }
  return {};
}

slide::Status parseSensitivityParameter(std::string_view name,
                                        SensitivityParameter &parameter)
{
  const std::string canonical = ParameterSet::canonicalName(name);
  for (const auto candidate : supported_sensitivity_parameters)
    if (canonical == sensitivityParameterName(candidate)) {
      parameter = candidate;
      return slide::Status::Success;
    }
  return slide::Status::Invalid_parameters;
}

real_t sensitivityParameterValue(const SpmFactoryInput &input,
                                 SensitivityParameter parameter)
{
  const auto neg = domain_index(Domain::neg);
  const auto pos = domain_index(Domain::pos);
  switch (parameter) {
  case SensitivityParameter::negative_diffusivity:
    return input.design.electrode[neg].active_material.D_s.reference_value;
  case SensitivityParameter::positive_diffusivity:
    return input.design.electrode[pos].active_material.D_s.reference_value;
  case SensitivityParameter::negative_reaction_rate:
    return input.design.electrode[neg].active_material.k_ct.reference_value;
  case SensitivityParameter::positive_reaction_rate:
    return input.design.electrode[pos].active_material.k_ct.reference_value;
  case SensitivityParameter::contact_resistance:
    return input.initial_current_collector_resistance
           / input.design.electrode_area;
  case SensitivityParameter::nominal_capacity:
    return input.design.capacity_Ah;
  case SensitivityParameter::negative_minimum_stoichiometry:
    return input.design.electrode[neg].active_material.x_0;
  case SensitivityParameter::negative_maximum_stoichiometry:
    return input.design.electrode[neg].active_material.x_100;
  case SensitivityParameter::positive_minimum_stoichiometry:
    return input.design.electrode[pos].active_material.x_100;
  case SensitivityParameter::positive_maximum_stoichiometry:
    return input.design.electrode[pos].active_material.x_0;
  }
  return std::numeric_limits<real_t>::quiet_NaN();
}

slide::Status solveCcForwardSensitivities(
  const SpmFactoryInput &input,
  int nch,
  real_t control_magnitude,
  bool control_is_c_rate,
  Direction direction,
  real_t duration,
  real_t sample_step,
  std::span<const SensitivityParameter>
    parameters,
  ForwardSensitivitySolution &output)
try {
  if (parameters.empty() || !(direction == Direction::charge || direction == Direction::discharge)
      || !(is_finite(control_magnitude) && control_magnitude > 0.0)
      || !(is_finite(duration) && duration >= 0.0)
      || !(is_finite(sample_step) && sample_step > 0.0)
      || std::find(registered_spm_nch.begin(), registered_spm_nch.end(), nch)
           == registered_spm_nch.end())
    return slide::Status::Invalid_parameters;
  for (std::size_t i = 0; i < parameters.size(); ++i)
    if (std::find(supported_sensitivity_parameters.begin(),
                  supported_sensitivity_parameters.end(),
                  parameters[i])
          == supported_sensitivity_parameters.end()
        || std::find(parameters.begin(), parameters.begin() + i, parameters[i])
             != parameters.begin() + i)
      return slide::Status::Invalid_parameters;

  SpmModelOptions options;
  options.nch = nch;
  SpmBatch validation;
  if (buildSpmBatch(input, options, 1, validation) != slide::Status::Success)
    return slide::Status::Invalid_parameters;
  real_t step_ratio = duration / sample_step;
  const real_t nearest_ratio = std::round(step_ratio);
  const real_t snap_tolerance =
    16.0 * std::numeric_limits<real_t>::epsilon()
    * std::max(real_t{ 1.0 }, std::abs(step_ratio));
  // Never snap a positive sub-step duration down to zero. For positive
  // integral ratios, snapping prevents quotient rounding from manufacturing a
  // final zero-length step (for example, a rounded duration/7 interval).
  if (nearest_ratio >= 1.0
      && std::abs(step_ratio - nearest_ratio) <= snap_tolerance)
    step_ratio = nearest_ratio;
  const real_t raw_steps = duration == 0.0
                             ? 0.0
                             : std::max(real_t{ 1.0 }, std::ceil(step_ratio));
  if (!is_finite(raw_steps)
      || raw_steps
           >= static_cast<real_t>(std::numeric_limits<std::size_t>::max()))
    return slide::Status::Invalid_parameters;
  const auto sample_count = static_cast<std::size_t>(raw_steps) + 1;
  const auto real_vector_limit = std::vector<real_t>{}.max_size();
  if (sample_count > real_vector_limit
      || parameters.size() > std::vector<SensitivityParameter>{}.max_size()
      || sample_count > real_vector_limit / parameters.size())
    return slide::Status::Invalid_parameters;
  ForwardSensitivitySolution candidate;
  candidate.time.resize(sample_count);
  candidate.terminal_voltage.resize(sample_count);
  candidate.parameters.assign(parameters.begin(), parameters.end());
  candidate.derivative.resize(sample_count * parameters.size());
  std::vector<real_t> one_derivative(sample_count);

  for (std::size_t p = 0; p < parameters.size(); ++p) {
    slide::Status status;
    if (nch == 5)
      status = solveOne<5>(input, control_magnitude, control_is_c_rate, direction, duration, sample_step, parameters[p], candidate.time, candidate.terminal_voltage, one_derivative, p == 0);
    else if (nch == 8)
      status = solveOne<8>(input, control_magnitude, control_is_c_rate, direction, duration, sample_step, parameters[p], candidate.time, candidate.terminal_voltage, one_derivative, p == 0);
    else {
      assert(nch == 12); // validated against registered_spm_nch above
      status = solveOne<12>(input, control_magnitude, control_is_c_rate, direction, duration, sample_step, parameters[p], candidate.time, candidate.terminal_voltage, one_derivative, p == 0);
    }
    if (status != slide::Status::Success)
      return status;
    for (std::size_t sample = 0; sample < sample_count; ++sample)
      candidate.derivative[sample * parameters.size() + p] =
        one_derivative[sample];
  }
  output = std::move(candidate);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return slide::Status::Numerical_failure;
}

} // namespace slide::core
