/**
 * @file SpmScalarKernels.hpp
 * @brief Scalar-only SPM physics shared by CPU, CUDA, Dual, and WASM paths.
 *
 * This header deliberately owns expression trees, not orchestration.  Callers retain
 * indexing, caches, curve-segment selection, validation, Status publication, and rollback.
 * Keeping the functions allocation-free and container-free makes the same definitions usable
 * from ordinary C++, nvcc device code, and forward-mode scalar types.
 */

#pragma once

#include <cmath>

#if defined(__CUDACC__)
#define SLIDE_SPM_HOST_DEVICE __host__ __device__
#else
#define SLIDE_SPM_HOST_DEVICE
#endif

namespace slide::core::spm_scalar {
namespace detail {

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE inline Scalar exponential(Scalar value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::exp(value);
#else
    using std::exp;
    return exp(value); // ADL selects slide::core::exp for Dual.
#endif
  }

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE inline Scalar exponentialMinusOne(Scalar value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::expm1(value);
#else
    using std::expm1;
    return expm1(value); // ADL selects slide::core::expm1 for Dual.
#endif
  }

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE inline Scalar squareRoot(Scalar value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::sqrt(value);
#else
    using std::sqrt;
    return sqrt(value); // ADL selects slide::core::sqrt for Dual.
#endif
  }

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE inline Scalar inverseHyperbolicSine(Scalar value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::asinh(value);
#else
    using std::asinh;
    return asinh(value); // ADL selects slide::core::asinh for Dual.
#endif
  }

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE inline double primalMagnitude(Scalar value) noexcept
  {
    const double primal = static_cast<double>(value);
    return primal < 0.0 ? -primal : primal;
  }

} // namespace detail

template <class Scalar>
SLIDE_SPM_HOST_DEVICE inline Scalar arrheniusFactor(
  Scalar reference_temperature,
  Scalar temperature,
  Scalar gas_constant) noexcept
{
  return (Scalar{ 1.0 } / reference_temperature
          - Scalar{ 1.0 } / temperature)
         / gas_constant;
}

template <class Reference, class Activation, class Arrhenius>
SLIDE_SPM_HOST_DEVICE inline auto activatedValue(
  Reference reference,
  Activation activation,
  Arrhenius arrhenius) noexcept
{
  return reference * detail::exponential(activation * arrhenius);
}

template <class Area, class ElectronCount, class Faraday, class Thickness>
SLIDE_SPM_HOST_DEVICE inline auto fluxDenominator(
  Area specific_area,
  ElectronCount electron_count,
  Faraday faraday,
  Thickness thickness) noexcept
{
  return specific_area * electron_count * faraday * thickness;
}

template <class Sign, class Current, class Denominator>
SLIDE_SPM_HOST_DEVICE inline auto molarFlux(
  Sign sign,
  Current current_density,
  Denominator denominator) noexcept
{
  return sign * current_density / denominator;
}

template <class State, class Diffusivity, class Eigenvalue, class InputCoefficient, class Flux>
SLIDE_SPM_HOST_DEVICE inline auto diffusionRate(
  State state,
  Diffusivity diffusivity,
  Eigenvalue eigenvalue,
  InputCoefficient input,
  Flux flux) noexcept
{
  return diffusivity * eigenvalue * state + input * flux;
}

template <class Scalar>
SLIDE_SPM_HOST_DEVICE inline Scalar modalPhi1(Scalar x) noexcept
{
  return detail::primalMagnitude(x) < 1e-7
           ? Scalar{ 1.0 }
               + x * (0.5 + x * (1.0 / 6.0 + x / 24.0))
           : detail::exponentialMinusOne(x) / x;
}

template <class State, class Diffusivity, class Eigenvalue, class Time, class InputCoefficient, class Flux>
SLIDE_SPM_HOST_DEVICE inline auto advanceModal(
  State state,
  Diffusivity diffusivity,
  Eigenvalue eigenvalue,
  Time dt,
  InputCoefficient input,
  Flux flux) noexcept
{
  const auto x = diffusivity * eigenvalue * dt;
  const auto phi1 = modalPhi1(x);
  return detail::exponential(x) * state + dt * phi1 * input * flux;
}

template <class State, class Coefficient, class Flux, class Diffusivity>
SLIDE_SPM_HOST_DEVICE inline State concentrationOutput(
  int modes,
  const State *first_mode,
  int mode_stride,
  const Coefficient *coefficients,
  Coefficient feedthrough,
  Flux flux,
  Diffusivity diffusivity) noexcept
{
  State result{};
  for (int mode = 0; mode < modes; ++mode)
    result += coefficients[mode] * first_mode[mode * mode_stride];
  return result + feedthrough * flux / diffusivity;
}

template <class Concentration, class Maximum>
SLIDE_SPM_HOST_DEVICE inline auto surfaceStoichiometry(
  Concentration concentration,
  Maximum maximum_concentration) noexcept
{
  return concentration / maximum_concentration;
}

template <class ReactionRate, class ElectronCount, class Faraday, class ElectrolyteConcentration, class Concentration, class MaximumConcentration>
SLIDE_SPM_HOST_DEVICE inline auto exchangeCurrent(
  ReactionRate reaction_rate,
  ElectronCount electron_count,
  Faraday faraday,
  ElectrolyteConcentration electrolyte_concentration,
  Concentration concentration,
  MaximumConcentration maximum_concentration) noexcept
{
  return reaction_rate * electron_count * faraday
         * detail::squareRoot(electrolyte_concentration * concentration
                              * (maximum_concentration - concentration));
}

template <class Sign, class Current, class Area, class Thickness, class ExchangeCurrent>
SLIDE_SPM_HOST_DEVICE inline auto activationArgument(
  Sign sign,
  Current current_density,
  Area specific_area,
  Thickness thickness,
  ExchangeCurrent exchange_current) noexcept
{
  return (0.5 * sign) * current_density
         / (specific_area * thickness * exchange_current);
}

template <class Temperature, class GasConstant, class ElectronCount, class Faraday, class Argument>
SLIDE_SPM_HOST_DEVICE inline auto activationOverpotential(
  Temperature temperature,
  GasConstant gas_constant,
  ElectronCount electron_count,
  Faraday faraday,
  Argument argument) noexcept
{
  return 2.0 * gas_constant * temperature / (electron_count * faraday)
         * detail::inverseHyperbolicSine(argument);
}

template <class Query, class Coordinate, class Value>
SLIDE_SPM_HOST_DEVICE inline auto linearInterpolate(
  Query query,
  Coordinate x0,
  Coordinate x1,
  Value y0,
  Value y1) noexcept
{
  return y0 + (y1 - y0) * (query - x0) / (x1 - x0);
}

template <class SpecificArea, class ElectrodeArea, class Thickness>
SLIDE_SPM_HOST_DEVICE inline auto activeArea(
  SpecificArea specific_area,
  ElectrodeArea electrode_area,
  Thickness thickness) noexcept
{
  return specific_area * electrode_area * thickness;
}

template <class SeiThickness, class SeiResistivityArea, class NegativeSpecificResistance, class PositiveSpecificResistance, class CollectorResistanceArea, class NegativeArea, class PositiveArea, class ElectrodeArea>
SLIDE_SPM_HOST_DEVICE inline auto seriesResistance(
  SeiThickness sei_thickness,
  SeiResistivityArea sei_resistivity_area,
  NegativeSpecificResistance negative_specific_resistance,
  PositiveSpecificResistance positive_specific_resistance,
  CollectorResistanceArea collector_resistance_area,
  NegativeArea negative_area,
  PositiveArea positive_area,
  ElectrodeArea electrode_area) noexcept
{
  return sei_thickness * sei_resistivity_area / negative_area
         + negative_specific_resistance / negative_area
         + positive_specific_resistance / positive_area
         + collector_resistance_area / electrode_area;
}

template <class NegativeOcv, class PositiveOcv, class Temperature, class ReferenceTemperature, class EntropicCoefficient>
SLIDE_SPM_HOST_DEVICE inline auto cellOpenCircuitVoltage(
  NegativeOcv negative_ocv,
  PositiveOcv positive_ocv,
  Temperature temperature,
  ReferenceTemperature reference_temperature,
  EntropicCoefficient entropic_coefficient) noexcept
{
  return positive_ocv - negative_ocv
         + (temperature - reference_temperature) * entropic_coefficient;
}

template <class Ocv, class NegativeOverpotential, class PositiveOverpotential, class Resistance, class Current>
SLIDE_SPM_HOST_DEVICE inline auto terminalVoltage(
  Ocv ocv,
  NegativeOverpotential negative_overpotential,
  PositiveOverpotential positive_overpotential,
  Resistance resistance,
  Current current) noexcept
{
  return ocv + positive_overpotential - negative_overpotential
         - resistance * current;
}

} // namespace slide::core::spm_scalar

#undef SLIDE_SPM_HOST_DEVICE
