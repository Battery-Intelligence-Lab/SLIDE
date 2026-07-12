/**
 * @file SpmScalarKernels.hpp
 * @brief Scalar-only SPM physics shared by CPU, CUDA, Dual, and ordinary-C++ paths.
 *
 * M0.5 / PC-10 hot-path contract: this header is the sole owner of SPM scalar
 * expression trees. Adapters own only indexing, caches, curve-segment selection,
 * validation, Status publication, and rollback. The allocation-free, container-free
 * definitions instantiate for CPU, CUDA, and Dual; WASM uses the ordinary C++
 * instantiation, with independent wasm/native validation deferred to M10.
 */

#pragma once

#include <cmath>
#include <type_traits>

#if defined(__CUDACC__)
#define SLIDE_SPM_HOST_DEVICE __host__ __device__
#else
#define SLIDE_SPM_HOST_DEVICE
#endif

#if defined(__CUDACC__)
#define SLIDE_SPM_FORCE_INLINE __forceinline__
#elif defined(_MSC_VER)
#define SLIDE_SPM_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SLIDE_SPM_FORCE_INLINE inline __attribute__((always_inline))
#else
#define SLIDE_SPM_FORCE_INLINE inline
#endif

namespace slide::core::spm_scalar {
namespace detail {

  SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE double exponential(double value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::exp(value);
#else
    return std::exp(value);
#endif
  }

  template <class Scalar>
  inline Scalar exponential(Scalar value) noexcept
  {
    using std::exp;
    return exp(value); // ADL selects slide::core::exp for Dual.
  }

  SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE double exponentialMinusOne(double value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::expm1(value);
#else
    return std::expm1(value);
#endif
  }

  template <class Scalar>
  inline Scalar exponentialMinusOne(Scalar value) noexcept
  {
    using std::expm1;
    return expm1(value); // ADL selects slide::core::expm1 for Dual.
  }

  SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE double squareRoot(double value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::sqrt(value);
#else
    return std::sqrt(value);
#endif
  }

  template <class Scalar>
  inline Scalar squareRoot(Scalar value) noexcept
  {
    using std::sqrt;
    return sqrt(value); // ADL selects slide::core::sqrt for Dual.
  }

  SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE double inverseHyperbolicSine(double value) noexcept
  {
#if defined(__CUDA_ARCH__)
    return ::asinh(value);
#else
    return std::asinh(value);
#endif
  }

  template <class Scalar>
  inline Scalar inverseHyperbolicSine(Scalar value) noexcept
  {
    using std::asinh;
    return asinh(value); // ADL selects slide::core::asinh for Dual.
  }

  template <class Scalar>
  SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE double primalMagnitude(Scalar value) noexcept
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

/**
 * Modal RHS expression kernel for vectorized consumers whose optimizer changes
 * loop IR across an inline-function boundary. Each argument is evaluated once.
 */
#define SLIDE_SPM_DIFFUSION_RATE(state, diffusivity, eigenvalue, input, flux) \
  ((diffusivity) * (eigenvalue) * (state) + (input) * (flux))

template <class State, class Diffusivity, class Eigenvalue, class InputCoefficient, class Flux>
SLIDE_SPM_HOST_DEVICE inline auto diffusionRate(
  State state,
  Diffusivity diffusivity,
  Eigenvalue eigenvalue,
  InputCoefficient input,
  Flux flux) noexcept
{
  return SLIDE_SPM_DIFFUSION_RATE(
    state, diffusivity, eigenvalue, input, flux);
}

/**
 * Exact modal update, expanded into the consumer's loop.
 *
 * Release builds intentionally permit reassociation and vectorisation.  A normal inline-call
 * boundary changed those choices under Clang -Ofast even with always_inline, violating the
 * recorded-bit contract.  This one hygienically-prefixed statement kernel keeps the expression
 * tree in one source while preserving each backend's existing math entry points and loop IR.
 * `destination` and `dt` are each evaluated twice and must be side-effect-free;
 * callable arguments must be function-name tokens, and the remaining arguments values.
 */
#define SLIDE_SPM_ADVANCE_MODAL_IN_PLACE(                                                            \
  destination, diffusivity, eigenvalue, dt, input, flux, abs_function, exp_function, expm1_function) \
  do {                                                                                               \
    const auto slide_spm_modal_x = (diffusivity) * (eigenvalue) * (dt);                              \
    const auto slide_spm_modal_phi1 =                                                                \
      abs_function(slide_spm_modal_x) < 1e-7                                                         \
        ? decltype(slide_spm_modal_x){ 1.0 }                                                         \
            + slide_spm_modal_x                                                                      \
                * (0.5                                                                               \
                   + slide_spm_modal_x                                                               \
                       * (1.0 / 6.0 + slide_spm_modal_x / 24.0))                                     \
        : expm1_function(slide_spm_modal_x) / slide_spm_modal_x;                                     \
    (destination) = exp_function(slide_spm_modal_x) * (destination)                                  \
                    + (dt) * slide_spm_modal_phi1 * (input) * (flux);                                \
  } while (false)

#define SLIDE_SPM_ADVANCE_MODAL_STD(destination, diffusivity, eigenvalue, dt, input, flux) \
  SLIDE_SPM_ADVANCE_MODAL_IN_PLACE(                                                        \
    destination, diffusivity, eigenvalue, dt, input, flux, std::abs, std::exp, std::expm1)

#define SLIDE_SPM_ADVANCE_MODAL_ADL(destination, diffusivity, eigenvalue, dt, input, flux) \
  SLIDE_SPM_ADVANCE_MODAL_IN_PLACE(                                                        \
    destination, diffusivity, eigenvalue, dt, input, flux, ::slide::core::spm_scalar::detail::primalMagnitude, exp, expm1)

#define SLIDE_SPM_ADVANCE_MODAL_CUDA(destination, diffusivity, eigenvalue, dt, input, flux) \
  SLIDE_SPM_ADVANCE_MODAL_IN_PLACE(                                                         \
    destination, diffusivity, eigenvalue, dt, input, flux, ::fabs, ::exp, ::expm1)

template <class State, class Diffusivity, class Eigenvalue, class Time, class InputCoefficient, class Flux>
SLIDE_SPM_HOST_DEVICE SLIDE_SPM_FORCE_INLINE auto advanceModal(
  State state,
  Diffusivity diffusivity,
  Eigenvalue eigenvalue,
  Time dt,
  InputCoefficient input,
  Flux flux) noexcept
{
  State result = state;
#if defined(__CUDACC__)
  SLIDE_SPM_ADVANCE_MODAL_CUDA(result, diffusivity, eigenvalue, dt, input, flux);
#else
  using Scalar = std::remove_cv_t<decltype(diffusivity * eigenvalue * dt)>;
  if constexpr (std::is_same_v<Scalar, double>)
    SLIDE_SPM_ADVANCE_MODAL_STD(result, diffusivity, eigenvalue, dt, input, flux);
  else
    SLIDE_SPM_ADVANCE_MODAL_ADL(result, diffusivity, eigenvalue, dt, input, flux);
#endif
  return result;
}

template <class ModalSum, class Coefficient, class Flux, class Diffusivity>
SLIDE_SPM_HOST_DEVICE inline auto concentrationOutput(
  ModalSum modal_sum,
  Coefficient feedthrough,
  Flux flux,
  Diffusivity diffusivity) noexcept
{
  return modal_sum + feedthrough * flux / diffusivity;
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
  return concentrationOutput(result, feedthrough, flux, diffusivity);
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

#undef SLIDE_SPM_FORCE_INLINE
#undef SLIDE_SPM_HOST_DEVICE
