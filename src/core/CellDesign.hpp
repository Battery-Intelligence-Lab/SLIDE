/**
 * @file CellDesign.hpp
 * @brief Cold physical description hierarchy and hot electrode parameters (PLAN.md §3.3).
 */

#pragma once

#include "StateArena.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace slide::core {

//!< v4's canonical domain order. It intentionally differs from the legacy unscoped enum;
//!< conversion belongs in the legacy factory and is never implicit.
enum class Domain : unsigned char { neg = 0,
                                    pos = 1 };

constexpr std::size_t domain_index(Domain domain)
{
  return static_cast<std::size_t>(domain);
}

inline constexpr std::array domains{ Domain::neg, Domain::pos };

constexpr Domain opposite(Domain domain)
{
  return domain == Domain::neg ? Domain::pos : Domain::neg;
}

constexpr int molar_flux_sign(Domain domain)
{
  return domain == Domain::neg ? 1 : -1;
}

template <class T>
using PerDomain = std::array<T, 2>;

template <class T>
constexpr T &domain_value(PerDomain<T> &values, Domain domain)
{
  return values[domain_index(domain)];
}

template <class T>
constexpr const T &domain_value(const PerDomain<T> &values, Domain domain)
{
  return values[domain_index(domain)];
}

//!< Value-semantic cold-path curve. build() resamples it to the hot uniform LUT form.
struct OCVCurve
{
  std::vector<real_t> stoichiometry{};
  std::vector<real_t> value{};

  friend bool operator==(const OCVCurve &, const OCVCurve &) = default;
};

struct Arrhenius
{
  real_t reference_value{};
  real_t activation_energy{};
  real_t reference_temperature{ 298.15 };
};

struct ActiveMaterial
{
  OCVCurve ocv{};
  real_t cs_max{};
  real_t x_0{};
  real_t x_100{};
  Arrhenius D_s{};
  Arrhenius k_ct{};
};

struct StressParams
{
  real_t youngs_modulus{};
  real_t poisson_ratio{};
  real_t partial_molar_volume{};
};

enum class AgingMechanismKind : unsigned char {
  sei_kinetic,
  sei_diffusion,
  lam_stress,
  crack_surface,
  lithium_plating,
  semi_empirical
};

struct AgingMechanismSpec
{
  AgingMechanismKind kind{};
  std::string name{};
  std::vector<real_t> coefficients{};
};

struct ElectrodeDesign
{
  ActiveMaterial active_material{};
  real_t thickness{};
  real_t porosity{};
  real_t active_fraction{};
  real_t particle_radius{};
  StressParams stress{};
  std::vector<AgingMechanismSpec> aging{};
};

struct SeparatorDesign
{
  real_t thickness{};
  real_t porosity{};
};

struct ElectrolyteDesign
{
  real_t concentration{};
  real_t diffusivity{};
  real_t transference_number{};
};

struct ThermalDesign
{
  real_t density{};
  real_t heat_capacity{};
  real_t volume{};
  real_t surface_area{};
  real_t h_conv{};
  real_t reference_temperature{ 298.15 };
  real_t environment_temperature{ 298.15 };
};

struct CellDesign
{
  PerDomain<ElectrodeDesign> electrode{};
  SeparatorDesign separator{};
  ElectrolyteDesign electrolyte{};
  ThermalDesign thermal{};
  real_t capacity_Ah{};
  real_t electrode_area{};
};

//!< Hot, build-resolved parameters used by electrode kernels. Curves/functions are represented
//!< by separate compiled LUT handles; only raw SI scalars survive here.
struct ElectrodeParams
{
  real_t cs_max{};
  real_t x_0{};
  real_t x_100{};
  real_t thickness{};
  real_t active_fraction{};
  real_t particle_radius{};
  real_t specific_surface_area{};
  real_t diffusion_ref{};
  real_t diffusion_activation{};
  real_t reaction_rate_ref{};
  real_t reaction_activation{};
};

} // namespace slide::core
