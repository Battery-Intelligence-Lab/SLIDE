/**
 * @file SpmStress.hpp
 * @brief Shared Dai/Laresgoiti stress observables and explicit history state.
 */

#pragma once

#include "AgeingKernel.hpp"
#include "CompiledCurve.hpp"
#include "SpmObservables.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <span>

namespace slide::core {

/** Previous-step values consumed by crack/LAM rate laws; all are checkpointed arena rows. */
struct StressHistoryLayout
{
  PerDomain<StateSlice> previous_dai{};
  StateSlice previous_laresgoiti_negative{};
  StateSlice interval{};
};

inline StressHistoryLayout declareStressHistory(BatchBuilder &builder)
{
  StressHistoryLayout layout;
  domain_value(layout.previous_dai, Domain::neg) = builder.declare({ "previous_dai_stress_neg", 1, Unit::Pa, StateRole::algebraic });
  domain_value(layout.previous_dai, Domain::pos) = builder.declare({ "previous_dai_stress_pos", 1, Unit::Pa, StateRole::algebraic });
  layout.previous_laresgoiti_negative = builder.declare({ "previous_laresgoiti_stress_neg", 1, Unit::Pa, StateRole::algebraic });
  layout.interval = builder.declare({ "stress_interval", 1, Unit::s, StateRole::algebraic });
  return layout;
}

template <int NCH>
struct SpmStressParams
{
  static constexpr int full_nodes = 2 * NCH + 3;
  std::array<real_t, NCH> x_inner{};
  std::array<std::array<real_t, full_nodes>, full_nodes> integration{};
  PerDomain<real_t> partial_molar_volume{};
  PerDomain<real_t> youngs_modulus{};
  PerDomain<real_t> poisson_ratio{};
  IndexedPiecewiseLinear laresgoiti_negative{};
};

template <int NCH>
[[nodiscard]] inline slide::Status validateSpmStressParams(const SpmStressParams<NCH> &p)
{
  if (!p.laresgoiti_negative.valid())
    return slide::Status::Invalid_parameters;
  for (const Domain domain : domains) {
    [[maybe_unused]] const auto d = domain_index(domain);
    if (!(is_finite(p.partial_molar_volume[d])
          && is_finite(p.youngs_modulus[d]) && p.youngs_modulus[d] > 0.0
          && is_finite(p.poisson_ratio[d]) && p.poisson_ratio[d] < 1.0))
      return slide::Status::Invalid_parameters;
  }
  return slide::Status::Success;
}

template <class Real>
struct BasicSpmStress
{
  PerDomain<std::span<Real>> dai_maximum_hydrostatic{};
  std::span<Real> laresgoiti_negative{};
};

template <class Real = real_t>
class SpmStressScratch
{
public:
  explicit SpmStressScratch(int n_lanes)
    : storage_{ n_lanes }
  {}

  BasicSpmStress<Real> view()
  {
    BasicSpmStress<Real> result;
    result.dai_maximum_hydrostatic[domain_index(Domain::neg)] = storage_.field(0);
    result.dai_maximum_hydrostatic[domain_index(Domain::pos)] = storage_.field(1);
    result.laresgoiti_negative = storage_.field(2);
    return result;
  }

private:
  detail::AgeingScratchStorage<Real, 3> storage_;
};

/** Reconstruct maximum Dai hydrostatic stress and Laresgoiti graphite stress per lane. */
template <int NCH, class Real>
[[nodiscard]] slide::Status computeSpmStress(const SpmStressParams<NCH> &p,
                                             const BasicSpmObservables<Real> &observables,
                                             int n_lanes,
                                             BasicSpmStress<Real>
                                               output)
{
  constexpr int half_nodes = NCH + 2;
  constexpr int full_nodes = 2 * NCH + 3;
  constexpr int centre = NCH + 1;
  const auto L = static_cast<std::size_t>(n_lanes);
  for (const Domain domain : domains) {
    [[maybe_unused]] const auto d = domain_index(domain);
    assert(observables.concentration[d].size() == static_cast<std::size_t>(half_nodes) * L);
    assert(output.dai_maximum_hydrostatic[d].size() == L);
  }
  assert(observables.surface_stoichiometry[domain_index(Domain::neg)].size() == L
         && output.laresgoiti_negative.size() == L);

  std::array<real_t, half_nodes> positive_x{};
  positive_x[0] = 1.0;
  for (int i = 0; i < NCH; ++i)
    positive_x[static_cast<std::size_t>(i + 1)] = p.x_inner[static_cast<std::size_t>(i)];
  positive_x[static_cast<std::size_t>(NCH + 1)] = 0.0;
  std::array<real_t, full_nodes> full_x{};
  full_x[centre] = 0.0;
  for (int i = 0; i < NCH + 1; ++i) {
    full_x[static_cast<std::size_t>(i)] = -positive_x[static_cast<std::size_t>(i)];
    full_x[static_cast<std::size_t>(NCH + 2 + i)] = positive_x[static_cast<std::size_t>(NCH - i)];
  }

  for (int lane = 0; lane < n_lanes; ++lane) {
    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      std::array<Real, full_nodes> concentration{};
      concentration[centre] = observables.concentration[d][static_cast<std::size_t>(NCH + 1) * L + lane];
      for (int i = 0; i < NCH + 1; ++i) {
        concentration[static_cast<std::size_t>(i)] = observables.concentration[d][static_cast<std::size_t>(i) * L + lane];
        concentration[static_cast<std::size_t>(NCH + 2 + i)] = observables.concentration[d][static_cast<std::size_t>(NCH - i) * L + lane];
      }

      std::array<Real, full_nodes> integral{};
      for (int row = 0; row < full_nodes; ++row)
        for (int column = 0; column < full_nodes; ++column)
          integral[static_cast<std::size_t>(row)] += p.integration[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                                                     * concentration[static_cast<std::size_t>(column)]
                                                     * full_x[static_cast<std::size_t>(column)]
                                                     * full_x[static_cast<std::size_t>(column)];

      const Real total = integral.back() - integral[centre];
      const Real factor = p.partial_molar_volume[d] * p.youngs_modulus[d]
                          / (Real{ 1 } - p.poisson_ratio[d]);
      Real maximum{};
      for (int i = 0; i < half_nodes; ++i) {
        Real radial{}, tangential{};
        if (i == 0) {
          radial = Real{ 2 } * factor / Real{ 9 }
                   * (Real{ 3 } * total - concentration[centre]);
          tangential = radial;
        } else {
          const Real x = full_x[static_cast<std::size_t>(centre + i)];
          const Real partial = (integral[static_cast<std::size_t>(centre + i)]
                                - integral[centre])
                               / (x * x * x);
          radial = Real{ 2 } * factor / Real{ 3 } * (total - partial);
          // Preserve the legacy node association exactly (Cell_SPM_degradation.cpp:695).
          const Real local = observables.concentration[d][static_cast<std::size_t>(i) * L + lane];
          tangential = factor / Real{ 3 }
                       * (Real{ 2 } * total + partial - local);
        }
        const Real hydrostatic = (radial + Real{ 2 } * tangential) / Real{ 3 };
        if (std::abs(primal_value(hydrostatic)) > std::abs(primal_value(maximum)))
          maximum = hydrostatic;
      }
      if (!is_finite_primal(maximum))
        return slide::Status::Numerical_failure;
      output.dai_maximum_hydrostatic[d][static_cast<std::size_t>(lane)] = maximum;
    }

    const Real z_negative = observables.surface_stoichiometry[domain_index(Domain::neg)][static_cast<std::size_t>(lane)];
    output.laresgoiti_negative[static_cast<std::size_t>(lane)] = p.laresgoiti_negative.eval(z_negative);
  }
  return slide::Status::Success;
}

} // namespace slide::core
