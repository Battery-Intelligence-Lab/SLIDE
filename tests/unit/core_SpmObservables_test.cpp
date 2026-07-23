/**
 * @file core_SpmObservables_test.cpp
 * @brief Phase-1 correctness gate for the SPM particle-concentration observable.
 *
 * The observable reconstructs surface, interior and centre concentrations from arena z-rows.
 * The centre Cc/cc_coeff output path carried the historical nch!=5 bug (PLAN.md §2.2).
 *
 * REGISTERED BANDS (written before the first run):
 *  - Legacy parity: every output node for both electrodes has relative error <= 1e-12.
 *    OUTCOME 2026-07-10: Debug max relative error = 0 (bit-identical); Release = 2.665e-15.
 *  - Independent uniform-profile round trip: every output node in eight distinct lanes has
 *    relative error <= 1e-8. OUTCOME 2026-07-10: max relative error = 9.027e-15;
 *    positive-surface lane spread = 7193.9 mol m^-3.
 */

#include "../../src/slide.hpp"
#include "../../src/core/BatchBuilder.hpp"
#include "../../src/core/SpmObservables.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <span>

using namespace slide;

namespace {

constexpr std::size_t coreIndex(slide::Domain domain)
{
  return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
}

template <int NCH>
core::SpmConcentrationParams<NCH> captureParams(Model_SPM<NCH> &model, double Tref)
{
  core::SpmConcentrationParams<NCH> p;
  p.T_ref = Tref;
  p.D_T[coreIndex(pos)] = 29000.0;
  p.D_T[coreIndex(neg)] = 35000.0 / 5.0;
  for (auto dom : { pos, neg }) {
    const auto d = coreIndex(dom);
    p.R[d] = (dom == pos) ? model.Rp : model.Rn;
    for (int node = 0; node < NCH + 1; ++node) {
      p.Dout[d][node] = model.D[dom](node);
      for (int mode = 0; mode < NCH; ++mode)
        p.C[d][node][mode] = model.C[dom](node, mode);
    }
  }
  for (int node = 0; node < NCH + 1; ++node)
    p.Cc[node] = model.Cc(node);
  p.cc_coeff = model.cc_coeff;
  return p;
}

} // namespace

TEST_CASE("SPM transport cache hits and invalidations match the cache-free path",
          "[core][observables][cache]")
{
  // Registered before the first run: the cold miss, unchanged-key hit, and
  // temperature/D-reference/area/thickness invalidations must reproduce both
  // cache-free transport outputs bit-for-bit in every lane and domain.
  constexpr int NCH = 1;
  constexpr int lanes = 2;
  constexpr auto neg = core::domain_index(core::Domain::neg);
  constexpr auto pos = core::domain_index(core::Domain::pos);
  core::SpmConcentrationParams<NCH> params;
  params.T_ref = 298.15;
  params.D_T[neg] = 7000.0;
  params.D_T[pos] = 29000.0;

  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(lanes);
  for (int lane = 0; lane < lanes; ++lane)
    arena.at(layout.temperature, 0, lane) =
      295.0 + 5.0 * static_cast<double>(lane);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane) {
      const auto d_value = static_cast<double>(d);
      const auto lane_value = static_cast<double>(lane);
      arena.at(layout.diffusion_coefficient[d], 0, lane) =
        1.0e-14 * (1.0 + d_value + 0.1 * lane_value);
      arena.at(layout.specific_surface_area[d], 0, lane) =
        2.0 + d_value + 0.25 * lane_value;
      arena.at(layout.electrode_thickness[d], 0, lane) =
        1.0e-4 * (1.0 + d_value + 0.5 * lane_value);
    }
  }

  constexpr std::array current_density{ -2.5, 4.0 };
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.25, .i_app = current_density };
  core::SpmTransportCache cache{ lanes };
  core::PerDomain<std::array<double, lanes>> cached_diffusivity{};
  core::PerDomain<std::array<double, lanes>> cached_flux{};
  core::PerDomain<std::array<double, lanes>> direct_diffusivity{};
  core::PerDomain<std::array<double, lanes>> direct_flux{};

  const auto require_cache_matches = [&] {
    for (const auto domain : core::domains) {
      const auto d = core::domain_index(domain);
      cached_diffusivity[d].fill(-101.0);
      cached_flux[d].fill(202.0);
      direct_diffusivity[d].fill(-303.0);
      direct_flux[d].fill(404.0);
    }
    core::computeSpmTransport(
      params,
      state,
      layout,
      ctx,
      { std::span{ direct_diffusivity[neg] },
        std::span{ direct_diffusivity[pos] } },
      { std::span{ direct_flux[neg] }, std::span{ direct_flux[pos] } });
    core::computeSpmTransport(
      params,
      state,
      layout,
      ctx,
      { std::span{ cached_diffusivity[neg] },
        std::span{ cached_diffusivity[pos] } },
      { std::span{ cached_flux[neg] }, std::span{ cached_flux[pos] } },
      &cache);
    for (const auto domain : core::domains) {
      const auto d = core::domain_index(domain);
      for (int lane = 0; lane < lanes; ++lane) {
        const auto i = static_cast<std::size_t>(lane);
        CHECK(std::bit_cast<std::uint64_t>(cached_diffusivity[d][i])
              == std::bit_cast<std::uint64_t>(direct_diffusivity[d][i]));
        CHECK(std::bit_cast<std::uint64_t>(cached_flux[d][i])
              == std::bit_cast<std::uint64_t>(direct_flux[d][i]));
      }
    }
  };

  require_cache_matches(); // cold miss and store
  require_cache_matches(); // unchanged-key hit
  arena.at(layout.temperature, 0, 1) += 0.5;
  require_cache_matches();
  arena.at(layout.diffusion_coefficient[neg], 0, 0) *= 1.25;
  require_cache_matches();
  arena.at(layout.specific_surface_area[pos], 0, 1) *= 1.5;
  require_cache_matches();
  arena.at(layout.electrode_thickness[neg], 0, 1) *= 0.75;
  require_cache_matches();
}

TEST_CASE("SPM concentration observable reproduces legacy Cell_SPM::getC", "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);

  Cell_SPM cell;
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(20.0, false, false) == Status::Success);

  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  auto *model = Model_SPM<>::makeModel();
  const auto params = captureParams<NCH>(*model, Tref);
  const std::array legacy{ cell.getC(pos), cell.getC(neg) };

  auto &legacy_state = cell.getStateObj();
  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1);
  for (int mode = 0; mode < NCH; ++mode) {
    arena.at(core::domain_value(layout.z, core::Domain::pos), mode, 0) = legacy_state.zp(mode);
    arena.at(core::domain_value(layout.z, core::Domain::neg), mode, 0) = legacy_state.zn(mode);
  }
  arena.at(layout.temperature, 0, 0) = legacy_state.T();
  for (const auto legacy_domain : { pos, neg }) {
    const auto d = coreIndex(legacy_domain);
    arena.at(layout.diffusion_coefficient[d], 0, 0) = legacy_state.D(legacy_domain);
    arena.at(layout.specific_surface_area[d], 0, 0) = legacy_state.a(legacy_domain);
    arena.at(layout.electrode_thickness[d], 0, 0) = legacy_state.thick(legacy_domain);
  }

  constexpr double electrode_area = 0.1 * 0.2 * 31;
  const std::array<double, 1> iapp{ cell.I() / electrode_area };
  std::array<double, NCH + 2> cp{}, cn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  core::computeSpmConcentrations(params, state, layout, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });

  const std::array output{ cp, cn };
  double max_rel = 0.0;
  for (int dom = 0; dom < 2; ++dom)
    for (int node = 0; node < NCH + 2; ++node) {
      REQUIRE(std::isfinite(output[dom][node]));
      max_rel = std::max(max_rel, std::abs(output[dom][node] - legacy[dom][node]) / std::abs(legacy[dom][node]));
    }

  std::printf("SpmObservables parity: max_rel(all nodes)=%.3e  csp=%.12g  csn=%.12g\n",
              max_rel,
              cp[0],
              cn[0]);
  CAPTURE(max_rel);
  REQUIRE(max_rel <= 1e-12);

  // Degradation changes D, a and thickness in State_SPM. The observable must therefore read
  // them from the trial-vector arena, not from immutable archetype parameters.
  const auto pos_index = coreIndex(pos);
  const double baseline_surface = cp[0];
  double modal_surface = 0.0;
  for (int mode = 0; mode < NCH; ++mode)
    modal_surface += params.C[pos_index][0][mode]
                     * arena.at(layout.z[pos_index], mode, 0);
  const double expected_half_correction = modal_surface
                                          + 0.5 * (baseline_surface - modal_surface);

  arena.at(layout.diffusion_coefficient[pos_index], 0, 0) *= 2.0;
  core::computeSpmConcentrations(params, state, layout, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });
  REQUIRE(std::abs(cp[0] - expected_half_correction)
          <= 1e-12 * std::abs(expected_half_correction));

  arena.at(layout.diffusion_coefficient[pos_index], 0, 0) = legacy_state.D(pos);
  arena.at(layout.specific_surface_area[pos_index], 0, 0) *= 2.0;
  core::computeSpmConcentrations(params, state, layout, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });
  REQUIRE(std::abs(cp[0] - expected_half_correction)
          <= 1e-12 * std::abs(expected_half_correction));
}

TEST_CASE("SPM concentration observable uniform round trip on heterogeneous lanes",
          "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr int lanes = 8;

  Cell_SPM cell;
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(0.0, false, false) == Status::Success);
  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  auto *model = Model_SPM<>::makeModel();
  const auto params = captureParams<NCH>(*model, Tref);
  const std::array<double, 2> cmax{ 51385.0, 30555.0 };

  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(lanes);
  const std::array slices{ core::domain_value(layout.z, core::Domain::pos),
                           core::domain_value(layout.z, core::Domain::neg) };

  std::array<std::array<double, lanes>, 2> expected{};
  for (int dom = 0; dom < 2; ++dom)
    for (int lane = 0; lane < lanes; ++lane) {
      expected[dom][lane] = (0.30 + 0.02 * lane) * cmax[dom];
      double uniform_mode = 0.0;
      for (int node = 0; node < NCH; ++node)
        uniform_mode += model->V[dom](model->zero, node)
                        * (params.R[coreIndex(static_cast<slide::Domain>(dom))]
                           * expected[dom][lane] * model->xch(node));
      for (int mode = 0; mode < NCH; ++mode)
        arena.at(slices[dom], mode, lane) = (mode == model->zero) ? uniform_mode : 0.0;
    }
  for (int lane = 0; lane < lanes; ++lane)
    arena.at(layout.temperature, 0, lane) = cell.T() + 1.5 * lane;
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane) {
      arena.at(layout.diffusion_coefficient[d], 0, lane) = cell.getStateObj().D(domain == core::Domain::pos ? pos : neg);
      arena.at(layout.specific_surface_area[d], 0, lane) = 1.0;
      arena.at(layout.electrode_thickness[d], 0, lane) = 1.0;
    }
  }

  std::array<double, lanes> iapp{};
  std::array<double, (NCH + 2) * lanes> cp{}, cn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  core::computeSpmConcentrations(params, state, layout, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });

  const std::array output{ cp, cn };
  double max_rel = 0.0;
  for (int dom = 0; dom < 2; ++dom)
    for (int node = 0; node < NCH + 2; ++node)
      for (int lane = 0; lane < lanes; ++lane) {
        const auto index = static_cast<std::size_t>(node) * lanes + lane;
        max_rel = std::max(max_rel, std::abs(output[dom][index] - expected[dom][lane]) / expected[dom][lane]);
      }

  double spread = 0.0;
  for (int lane = 0; lane < lanes; ++lane)
    spread = std::max(spread, std::abs(cp[lane] - cp[0]));
  std::printf("SpmObservables round-trip: max_rel(all nodes)=%.3e  lane_spread(csp)=%.6g\n",
              max_rel,
              spread);

  CAPTURE(max_rel, spread);
  REQUIRE(spread > 0.0);
  REQUIRE(max_rel <= 1e-8);
}
