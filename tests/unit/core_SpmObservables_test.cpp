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
#include <cmath>
#include <cstdio>
#include <span>

using namespace slide;

namespace {

constexpr std::size_t coreIndex(slide::Domain domain)
{
  return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
}

template <int NCH>
core::SpmConcentrationParams<NCH> captureParams(Cell_SPM &cell, Model_SPM<NCH> &model,
                                                double Tref)
{
  auto &state = cell.getStateObj();
  core::SpmConcentrationParams<NCH> p;
  p.T_ref = Tref;
  p.elec_surf = 0.1 * 0.2 * 31; // Geometry_SPM::Acell * 31; protected in legacy
  p.D_T[coreIndex(pos)] = 29000.0;
  p.D_T[coreIndex(neg)] = 35000.0 / 5.0;
  for (auto dom : { pos, neg }) {
    const auto d = coreIndex(dom);
    p.D0[d] = state.D(dom);
    p.a[d] = state.a(dom);
    p.thick[d] = state.thick(dom);
    p.sgn[d] = sign(dom);
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

TEST_CASE("SPM concentration observable reproduces legacy Cell_SPM::getC", "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);

  Cell_SPM cell;
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(20.0, false, false) == Status::Success);

  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  auto *model = Model_SPM<>::makeModel();
  const auto params = captureParams<NCH>(cell, *model, Tref);
  const std::array legacy{ cell.getC(pos), cell.getC(neg) };

  auto &legacy_state = cell.getStateObj();
  core::BatchBuilder builder;
  const auto zp = builder.declare({ "zp", NCH, core::Unit::none });
  const auto zn = builder.declare({ "zn", NCH, core::Unit::none });
  const auto temperature = builder.declare({ "T", 1, core::Unit::K });
  auto arena = builder.build(1);
  for (int mode = 0; mode < NCH; ++mode) {
    arena.at(zp, mode, 0) = legacy_state.zp(mode);
    arena.at(zn, mode, 0) = legacy_state.zn(mode);
  }
  arena.at(temperature, 0, 0) = legacy_state.T();

  const std::array<double, 1> iapp{ cell.I() / params.elec_surf };
  std::array<double, NCH + 2> cp{}, cn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  core::computeSpmConcentrations(params, state, zp, zn, temperature, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });

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
  const auto params = captureParams<NCH>(cell, *model, Tref);
  const std::array<double, 2> cmax{ 51385.0, 30555.0 };

  core::BatchBuilder builder;
  const auto zp = builder.declare({ "zp", NCH, core::Unit::none });
  const auto zn = builder.declare({ "zn", NCH, core::Unit::none });
  const auto temperature = builder.declare({ "T", 1, core::Unit::K });
  auto arena = builder.build(lanes);
  const std::array slices{ zp, zn };

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
    arena.at(temperature, 0, lane) = cell.T() + 1.5 * lane;

  std::array<double, lanes> iapp{};
  std::array<double, (NCH + 2) * lanes> cp{}, cn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  core::computeSpmConcentrations(params, state, zp, zn, temperature, ctx, { std::span<core::real_t>{ cn }, std::span<core::real_t>{ cp } });

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
