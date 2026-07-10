/**
 * @file core_ChebyshevTransient_test.cpp
 * @brief P1-G3(b): constant-flux spherical-diffusion transient vs analytic series.
 *
 * This oracle exercises Model_SPM A/B/C/D jointly and the centre Cc/cc_coeff path through
 * slide::core::computeSpmConcentrations. For dimensionless radius x and time tau=D*t/R^2,
 * with the SLIDE convention D*dc/dr|R = -j (positive j depletes the particle) and initially
 * uniform c0, Carslaw-Jaeger/Crank gives
 *
 * (c-c0)/(jR/D) = -[3tau + (5x^2-3)/10
 *                  - (2/x) sum sin(mu_k*x)/(mu_k^2 sin(mu_k)) exp(-mu_k^2 tau)],
 * tan(mu_k)=mu_k. At x=0, sin(mu*x)/x is evaluated by its limit mu.
 *
 * REGISTERED BEFORE FIRST RUN: at nch={5,8,12}, tau={0.2,1.0}, both electrodes and every
 * surface/interior/centre node, relative error <= 1e-6. The imposed gradient is nonzero and
 * centre != surface. This is independent of legacy trajectories and validates the output maps.
 * OUTCOME 2026-07-10: band holds. Worst max relative error at tau=0.2 is 8.254e-7,
 * 2.356e-11 and 3.949e-12 for nch=5,8,12; across Debug/Release tau=1.0 is <=4.661e-13.
 */

#include "../../src/cells/Cell_SPM/Model_SPM.hpp"
#include "../../src/core/BatchBuilder.hpp"
#include "../../src/core/SpmObservables.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <span>
#include <tuple>

namespace {

template <int N>
struct NchTag
{
  static constexpr int value = N;
};

using NchValues = std::tuple<NchTag<5>, NchTag<8>, NchTag<12>>;

double mu_root(int k)
{
  const auto g = [](double mu) { return std::sin(mu) - mu * std::cos(mu); };
  double lo = k * std::numbers::pi;
  double hi = (k + 0.5) * std::numbers::pi;
  double glo = g(lo);
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double mid = 0.5 * (lo + hi);
    const double gmid = g(mid);
    if (glo * gmid <= 0.0) {
      hi = mid;
    } else {
      lo = mid;
      glo = gmid;
    }
  }
  return 0.5 * (lo + hi);
}

double analytic_increment(double x, double tau)
{
  double series = 0.0;
  for (int k = 1; k <= 64; ++k) {
    const double mu = mu_root(k);
    const double radial = (x == 0.0) ? mu : std::sin(mu * x) / x;
    series += radial / (mu * mu * std::sin(mu)) * std::exp(-mu * mu * tau);
  }
  return -(3.0 * tau + (5.0 * x * x - 3.0) / 10.0 - 2.0 * series);
}

template <int NCH>
slide::core::SpmConcentrationParams<NCH> make_params(const slide::Model_SPM<NCH> &model,
                                                     double diffusivity)
{
  slide::core::SpmConcentrationParams<NCH> p;
  p.F = 1.0;
  p.Rg = 1.0;
  p.n = 1.0;
  p.T_ref = 1.0;
  p.cc_coeff = model.cc_coeff;
  for (int node = 0; node < NCH + 1; ++node)
    p.Cc[node] = model.Cc(node);

  for (int dom = 0; dom < 2; ++dom) {
    p.R[dom] = (dom == slide::pos) ? model.Rp : model.Rn;
    p.D0[dom] = diffusivity;
    p.D_T[dom] = 0.0;
    p.a[dom] = 1.0;
    p.thick[dom] = 1.0;
    p.sgn[dom] = 1;
    for (int node = 0; node < NCH + 1; ++node) {
      p.Dout[dom][node] = model.D[dom](node);
      for (int mode = 0; mode < NCH; ++mode)
        p.C[dom][node][mode] = model.C[dom](node, mode);
    }
  }
  return p;
}

} // namespace

TEMPLATE_LIST_TEST_CASE("Chebyshev constant-flux transient matches analytic spherical series",
                        "[core][Chebyshev][P1-G3]", NchValues)
{
  constexpr int NCH = TestType::value;
  constexpr double diffusivity = 1e-14;
  constexpr double flux = 1e-9;
  constexpr double initial_concentration = 10.0;
  constexpr std::array tau_values{ 0.2, 1.0 };

  slide::Model_SPM<NCH> model;
  const auto params = make_params(model, diffusivity);

  slide::core::BatchBuilder builder;
  const auto zp = builder.declare({ "zp", NCH, slide::core::Unit::none });
  const auto zn = builder.declare({ "zn", NCH, slide::core::Unit::none });
  const auto temperature = builder.declare({ "T", 1, slide::core::Unit::K });
  auto arena = builder.build(1);
  const std::array slices{ zp, zn };
  arena.at(temperature, 0, 0) = 1.0;

  const std::array<double, 1> iapp{ flux };
  const slide::core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  std::array<double, NCH + 2> cp{}, cn{};

  for (const double tau : tau_values) {
    for (int dom = 0; dom < 2; ++dom) {
      const double radius = params.R[dom];
      Eigen::Vector<double, NCH> u;
      for (int node = 0; node < NCH; ++node)
        u(node) = radius * initial_concentration * model.xch(node);
      const Eigen::Vector<double, NCH> z0 = model.V[dom] * u;
      const double physical_time = tau * radius * radius / diffusivity;

      for (int mode = 0; mode < NCH; ++mode) {
        const double rate = diffusivity * model.A[dom](mode);
        const double source = model.B[dom](mode) * flux;
        const double z = (rate == 0.0)
                           ? z0(mode) + source * physical_time
                           : std::exp(rate * physical_time) * z0(mode)
                               + source * std::expm1(rate * physical_time) / rate;
        arena.at(slices[dom], mode, 0) = z;
      }
    }

    const slide::core::ConstBatchView state{
      slide::core::BatchShape::from(arena),
      std::span<const slide::core::real_t>{ arena.raw() }
    };
    slide::core::computeSpmConcentrations(
      params, state, zp, zn, temperature, ctx, { std::span<slide::core::real_t>{ cp }, std::span<slide::core::real_t>{ cn } });

    const std::array output{ cp, cn };
    double max_relative_error = 0.0;
    for (int dom = 0; dom < 2; ++dom) {
      const double radius = params.R[dom];
      const double scale = flux * radius / diffusivity;
      for (int node = 0; node < NCH + 2; ++node) {
        const double x = (node == 0)         ? 1.0
                         : (node == NCH + 1) ? 0.0
                                             : model.xch(node - 1);
        const double expected = analytic_increment(x, tau);
        const double actual = (output[dom][node] - initial_concentration) / scale;
        const double relative_error = std::abs(actual - expected) / std::abs(expected);
        max_relative_error = std::max(max_relative_error, relative_error);
        CAPTURE(NCH, tau, dom, node, x, actual, expected, relative_error);
        REQUIRE(std::isfinite(actual));
        REQUIRE(relative_error <= 1e-6);
      }
      REQUIRE(output[dom][0] != output[dom][NCH + 1]);
    }

    std::printf("Cheb-transient nch=%d tau=%.1f max_rel=%.3e\n",
                NCH,
                tau,
                max_relative_error);
  }
}
