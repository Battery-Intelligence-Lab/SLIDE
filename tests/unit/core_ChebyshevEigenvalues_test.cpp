/**
 * @file core_ChebyshevEigenvalues_test.cpp
 * @brief P1-G3 (part a): ANALYTIC eigenvalue oracle for the Chebyshev solid-diffusion operator.
 *
 * Independent-mathematics validation of slide::Model_SPM's spectral discretisation, per the
 * math audit (.claude/reports/chebyshev-math-audit-2026-07-09.md §3). The folded + surface-
 * condensed dimensionless operator A_ (Model_SPM.hpp:123) discretises
 *     u_t = u_xx,   u(0) = 0,   u_x(1) = u(1)   (homogeneous / zero-flux surface)
 * on x in [0,1] with u = r*c ODD. Its eigenpairs are u_k = sin(mu_k * x), lambda_k = -mu_k^2,
 * where mu_k are the positive roots of  tan(mu) = mu :
 *     mu_0 = 0  (the MASS mode -> the single zero eigenvalue Model_SPM forces to 0),
 *     mu_1 ~ 4.4934,  mu_2 ~ 7.7253,  mu_3 ~ 10.9041, ...
 * Model_SPM stores A[dom] = eigenvalues of  A_ / R^2  (Model_SPM.hpp:128,138), so
 *     A[dom]_k * R^2   must equal   -mu_k^2.
 *
 * This oracle tests the OPERATOR alone — independent of any trajectory, initial condition, or
 * output map (B, C, D). It is therefore different mathematics from both the parity harness
 * (PLAN.md §5.2) and the Carslaw & Jaeger transient series (P1-G3 part b, separate file). It
 * also structurally forbids a spurious second near-zero eigenvalue (the findZeroEigenvalue
 * risk class, audit C3): a mislabelled zero mode would blow the spectral match wide open.
 *
 * ===========================================================================================
 * REGISTERED BANDS (written BEFORE first run, CLAUDE.md §3):
 *   R1 root self-check : |sin(mu_k) - mu_k*cos(mu_k)| < 1e-12   (validates the root finder).
 *   R2 mass mode       : A[dom](model.zero) == 0.0 EXACTLY, both electrodes.
 *   R3 reality/sign    : every nonzero eigenvalue strictly < 0.
 *   R4 spectral match  : for k = 1 .. ceil(nch/2), rel err |lambda_k*R^2 + mu_k^2|/mu_k^2 < 1e-10
 *                        at nch in {5, 8, 12}. [Band ASSUMED per audit §3; the FULL per-mode
 *                        error table is printed so the true resolved-mode count is recorded,
 *                        not guessed. If a high mode misses 1e-10 the run FALSIFIES the assumed
 *                        band -- a deliverable -- and the band is reset to the measured drift.]
 *
 *   OUTCOME 2026-07-09: R1/R2/R3 CONFIRMED. **R4 (1e-10 @ k<=ceil(nch/2)) FALSIFIED** -- the
 *   operator IS exactly the tan(mu)=mu spectrum (mu_1 matches to 3.6e-14 at nch=12), but the
 *   accuracy is the honest Chebyshev spectral-convergence curve, not a flat 1e-10. Measured
 *   fundamental (mu_1) rel err: nch=5 -> 1.98e-5, nch=8 -> 2.25e-10, nch=12 -> 3.63e-14
 *   (~1.5 digits per added node). Modes resolved to rel < 1e-3: 1 (nch=5), 3 (nch=8), 6 (nch=12)
 *   -- i.e. ~ceil(nch/2) modes, so Fable's MODE COUNT was right but the TOLERANCE was ~1e-3,
 *   not 1e-10. Notable [confirmed]: at the production default nch=5, even the fundamental
 *   diffusion eigenvalue is only accurate to ~2e-5 (fine for voltage, now quantified).
 *   Bands RE-REGISTERED from this run (with ~3x margin), assertions below; the ASSUMED band is
 *   recorded FALSIFIED here per CLAUDE.md §3, not silently loosened. PLAN.md P1-G3 updated.
 * ===========================================================================================
 */

#include "../../src/cells/Cell_SPM/Model_SPM.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

namespace {

template <int N>
struct NchTag { static constexpr int value = N; };

//!< nch values registered for P1-G3 (audit §3).
using NchValues = std::tuple<NchTag<5>, NchTag<8>, NchTag<12>>;

/**
 * k-th positive root of tan(mu) = mu, i.e. the root of g(mu) = sin(mu) - mu*cos(mu) = 0 in
 * (k*pi, (k+1/2)*pi). Bracketed bisection (robust, cannot leave the interval) then Newton
 * polish (g'(mu) = mu*sin(mu)) to machine precision. No hardcoded constants — derived, so the
 * test carries its own oracle (CLAUDE.md §2) and R1 checks it.
 */
double mu_root(int k)
{
  const auto g = [](double m) { return std::sin(m) - m * std::cos(m); };
  double a = k * std::numbers::pi;
  double b = (k + 0.5) * std::numbers::pi;
  double ga = g(a);
  for (int it = 0; it < 100; ++it) {
    const double m = 0.5 * (a + b);
    const double gm = g(m);
    if (ga * gm <= 0.0) { b = m; } else { a = m; ga = gm; }
  }
  double m = 0.5 * (a + b);
  for (int it = 0; it < 6; ++it) {
    const double gp = m * std::sin(m);
    if (std::abs(gp) < 1e-30) break;
    m -= g(m) / gp;
  }
  return m;
}

//!< Fundamental (mu_1) rel-err band per nch, RE-REGISTERED 2026-07-09 from the falsification run
//!< (measured 1.98e-5 / 2.25e-10 / 3.63e-14 at nch=5/8/12) with ~3x margin.
constexpr double fundamental_band(int nch)
{
  if (nch <= 5) return 5e-5;
  if (nch <= 8) return 1e-9;
  return 1e-12; // nch >= 12
}

//!< Floor on the number of nonzero modes resolved to rel < 1e-3 (measured 1/3/6 at nch=5/8/12);
//!< a regression that breaks the discretisation drops this count. This is Fable's ~ceil(nch/2)
//!< usable-mode intuition at the CORRECTED tolerance (1e-3, not 1e-10).
constexpr int resolved_floor(int nch)
{
  if (nch <= 5) return 1;
  if (nch <= 8) return 3;
  return 6; // nch >= 12
}

//!< Fundamental-mode rel err |min|lambda|*R^2 - mu_1^2| / mu_1^2 for one Model_SPM<nch> (pos).
template <int nch>
double fundamental_relerr()
{
  slide::Model_SPM<nch> model;
  const double R2 = model.Rp * model.Rp;
  double smallest = 1e300;
  for (int i = 0; i < nch; ++i) {
    if (i == model.zero) continue;
    smallest = std::min(smallest, std::abs(model.A[slide::pos](i)) * R2);
  }
  const double mu1 = mu_root(1);
  return std::abs(smallest - mu1 * mu1) / (mu1 * mu1);
}

} // namespace

TEMPLATE_LIST_TEST_CASE("Chebyshev operator eigenvalues == analytic roots of tan(mu)=mu",
                        "[core][Chebyshev][P1-G3]", NchValues)
{
  constexpr int nch = TestType::value;
  slide::Model_SPM<nch> model;

  for (int dom = 0; dom < 2; ++dom) {
    const double R = (dom == slide::pos) ? model.Rp : model.Rn;
    const double R2 = R * R;

    // R2 mass mode + R3 reality/sign: collect the nonzero eigenvalues (dimensionless: * R^2),
    // magnitude-sorted so index k-1 corresponds to mu_k.
    std::vector<double> lam;
    lam.reserve(nch - 1);
    for (int i = 0; i < nch; ++i) {
      if (i == model.zero) {
        REQUIRE(model.A[dom](i) == 0.0); // R2
        continue;
      }
      CAPTURE(nch, dom, i, model.A[dom](i));
      REQUIRE(model.A[dom](i) < 0.0);    // R3
      lam.push_back(model.A[dom](i) * R2);
    }
    std::sort(lam.begin(), lam.end(),
              [](double x, double y) { return std::abs(x) < std::abs(y); });

    // Full convergence table printed (ALL nonzero modes) for the record; the two evidence-based
    // bands (R4a fundamental, R4b resolved-mode count) are asserted after the loop.
    double fundamental = 1e300;
    int resolved = 0;
    for (int k = 1; k <= nch - 1; ++k) {
      const double mu = mu_root(k);
      REQUIRE(std::abs(std::sin(mu) - mu * std::cos(mu)) < 1e-12); // R1 (fatal: solver self-check)

      const double expected = -mu * mu;
      const double got = lam[static_cast<std::size_t>(k - 1)];
      const double relerr = std::abs(got - expected) / std::abs(expected);
      if (k == 1) fundamental = relerr;
      if (relerr < 1e-3) ++resolved;

      std::printf("Cheb-eig nch=%2d dom=%d k=%d  mu=%.12f  lambda*R^2=% .12g  -mu^2=% .12g  rel=%.3e\n",
                  nch, dom, k, mu, got, expected, relerr);
    }

    CAPTURE(nch, dom, fundamental, resolved);
    REQUIRE(fundamental < fundamental_band(nch));  // R4a: fundamental converges spectrally in nch
    REQUIRE(resolved >= resolved_floor(nch));      // R4b: ~ceil(nch/2) modes usable to rel < 1e-3
  }
}

//!< R4c (structural, nch-independent): the fundamental converges MONOTONICALLY with nch — the
//!< signature of spectral convergence. Catches a regression that breaks convergence even if the
//!< absolute per-nch bands were to drift. Cross-nch comparison needs all three instantiated here.
TEST_CASE("Chebyshev fundamental eigenvalue converges spectrally with nch", "[core][Chebyshev][P1-G3]")
{
  const double e5 = fundamental_relerr<5>();
  const double e8 = fundamental_relerr<8>();
  const double e12 = fundamental_relerr<12>();
  std::printf("Cheb-eig fundamental convergence: nch5=%.3e  nch8=%.3e  nch12=%.3e\n", e5, e8, e12);

  CAPTURE(e5, e8, e12);
  REQUIRE(e5 > e8);        // strictly decreasing = spectral convergence
  REQUIRE(e8 > e12);
  REQUIRE(e12 < 1e-12);    // resolved to ~machine precision by nch=12
}
