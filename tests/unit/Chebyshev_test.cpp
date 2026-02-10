/**
 * @file Chebyshev_test.cpp
 * @brief Unit tests for the Chebyshev spectral discretisation of the solid diffusion PDE
 *
 * Tests Model_SPM construction and state-space model correctness for multiple nch values,
 * verifying that the discretisation works for arbitrary nch (not just nch=5).
 *
 * Reference: https://github.com/davidhowey/Spectral_li-ion_SPM
 */

#include "../../src/cells/Cell_SPM/Model_SPM.hpp"
#include "../../src/settings/settings.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numeric>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

//!< Tag types to pass nch as a compile-time value to TEMPLATE_LIST_TEST_CASE
template <int N>
struct NchTag
{
  static constexpr int value = N;
};

using NchValues = std::tuple<NchTag<3>, NchTag<5>, NchTag<7>, NchTag<10>>;

//!< Physical constants used in the round-trip tests (Kokam NMC defaults)
constexpr double Rp = 8.5e-6;    //!< cathode particle radius [m]
constexpr double Rn = 12.5e-6;   //!< anode particle radius [m]
constexpr double CmaxP = 51385;  //!< max Li concentration cathode [mol m-3]
constexpr double CmaxN = 30555;  //!< max Li concentration anode [mol m-3]
constexpr double D_ref = 8e-14;  //!< reference diffusion constant [m^2 s-1]

} // namespace

// ---------------------------------------------------------------------------
// Test 1: Chebyshev node positions
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Chebyshev nodes are in (0,1) and decreasing", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  slide::Model_SPM<nch> model;

  REQUIRE(model.xch.size() == nch);

  for (int i = 0; i < nch; ++i) {
    CAPTURE(i, nch, model.xch(i));
    REQUIRE(model.xch(i) > 0.0);
    REQUIRE(model.xch(i) < 1.0);
  }

  //!< Nodes should be strictly decreasing (surface → centre)
  for (int i = 0; i < nch - 1; ++i) {
    CAPTURE(i, model.xch(i), model.xch(i + 1));
    REQUIRE(model.xch(i) > model.xch(i + 1));
  }
}

// ---------------------------------------------------------------------------
// Test 2: Eigenvalues are real and non-positive
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Eigenvalues are real and non-positive", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  slide::Model_SPM<nch> model;

  for (int dom = 0; dom < 2; ++dom) {
    for (int i = 0; i < nch; ++i) {
      CAPTURE(dom, i, nch, model.A[dom](i));
      REQUIRE(model.A[dom](i) <= 0.0);
    }
  }
}

// ---------------------------------------------------------------------------
// Test 3: Exactly one zero eigenvalue per electrode, at the stored index
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Exactly one zero eigenvalue at stored index", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  slide::Model_SPM<nch> model;

  int zero_count_pos = 0, zero_count_neg = 0;
  for (int i = 0; i < nch; ++i) {
    if (model.A[slide::pos](i) == 0.0) zero_count_pos++;
    if (model.A[slide::neg](i) == 0.0) zero_count_neg++;
  }
  REQUIRE(zero_count_pos == 1);
  REQUIRE(zero_count_neg == 1);

  REQUIRE(model.A[slide::pos](model.zero) == 0.0);
  REQUIRE(model.A[slide::neg](model.zero) == 0.0);
}

// ---------------------------------------------------------------------------
// Test 4: Uniform concentration round-trip (setC → getC logic)
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Uniform concentration round-trip", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  constexpr int N = nch + 1;
  slide::Model_SPM<nch> model;

  const double lifrac = 0.5;

  for (int dom = 0; dom < 2; ++dom) {
    const double R = (dom == slide::pos) ? Rp : Rn;
    const double Cmax = (dom == slide::pos) ? CmaxP : CmaxN;
    const double concentration = lifrac * Cmax;

    //!< Simulate setC: compute z for uniform concentration
    std::array<double, nch> z{};
    {
      double zu = 0.0;
      for (int i = 0; i < nch; ++i) {
        const double u = R * concentration * model.xch(i);
        zu += model.V[dom](model.zero, i) * u;
      }
      z[model.zero] = zu;
    }

    //!< Simulate getC: recover concentration at surface + inner nodes (zero current)
    std::array<double, N> c{};
    for (int i = 0; i < N; ++i) {
      double cpt = 0.0;
      for (int j = 0; j < nch; ++j)
        cpt += model.C[dom](i, j) * z[j];
      c[i] = cpt; //!< D term is zero since molarFlux = 0
    }

    //!< All node concentrations should equal the uniform concentration
    for (int i = 0; i < N; ++i) {
      CAPTURE(dom, i, nch, c[i], concentration);
      REQUIRE_THAT(c[i], WithinRel(concentration, 1e-8));
    }
  }
}

// ---------------------------------------------------------------------------
// Test 5: Centre concentration via Cc for uniform concentration
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Centre concentration via Cc for uniform c", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  constexpr int N = nch + 1;
  slide::Model_SPM<nch> model;

  const double lifrac = 0.5;
  const double concentration = lifrac * CmaxP;

  //!< Simulate setC for positive electrode
  std::array<double, nch> z{};
  {
    double zu = 0.0;
    for (int i = 0; i < nch; ++i) {
      const double u = Rp * concentration * model.xch(i);
      zu += model.V[slide::pos](model.zero, i) * u;
    }
    z[model.zero] = zu;
  }

  //!< Recover surface + inner node concentrations (zero current)
  std::array<double, N> c{};
  for (int i = 0; i < N; ++i) {
    double cpt = 0.0;
    for (int j = 0; j < nch; ++j)
      cpt += model.C[slide::pos](i, j) * z[j];
    c[i] = cpt;
  }

  //!< Compute centre concentration using Cc (molarFlux = 0)
  double c_centre_sum = 0.0;
  for (int i = 0; i < N; ++i)
    c_centre_sum += model.Cc(i) * c[i];
  const double c_centre = model.cc_coeff * c_centre_sum; //!< molarFlux*R/Dt = 0

  CAPTURE(nch, c_centre, concentration);
  REQUIRE_THAT(c_centre, WithinRel(concentration, 1e-6));
}

// ---------------------------------------------------------------------------
// Test 6: Zero time-derivative for uniform concentration at zero current
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Uniform concentration has zero time derivative", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  slide::Model_SPM<nch> model;

  const double lifrac = 0.5;
  const double concentration = lifrac * CmaxP;

  //!< Set z to uniform concentration
  std::array<double, nch> z{};
  {
    double zu = 0.0;
    for (int i = 0; i < nch; ++i) {
      const double u = Rp * concentration * model.xch(i);
      zu += model.V[slide::pos](model.zero, i) * u;
    }
    z[model.zero] = zu;
  }

  //!< dz/dt = D * A * z + B * j  (j = 0 for zero current)
  for (int k = 0; k < nch; ++k) {
    const double dzdt = D_ref * model.A[slide::pos](k) * z[k];
    CAPTURE(k, nch, dzdt, model.A[slide::pos](k), z[k]);
    REQUIRE_THAT(dzdt, WithinAbs(0.0, 1e-20));
  }
}

// ---------------------------------------------------------------------------
// Test 7: Q integration matrix first row is zero
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Q matrix first row is zero", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  constexpr int M = 2 * (nch + 1);
  slide::Model_SPM<nch> model;

  for (int j = 0; j <= M; ++j) {
    CAPTURE(j, nch);
    REQUIRE_THAT(model.Q(0, j), WithinAbs(0.0, 1e-14));
  }
}

// ---------------------------------------------------------------------------
// Test 8: Differentiation matrix accuracy on polynomial f(x) = x^2
//
// We test that the first-order differentiation matrix correctly computes
// f'(x) = 2x on the full Chebyshev nodes. Since the differentiation matrix
// is consumed by DN1/DN2 and not stored directly, we instead verify via
// the state-space model: for a known linear concentration profile c(r) = r/R
// (i.e. u(r) = r * c = r^2/R), the derivative dc/dr = 1/R should be recoverable.
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("State-space model recovers linear concentration profile", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  constexpr int N = nch + 1;
  slide::Model_SPM<nch> model;

  //!< Linear concentration c(r) = r / Rp at each inner node
  //!< Transformed: u(r_i) = r_i * c(r_i) = r_i^2 / Rp
  //!< In eigenspace: z = V^{-1} * u
  Eigen::Vector<double, nch> u_vec;
  for (int i = 0; i < nch; ++i) {
    const double r_i = model.xch(i) * Rp;
    u_vec(i) = r_i * r_i / Rp; //!< u_i = r_i^2 / Rp
  }

  //!< Transform to eigenspace: z = V_inv * u (V stores the inverse already)
  Eigen::Vector<double, nch> z_vec = model.V[slide::pos] * u_vec;

  //!< Recover concentration at all nodes: c = C * z (no flux term for this test)
  for (int i = 1; i < N; ++i) { //!< skip surface (row 0 has different formula)
    double c_recovered = 0.0;
    for (int j = 0; j < nch; ++j)
      c_recovered += model.C[slide::pos](i, j) * z_vec(j);

    //!< Expected: c(r_i) = r_i / Rp = xch(i-1), where i-1 maps C matrix row i to xch index
    const double c_expected = model.xch(i - 1); //!< r_i / Rp = xch(i-1)
    CAPTURE(i, nch, c_recovered, c_expected);
    REQUIRE_THAT(c_recovered, WithinRel(c_expected, 1e-6));
  }
}

// ---------------------------------------------------------------------------
// Test 9: nch=5 backward compatibility — setC → getC round-trip with Kokam NMC values
// ---------------------------------------------------------------------------
TEST_CASE("nch=5 backward compatibility", "[Chebyshev]")
{
  constexpr int nch = 5;
  constexpr int N = nch + 1;
  slide::Model_SPM<nch> model;

  //!< Zero eigenvalue should be correctly set
  REQUIRE(model.A[slide::pos](model.zero) == 0.0);
  REQUIRE(model.A[slide::neg](model.zero) == 0.0);

  //!< Kokam NMC at ~50% SOC
  const double xp_now = 0.983999588653496 + 0.5 * (0.400145394039564 - 0.983999588653496);
  const double xn_now = 0.029397569380507 + 0.5 * (0.932469496648387 - 0.029397569380507);
  const double conc_p = xp_now * CmaxP;
  const double conc_n = xn_now * CmaxN;

  //!< For each electrode: setC (uniform) → getC (zero current) must recover the uniform concentration
  for (int dom = 0; dom < 2; ++dom) {
    const double R = (dom == slide::pos) ? Rp : Rn;
    const double conc = (dom == slide::pos) ? conc_p : conc_n;

    //!< Simulate setC: compute z for uniform concentration
    std::array<double, nch> z{};
    {
      double zu = 0.0;
      for (int i = 0; i < nch; ++i) {
        const double u = R * conc * model.xch(i);
        zu += model.V[dom](model.zero, i) * u;
      }
      z[model.zero] = zu;
    }

    //!< Simulate getC: recover concentration at surface + inner nodes (zero current)
    std::array<double, N> c{};
    for (int i = 0; i < N; ++i) {
      double cpt = 0.0;
      for (int j = 0; j < nch; ++j)
        cpt += model.C[dom](i, j) * z[j];
      c[i] = cpt;
    }

    //!< All node concentrations should equal the uniform concentration
    for (int i = 0; i < N; ++i) {
      CAPTURE(dom, i, c[i], conc);
      REQUIRE_THAT(c[i], WithinRel(conc, 1e-8));
    }
  }
}

// ---------------------------------------------------------------------------
// Test 10: Model_SPM matrix dimensions are correct
// ---------------------------------------------------------------------------
TEMPLATE_LIST_TEST_CASE("Model_SPM matrix dimensions", "[Chebyshev]", NchValues)
{
  constexpr int nch = TestType::value;
  constexpr int N = nch + 1;
  constexpr int M = 2 * N;
  slide::Model_SPM<nch> model;

  for (int dom = 0; dom < 2; ++dom) {
    REQUIRE(model.A[dom].size() == nch);
    REQUIRE(model.B[dom].size() == nch);
    REQUIRE(model.C[dom].rows() == N);
    REQUIRE(model.C[dom].cols() == nch);
    REQUIRE(model.D[dom].size() == N);
    REQUIRE(model.V[dom].rows() == nch);
    REQUIRE(model.V[dom].cols() == nch);
  }

  REQUIRE(model.Cc.size() == N);
  REQUIRE(model.Q.rows() == M + 1);
  REQUIRE(model.Q.cols() == M + 1);
}
