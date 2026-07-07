/*
 * Module_p_phase0_test.cpp
 *
 * Regression tests for the Phase-0 solver/state bug fixes (PLAN.md §2.4, bugs A1-A7).
 * Each test registers its pass/fail band BEFORE the assertion and is written so that it
 * FAILS on the pre-fix code and PASSES after the fix.
 *
 * Note: the default Cell_ECM `ocv_coefs` evaluate to a nonsensical voltage (~ -55782 V) on
 * this tree (a pre-existing, out-of-scope bug in Cell_ECM.hpp). Tests that need the parallel
 * solver to converge override each cell's OCV via set_ocv_coefs({slope, intercept}) so the
 * cell voltage is a sane linear function of SOC and the tests isolate the solver bugs.
 *
 *  Created: Phase-0 Agent A
 */

#include "../../src/slide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <vector>

using Catch::Matchers::WithinAbs;
using namespace slide;

namespace {

//!< Test-only subclass to reach the protected pack-solver internals.
class Module_p_access : public Module_p
{
public:
  using Module_p::Module_p;                    //!< inherit constructors
  using Module_p::getVall;                     //!< expose for A6
  using Module_p::setCurrent_analytical_impl;  //!< expose for A5
};

//!< Give a Cell_ECM<N> a sane linear OCV = slope*SOC + intercept (overrides the broken default).
template <size_t N>
void setLinearOCV(Cell_ECM<N> *c, double slope = 1.0, double intercept = 3.3)
{
  c->set_ocv_coefs({ slope, intercept });
}

constexpr double T_ENV = settings::T_ENV;

} // namespace

// ---------------------------------------------------------------------------
// A1: Module::setStates rollback must restore children 0..i to THEIR OWN slice
//     of the original states. Pre-fix used SUs[i] for every j and never advanced
//     per-child offsets, so children 0..i-1 kept their (new) states -> corruption.
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A1_setStates_rollback_restores_all_children", "[Module_p][phase0]")
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>(), make<Cell_Bucket>() };
  auto c0 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto c1 = dynamic_cast<Cell_Bucket *>(cs[1].get());
  auto c2 = dynamic_cast<Cell_Bucket *>(cs[2].get());

  //!< Distinct, valid original SOCs so a wrong slice/offset would be observable.
  c0->setSOC(0.50, false);
  c1->setSOC(0.70, false);
  c2->setSOC(0.30, false);

  auto mp = make<Module_p>("A1", T_ENV, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, false, true);

  std::vector<double> sorig;
  mp->getStates(sorig);

  //!< per-cell state stride (State_ECM<0> = {T, SOC, I, time, Ah, Wh}); SOC index = 1.
  std::vector<double> one;
  c0->getStates(one);
  const size_t ns = one.size();
  const size_t iSOC = State_ECM<0>::i_SOC;

  //!< New target: children 0,1 valid-but-changed; child 2 invalid (SOC out of [0,1]) -> bad status.
  std::vector<double> s = sorig;
  s[0 * ns + iSOC] = 0.55;
  s[1 * ns + iSOC] = 0.65;
  s[2 * ns + iSOC] = 2.00; //!< invalid -> triggers isStatusBad on child 2

  //!< REGISTERED (before run): setStates returns a bad status AND children 0,1 are restored
  //!< to their ORIGINAL SOCs (0.50 and 0.70), not the new (0.55 / 0.65).  band |dSOC| < 1e-12.
  int n = 0;
  const Status st = mp->setStates(s, n, false, false);

  REQUIRE(isStatusBad(st));
  REQUIRE_THAT(c0->SOC(), WithinAbs(0.50, 1e-12)); // pre-fix: stuck at 0.55
  REQUIRE_THAT(c1->SOC(), WithinAbs(0.70, 1e-12)); // pre-fix: stuck at 0.65
  REQUIRE_THAT(c2->SOC(), WithinAbs(0.30, 1e-12)); // child that failed restores itself
}

// ---------------------------------------------------------------------------
// A2: solver state must be per-call, not function-static shared across instances.
//     Solve a 4-cell module then a 2-cell module (different child counts). Pre-fix
//     the 2-cell solve reuses the 4-cell static matrices -> wrong currents / no
//     convergence. Post-fix each module solves correctly, used alternately.
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A2_solver_state_not_shared_across_modules", "[Module_p][phase0]")
{
  using RCp = Cell_ECM<0>::R_C_pair;
  std::span<RCp> noRC{};

  Deep_ptr<StorageUnit> big[] = {
    make<Cell_Bucket>(16.0, 0.5, 0.01, noRC), make<Cell_Bucket>(16.0, 0.5, 0.01, noRC),
    make<Cell_Bucket>(16.0, 0.5, 0.01, noRC), make<Cell_Bucket>(16.0, 0.5, 0.01, noRC)
  };
  Deep_ptr<StorageUnit> small[] = { make<Cell_Bucket>(16.0, 0.5, 0.02, noRC),
                                    make<Cell_Bucket>(16.0, 0.5, 0.02, noRC) };
  for (auto &c : big) setLinearOCV(dynamic_cast<Cell_Bucket *>(c.get()));
  for (auto &c : small) setLinearOCV(dynamic_cast<Cell_Bucket *>(c.get()));

  auto mBig = make<Module_p>("big", T_ENV, true, false, std::size(big), 1, 1);
  auto mSmall = make<Module_p>("small", T_ENV, true, false, std::size(small), 1, 1);
  mBig->setSUs(big, false, true);
  mSmall->setSUs(small, false, true);

  constexpr double tol = 1e-6;

  //!< REGISTERED: identical cells -> equal split. mSmall(2 A)->1 A each; mBig(4 A)->1 A each;
  //!< correct when the two DIFFERENTLY-SIZED modules are solved alternately.  band 1e-6 A.
  //!< Solve the SMALL (2-cell) module FIRST, then the BIG (4-cell). Pre-fix, the first call
  //!< sizes the function-static solver vectors/matrix to 2; the second (4-cell) call then
  //!< indexes them out of bounds -> Eigen assertion abort. Post-fix (per-call locals) both
  //!< solve correctly.
  REQUIRE(isStatusSuccessful(mSmall->setCurrent(2.0, false, false)));
  REQUIRE(isStatusSuccessful(mBig->setCurrent(4.0, false, false))); // pre-fix: OOB abort here
  REQUIRE(isStatusSuccessful(mSmall->setCurrent(2.0, false, false)));
  REQUIRE(isStatusSuccessful(mBig->setCurrent(4.0, false, false)));

  for (auto &c : mBig->getSUs()) REQUIRE_THAT(c->I(), WithinAbs(1.0, tol));
  for (auto &c : mSmall->getSUs()) REQUIRE_THAT(c->I(), WithinAbs(1.0, tol));
  REQUIRE_THAT(mSmall->I(), WithinAbs(2.0, tol));
  REQUIRE_THAT(mBig->I(), WithinAbs(4.0, tol));
}

// ---------------------------------------------------------------------------
// A3: heterogeneous 4-branch parallel pack converges to the exact conductance-weighted
//     current split with the (now refactorised) Jacobian.  R = {0.01,0.02,0.04,0.08} ohm.
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A3_heterogeneous_branch_currents_exact", "[Module_p][phase0]")
{
  using RCp = Cell_ECM<0>::R_C_pair;
  std::span<RCp> noRC{};
  const double R[] = { 0.01, 0.02, 0.04, 0.08 };

  Deep_ptr<StorageUnit> cs[] = {
    make<Cell_Bucket>(16.0, 0.5, R[0], noRC), make<Cell_Bucket>(16.0, 0.5, R[1], noRC),
    make<Cell_Bucket>(16.0, 0.5, R[2], noRC), make<Cell_Bucket>(16.0, 0.5, R[3], noRC)
  };
  for (auto &c : cs) setLinearOCV(dynamic_cast<Cell_Bucket *>(c.get()));

  auto mp = make<Module_p>("A3", T_ENV, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, false, true);

  const double Itot = 4.0;

  //!< Hand-derived exact solution (equal OCV, no contact R): equal terminal voltage =>
  //!< I[k] = Itot * (1/R[k]) / sum(1/R).  ROUND: sum(1/R) = 187.5.
  double gsum = 0.0;
  for (double r : R) gsum += 1.0 / r;

  //!< REGISTERED (before run): status Success (=> iterations < maxIteration=50; no blow-up),
  //!< branch currents match exact to 1e-6 A, and sum(I) == Itot.  Expected iters <= 15.
  const Status st = mp->setCurrent(Itot, false, false);
  REQUIRE(isStatusSuccessful(st));

  double isum = 0.0;
  auto &sus = mp->getSUs();
  for (size_t k = 0; k < std::size(cs); k++) {
    const double Iexact = Itot * (1.0 / R[k]) / gsum;
    REQUIRE_THAT(sus[k]->I(), WithinAbs(Iexact, 1e-6));
    isum += sus[k]->I();
  }
  REQUIRE_THAT(isum, WithinAbs(Itot, 1e-6));
}

// ---------------------------------------------------------------------------
// A4: Qcontact must not double-count. Contact resistor i carries sum(I[i..N-1]);
//     pre-fix the accumulator was never reset across i so the heat grew quadratically.
//     The module is nested under a parent so its timeStep_CC does not run/reset the
//     top-level thermal model (which would zero Qcontact).
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A4_Qcontact_no_double_count", "[Module_p][phase0]")
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  for (auto &c : cs) setLinearOCV(dynamic_cast<Cell_Bucket *>(c.get()));

  auto child = make<Module_p>("child", T_ENV, true, false, std::size(cs), 1, 2); // cooltype 2 = open (sub-module)
  child->setSUs(cs, false, true);
  double Rc[] = { 0.02, 0.01 };
  child->setRcontact(Rc);

  //!< Nest into an HVAC parent so child->timeStep_CC() skips the top-level thermal solve
  //!< (which resets therm.Qcontact). The child keeps its accumulated Qcontact.
  Deep_ptr<StorageUnit> kids[] = { std::move(child) };
  auto parent = make<Module_s>("parent", T_ENV, true, false, 2, 1, 1); // cooltype 1 = HVAC top-level
  parent->setSUs(kids, false, true);
  auto childp = dynamic_cast<Module_p *>(parent->getSUs()[0].get());

  //!< Known branch currents I0 = 3 A, I1 = 1 A (set directly, checkV=false).
  childp->getSUs()[0]->setCurrent(3.0, false, false);
  childp->getSUs()[1]->setCurrent(1.0, false, false);

  //!< REGISTERED (before run): with nstep*dt = 1,
  //!<   Qcontact = Rc[0]*(I0+I1)^2 + Rc[1]*(I1)^2 = 0.02*16 + 0.01*1 = 0.33 J.
  //!<   (pre-fix double-counts to 0.02*16 + 0.01*25 = 0.57 J.)  band 1e-9 J.
  const double dt = 1.0;
  const int nstep = 1;
  childp->timeStep_CC(dt, nstep);

  const double expected = 0.02 * 16.0 + 0.01 * 1.0; // 0.33
  REQUIRE_THAT(childp->getThermQcontact(), WithinAbs(expected, 1e-9));
}

// ---------------------------------------------------------------------------
// A5: the analytical parallel path must not dereference an unchecked dynamic_cast.
//     Calling it with a non-ECM (Cell_SPM) child must return a failure Status, not crash.
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A5_analytical_rejects_non_ECM_child", "[Module_p][phase0]")
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_SPM>(), make<Cell_SPM>() };

  Module_p_access mp("A5", T_ENV, true, false, std::size(cs), 1, 1);
  mp.setSUs(cs, false, true);

  //!< REGISTERED (before run): analytical impl returns Status::Invalid_SUs (no crash / no
  //!< nullptr dereference) when children are not Cell_ECM<1>.
  const Status st = mp.setCurrent_analytical_impl(-2.0, false, false);
  REQUIRE(st == Status::Invalid_SUs);
}

// ---------------------------------------------------------------------------
// A6: Module_p::V() must equal getVall's terminal-voltage definition. V() uses
//     v[0] - I()*Rc[0]; getVall()[0] is the terminal voltage seen along branch 0, and
//     all getVall entries agree once the module is solved. This test verifies the two
//     agree (and that V() == getVall()[0] exactly).
// ---------------------------------------------------------------------------
TEST_CASE("phase0_A6_V_consistent_with_getVall", "[Module_p][phase0]")
{
  using RCp = Cell_ECM<0>::R_C_pair;
  std::span<RCp> noRC{};

  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(16.0, 0.5, 0.01, noRC),
                                 make<Cell_Bucket>(16.0, 0.5, 0.02, noRC),
                                 make<Cell_Bucket>(16.0, 0.5, 0.04, noRC) };
  for (auto &c : cs) setLinearOCV(dynamic_cast<Cell_Bucket *>(c.get()));

  Module_p_access mp("A6", T_ENV, true, false, std::size(cs), 1, 1);
  mp.setSUs(cs, false, true);
  double Rc[] = { 0.005, 0.003, 0.002 };
  mp.setRcontact(Rc);

  REQUIRE(isStatusSuccessful(mp.setCurrent(6.0, false, false)));

  std::vector<double> Vall(mp.getNSUs());
  mp.getVall(std::span<double>(Vall.data(), Vall.size()), false);

  //!< REGISTERED (before run): V() == getVall()[0] to 1e-12 V, and all getVall entries
  //!< (the per-branch terminal voltage) agree with V() to the solver tolerance (< 1e-4 V).
  REQUIRE_THAT(mp.V(), WithinAbs(Vall[0], 1e-12));
  for (double vj : Vall)
    REQUIRE_THAT(mp.V(), WithinAbs(vj, 1e-4));
}
