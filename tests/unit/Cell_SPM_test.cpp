/*
 * CellSPM_test.cpp
 *
 *  Created on: 8 Feb 2020
 *   Author(s): Jorn Reniers
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <span>

namespace slide::tests::unit {

bool test_constructor_SPM()
{
  Cell_SPM c1;
  assert(NEAR(c1.Cap(), 16));
  assert(NEAR(c1.Vmin(), 2.7));
  assert(NEAR(c1.Vmax(), 4.2));
  assert(NEAR(c1.I(), 0));
  assert(NEAR(c1.SOC(), 0.5));
  assert(NEAR(c1.T(), settings::T_ENV));

  return true;
}

bool test_getStates_SPM()
{
  //!< bool getStates(double s[], int nin, int&nout);
  Cell_SPM c1;
  const double tol = 0.001;
  const double tol2 = tol / 100.0;

  auto s = c1.getStateObj();

  //!< 0 to nch = zp
  //!< nch to 2*nch = zn
  assert(NEAR(s.zp(0), 0, tol2));
  assert(NEAR(s.zp(1), 0, tol2));
  assert(NEAR(s.zp(2), 0, tol2));
  assert(NEAR(s.zp(3), 0.47605127273, tol2));
  assert(NEAR(s.zp(4), 0, tol2));

  assert(NEAR(s.zn(0), 0, tol2));
  assert(NEAR(s.zn(1), 0, tol2));
  assert(NEAR(s.zn(2), 0, tol2));
  assert(NEAR(s.zn(3), 0.289437188135, tol2));
  assert(NEAR(s.zn(4), 0, tol2));

  assert(EQ(s.delta(), 1e-9));
  assert(EQ(s.LLI(), 0));
  assert(EQ(s.thickp(), 70 * 1e-6));
  assert(EQ(s.thickn(), 73.5 * 1e-6));
  assert(EQ(s.ep(), 0.5));
  assert(EQ(s.en(), 0.5));
  assert(EQ(s.ap(), 1.5 / (8.5 * 1e-6)));
  assert(EQ(s.an(), 1.5 / (1.25 * 1e-5)));
  assert(NEAR(s.CS(), 0.01 * 1.5 / (1.25 * 1e-5) * 0.62 * 73.5 * 1e-6, tol));
  assert(EQ(s.Dp(), 8 * 1e-14));
  assert(EQ(s.Dn(), 7 * 1e-14));

  assert(EQ(s.rDCp(), 0.0028));
  assert(EQ(s.rDCn(), 0.0028));
  assert(EQ(s.rDCcc(), 0.0002325));

  assert(EQ(s.delta_pl(), 0.0));
  assert(EQ(s.SOC(), 0.5));
  assert(EQ(s.T(), settings::T_ENV));
  assert(EQ(s.I(), 0.0));

  assert(NEAR(c1.SOC(), 0.5));
  assert(NEAR(c1.I(), 0));
  assert(NEAR(c1.T(), settings::T_ENV));


  double cps{}, cns{};
  c1.getCSurf(cps, cns, false);

  assert(NEAR(cps, 35421.3, 0.1)); //!< allow slightly larger error since we approximated all elements of the initial array and matrix
  assert(NEAR(cns, 14644.5, 0.1)); //!< concentration ~ 10,000 so 0.1 is still a relative error of e-5
  assert(NEAR(c1.getRdc(), 0.001253, tol));

  return true;
}

bool test_getV_SPM()
{
  Cell_SPM c1;
  double tol = 0.01;

  //!< normal cell, should give no errors
  double Vini = 3.68136;
  assert(NEAR(c1.V(), Vini, tol));
  assert(NEAR(c1.V(), Vini, tol));
  assert(NEAR(c1.V(), Vini, tol));

  //!< set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  assert(V > Vini);

  //!< set to discharge
  V = c1.V();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();

  return true;
}

bool test_setStates_SPM()
{
  Cell_SPM c1;

  //!< set valid new states
  double zp[settings::nch], zn[settings::nch];
  double T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, rp, rn, rcc, delta_pl, SOC, I;
  T = 45_degC;
  delta = 2e-9;
  LLI = 1;
  thickp = 60e-6;
  thickn = 70e-6;
  ep = 0.4;
  en = 0.4;
  ap = 3 * ep / 8.5e-6;
  an = 3 * en / 1.25e-5;
  CS = 0.1;
  Dp = 5e-14;
  Dn = 6e-14;
  rp = 0.005;
  rn = 0.005;
  rcc = 0.002;
  delta_pl = 0.01;
  SOC = 0.4;
  I = -1;

  auto &st = c1.getStateObj();

  for (int i = 0; i < State_SPM::nch; i++) {
    zp[i] = st.zp(i);
    zn[i] = st.zn(i);
  }

  std::vector<double> sini(st.size()); // Also includes cumulative.

  for (int i = 0; i < State_SPM::nch; i++) {
    sini[i + State_SPM::i_zp] = zp[i];
    sini[i + State_SPM::i_zn] = zn[i];
  }

  sini[State_SPM::i_delta] = delta;
  sini[State_SPM::i_LLI] = LLI;
  sini[State_SPM::i_thickp] = thickp;
  sini[State_SPM::i_thickn] = thickn;
  sini[State_SPM::i_ep] = ep;
  sini[State_SPM::i_en] = en;
  sini[State_SPM::i_ap] = ap;
  sini[State_SPM::i_an] = an;
  sini[State_SPM::i_CS] = CS;
  sini[State_SPM::i_Dp] = Dp;
  sini[State_SPM::i_Dn] = Dn;
  sini[State_SPM::i_rDCp] = rp;
  sini[State_SPM::i_rDCn] = rn;
  sini[State_SPM::i_rDCcc] = rcc;
  sini[State_SPM::i_delta_pl] = delta_pl;
  sini[State_SPM::i_SOC] = SOC;
  sini[State_SPM::i_T] = T;
  sini[State_SPM::i_I] = I;


  std::span<double> spn(sini);
  c1.setStates(spn, true, true); //!< this checks states are valid

  assert(NEAR(st.SOC(), SOC));
  assert(NEAR(st.I(), I));
  assert(NEAR(st.T(), T));
  assert(NEAR(st.delta(), delta));
  assert(NEAR(st.CS(), CS));
  assert(NEAR(st.delta_pl(), delta_pl));
  assert(NEAR(st.LLI(), LLI));
  assert(NEAR(st.thickp(), thickp));
  assert(NEAR(st.thickn(), thickn));
  assert(NEAR(st.ep(), ep));
  assert(NEAR(st.en(), en));
  assert(NEAR(st.ap(), ap));
  assert(NEAR(st.an(), an));
  assert(NEAR(st.Dp(), Dp));
  assert(NEAR(st.Dn(), Dn));
  assert(NEAR(st.rDCp(), rp));
  assert(NEAR(st.rDCn(), rn));
  assert(NEAR(st.rDCcc(), rcc));


  // //!< set invalid states
  // //!< test with Ap != 3*e/R
  // ap = 3 * ep / 8.5 * 1e-6;
  // sini[State_SPM::i_ap] = ap;
  // auto status = c1.setStates(spn, true, true); //!< this checks states are valid


  // if (isStatusSuccessful(status)) return false; // Failed the test if it does not throw!

  // sini[State_SPM::i_an] = 3 * ep / 8.5e-6;

  // //!< test with negative uniform concentration
  // zp[3] = -zp[3];
  // for (int i = 0; i < settings::nch; i++)
  //   sini[i] = zp[i];


  // status = c1.setStates(spn, true, true);       //!< this checks states are valid
  // if (isStatusSuccessful(status)) return false; // Failed the test if it does not throw!


  return true;
}

bool test_timeStep_CC_SPM()
{
  //!< bool timeStep_CC(double dt);
  Cell_SPM c1;
  double soc = c1.SOC();

  //!< set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() > V);
  assert(c1.SOC() > soc);

  //!< set to discharge
  V = c1.V();
  soc = c1.SOC();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() < V);
  assert(c1.SOC() < soc);

  return true;
}

int test_all_Cell_SPM()
{
  //!< calls all test-functions
  if (!TEST(test_constructor_SPM, "test_constructor_SPM")) return 1;
  if (!TEST(test_getStates_SPM, "test_getStates_SPM")) return 2;
  if (!TEST(test_getV_SPM, "test_getV_SPM")) return 3;
  if (!TEST(test_setStates_SPM, "test_setStates_SPM")) return 4;
  if (!TEST(test_timeStep_CC_SPM, "test_timeStep_CC_SPM")) return 5;

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Cell_SPM(); }