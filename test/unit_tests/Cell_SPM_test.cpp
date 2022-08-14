/*
 * CellSPM_test.cpp
 *
 *  Created on: 8 Feb 2020
 *   Author(s): Jorn Reniers
 */

#include "Cell_SPM.hpp"

#include "unit_tests.hpp"
#include "Cell_SPM.hpp"
#include "constants.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

namespace slide::unit_tests {

void test_constructor_SPM()
{
  // Cell_SPM();
  // Cell_SPM(double capin, double SOCin);

  Cell_SPM c1;
  assert(c1.Cap() == 16);
  assert(c1.Vmin() == 2.7);
  assert(c1.Vmax() == 4.2);
  assert(c1.I() == 0);
  assert(c1.SOC() == 0.5);
  assert(c1.T() == settings::T_ENV);
}

void test_getStates_SPM()
{
  // void getStates(double s[], int nin, int&nout);
  Cell_SPM c1;
  int nin = settings::CELL_NSTATE_MAX;
  double s[nin];
  int n;
  double tol = 0.001;
  double tol2 = tol / 100.0;

  c1.getStates(s, nin, n);

  // 0 to nch = zp
  // nch to 2*nch = zn
  assert(std::abs(s[0]) < tol2);
  assert(std::abs(s[1]) < tol2);
  assert(std::abs(s[2]) < tol2);
  assert(std::abs(s[3] - 0.47605127273) < tol2);
  assert(std::abs(s[4]) < tol2);
  assert(std::abs(s[settings::nch + 0]) < tol2);
  assert(std::abs(s[settings::nch + 1]) < tol2);
  assert(std::abs(s[settings::nch + 2]) < tol2);
  assert(std::abs(s[settings::nch + 3] - 0.289437188135) < tol2);
  assert(std::abs(s[settings::nch + 4]) < tol2);

  assert(s[2 * settings::nch + 0] == 1e-9);                                                           // delta
  assert(s[2 * settings::nch + 1] == 0);                                                              // lli
  assert(s[2 * settings::nch + 2] == 70 * 1e-6);                                                      // thickp
  assert(s[2 * settings::nch + 3] == 73.5 * 1e-6);                                                    // thickn
  assert(s[2 * settings::nch + 4] == 0.5);                                                            // ep
  assert(s[2 * settings::nch + 5] == 0.5);                                                            // en
  assert(s[2 * settings::nch + 6] == 1.5 / (8.5 * 1e-6));                                             // ap
  assert(s[2 * settings::nch + 7] == 1.5 / (1.25 * 1e-5));                                            // an
  assert(std::abs(s[2 * settings::nch + 8] - 0.01 * 1.5 / (1.25 * 1e-5) * 0.62 * 73.5 * 1e-6) < tol); // CS
  assert(s[2 * settings::nch + 9] == 8 * 1e-14);                                                      // Dp
  assert(s[2 * settings::nch + 10] == 7 * 1e-14);                                                     // Dn
  assert(s[2 * settings::nch + 11] == 0.0028);                                                        // rp
  assert(s[2 * settings::nch + 12] == 0.0028);                                                        // rn
  assert(s[2 * settings::nch + 13] == 0.0002325);                                                     // rcc
  assert(s[2 * settings::nch + 14] == 0);                                                             // delta_pl
  assert(s[2 * settings::nch + 15] == 0.5);                                                           // SOC
  assert(s[2 * settings::nch + 16] == settings::T_ENV);                                               // T
  assert(s[2 * settings::nch + 17] == 0);                                                             // I

  assert(c1.SOC() == 0.5);
  assert(c1.I() == 0);
  assert(c1.T() == settings::T_ENV);
  assert(std::abs(c1.getCp_surface(c1.getZp()) - 35421.3) < 0.1); // allow slightly larger error since we approximated all elements of the initial array and matrix
  assert(std::abs(c1.getCn_surface(c1.getZn()) - 14644.5) < 0.1); // concentration ~ 10,000 so 0.1 is still a relative error of e-5
  assert(c1.getDelta() == 1e-9);
  assert(std::abs(c1.getCS() - 0.0441 * 31.0 / 25.0) < tol);
  assert(c1.getDelta_pl() == 0);
  assert(c1.getLLI() == 0);
  assert(c1.getThickp() == 70 * 1e-6);
  assert(c1.getThickn() == 73.5 * 1e-6);
  assert(c1.getEp() == 0.5);
  assert(c1.getEn() == 0.5);
  assert(std::abs(c1.getAp() - 1.764705882352941 * 1e5) < tol);
  assert(std::abs(c1.getAn() - 1.2 * 1e5) < tol);
  assert(c1.getDp() == 8 * 1e-14);
  assert(c1.getDn() == 7 * 1e-14);
  assert(c1.getrdcp() == 0.0028);
  assert(c1.getrdcn() == 0.0028);
  assert(c1.getrdccc() == 0.0002325);
  assert(std::abs(c1.getRdc() - 0.001253) < tol);
}

void test_getV_SPM(bool testErrors)
{
  Cell_SPM c1;
  double tol = 0.01;

  // normal cell, should give no errors
  double Vini = 3.68136;
  assert(std::abs(c1.V(false) - Vini) < tol);
  assert(std::abs(c1.V(true) - Vini) < tol);
  assert(std::abs(c1.V() - Vini) < tol);

  // set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  assert(V > Vini);

  // set to discharge
  V = c1.V();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();
}

void test_setStates_SPM(bool testErrors)
{
  Cell_SPM c1;

  // set valid new states
  double zp[settings::nch], zn[settings::nch];
  double T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, rp, rn, rcc, delta_pl, SOC, I;
  T = 273 + 45;
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
  for (int i = 0; i < settings::nch; i++) {
    zp[i] = c1.getZp()[i];
    zn[i] = c1.getZn()[i];
  }

  double sini[2 * settings::nch + 18];
  for (int i = 0; i < settings::nch; i++) {
    sini[i] = zp[i];
    sini[settings::nch + i] = zn[i];
  }
  sini[2 * settings::nch + 0] = delta;
  sini[2 * settings::nch + 1] = LLI;
  sini[2 * settings::nch + 2] = thickp;
  sini[2 * settings::nch + 3] = thickn;
  sini[2 * settings::nch + 4] = ep;
  sini[2 * settings::nch + 5] = en;
  sini[2 * settings::nch + 6] = ap;
  sini[2 * settings::nch + 7] = an;
  sini[2 * settings::nch + 8] = CS;
  sini[2 * settings::nch + 9] = Dp;
  sini[2 * settings::nch + 10] = Dn;
  sini[2 * settings::nch + 11] = rp;
  sini[2 * settings::nch + 12] = rn;
  sini[2 * settings::nch + 13] = rcc;
  sini[2 * settings::nch + 14] = delta_pl;
  sini[2 * settings::nch + 15] = SOC;
  sini[2 * settings::nch + 16] = T;
  sini[2 * settings::nch + 17] = I;
  c1.setStates(sini, 2 * settings::nch + 18, true); // this checks states are valid

  assert(c1.SOC() == SOC);
  assert(c1.I() == I);
  assert(c1.T() == T);
  assert(c1.getDelta() == delta);
  assert(c1.getCS() == CS);
  assert(c1.getDelta_pl() == delta_pl);
  assert(c1.getLLI() == LLI);
  assert(c1.getThickp() == thickp);
  assert(c1.getThickn() == thickn);
  assert(c1.getEp() == ep);
  assert(c1.getEn() == en);
  assert(c1.getAp() == ap);
  assert(c1.getAn() == an);
  assert(c1.getDp() == Dp);
  assert(c1.getDn() == Dn);
  assert(c1.getrdcp() == rp);
  assert(c1.getrdcn() == rn);
  assert(c1.getrdccc() == rcc);

  // set invalid states
  if (testErrors) {

    // test with Ap != 3*e/R
    ap = 3 * ep / 8.5 * pow(10, -6);
    sini[2 * settings::nch + 6] = ap;
    try {
      c1.setStates(sini, 2 * settings::nch + 18, true); // this checks states are valid
      assert(false);
    } catch (...) {
    };
    sini[2 * settings::nch + 6] = 3 * ep / (8.5 * pow(10, -6));

    // test with negative uniform concentration
    zp[3] = -zp[3];
    for (int i = 0; i < settings::nch; i++)
      sini[i] = zp[i];
    try {
      c1.setStates(sini, 2 * settings::nch + 18, true); // this checks states are valid
      assert(false);
    } catch (...) {
    };
  }
}

void test_timeStep_CC_SPM()
{
  // void timeStep_CC(double dt);
  Cell_SPM c1;
  double soc = c1.SOC();

  // set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() > V);
  assert(c1.SOC() > soc);

  // set to discharge
  V = c1.V();
  soc = c1.SOC();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() < V);
  assert(c1.SOC() < soc);
}

void testCell_SPM(bool testErrors)
{
  /*
   * calls all test-functions
   *
   * IN
   * testErrors 	if true, we will also verify that things which have to go wrong, actually do go wrong
   * 					this will result in error messages being printed to the terminal
   * 					all the errors are caught, so the code should not crash
   * 				if false, we only test things which should go well
   */

  // if we test the errors, suppress error messages

  test_constructor_SPM();
  test_getStates_SPM();
  test_getV_SPM(testErrors);
  test_setStates_SPM(testErrors);
  test_timeStep_CC_SPM();
}
} // namespace slide::unit_tests