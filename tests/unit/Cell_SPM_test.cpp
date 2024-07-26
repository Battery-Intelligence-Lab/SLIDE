/*
 * Cell_SPM_test.cpp
 *
 *  Created on: 8 Feb 2020
 *   Author(s): Jorn Reniers
 */

#include "../../src/slide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <span>

using Catch::Matchers::WithinAbs;
using namespace slide;

constexpr double TOL_EQ = 1e-15;

TEST_CASE("test_constructor_SPM", "[CELL_SPM]")
{
  Cell_SPM c1;
  REQUIRE_THAT(c1.Cap(), WithinAbs(16, TOL_EQ));
  REQUIRE_THAT(c1.Vmin(), WithinAbs(2.7, TOL_EQ));
  REQUIRE_THAT(c1.Vmax(), WithinAbs(4.2, TOL_EQ));
  REQUIRE_THAT(c1.I(), WithinAbs(0, TOL_EQ));
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5, TOL_EQ));
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, TOL_EQ));
}

TEST_CASE("test_getStates_SPM", "[CELL_SPM]")
{
  Cell_SPM c1;

  constexpr double tol = 0.001;
  constexpr double tol2 = tol / 100.0;
  auto s = c1.getStateObj();

  // Check if zp and zn values are as expected.
  double expectedZp[] = { 0, 0, 0, 0.47605127273, 0 };
  double expectedZn[] = { 0, 0, 0, 0.289437188135, 0 };

  for (int i = 0; i < 5; i++) {
    REQUIRE_THAT(s.zp(i), WithinAbs(expectedZp[i], tol2));
    REQUIRE_THAT(s.zn(i), WithinAbs(expectedZn[i], tol2));
  }

  // Check if remaining states are as expected.
  REQUIRE_THAT(s.delta(), WithinAbs(1e-9, TOL_EQ));
  REQUIRE_THAT(s.LLI(), WithinAbs(0, TOL_EQ));
  REQUIRE_THAT(s.thickp(), WithinAbs(70e-6, TOL_EQ));
  REQUIRE_THAT(s.thickn(), WithinAbs(73.5e-6, TOL_EQ));
  REQUIRE_THAT(s.ep(), WithinAbs(0.5, TOL_EQ));
  REQUIRE_THAT(s.en(), WithinAbs(0.5, TOL_EQ));
  REQUIRE_THAT(s.ap(), WithinAbs(1.5 / (8.5 * 1e-6), TOL_EQ));
  REQUIRE_THAT(s.an(), WithinAbs(1.5 / (1.25 * 1e-5), TOL_EQ));
  REQUIRE_THAT(s.CS(), WithinAbs(0.01 * 1.5 / (1.25 * 1e-5) * 0.62 * 73.5 * 1e-6, TOL_EQ));
  REQUIRE_THAT(s.Dp(), WithinAbs(8e-14, TOL_EQ));
  REQUIRE_THAT(s.Dn(), WithinAbs(7e-14, TOL_EQ));
  REQUIRE_THAT(s.rDCp(), WithinAbs(0.0028, TOL_EQ));
  REQUIRE_THAT(s.rDCn(), WithinAbs(0.0028, TOL_EQ));
  REQUIRE_THAT(s.rDCcc(), WithinAbs(0.0002325, TOL_EQ));
  REQUIRE_THAT(s.delta_pl(), WithinAbs(0, TOL_EQ));
  REQUIRE_THAT(s.SOC(), WithinAbs(0.5, TOL_EQ));
  REQUIRE_THAT(s.T(), WithinAbs(settings::T_ENV, TOL_EQ));
  REQUIRE_THAT(s.I(), WithinAbs(0.0, TOL_EQ));


  // Check states over c1
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5, TOL_EQ));
  REQUIRE_THAT(c1.I(), WithinAbs(0, TOL_EQ));
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, TOL_EQ));

  // Check CSurf values
  double cps{}, cns{};
  c1.getCSurf(cps, cns, false);
  REQUIRE_THAT(cps, WithinAbs(35421.3, 0.1));
  REQUIRE_THAT(cns, WithinAbs(14644.5, 0.1));
  REQUIRE_THAT(c1.getRdc(), WithinAbs(0.001253, tol));
}

TEST_CASE("test_getV_SPM", "[CELL_SPM]")
{
  Cell_SPM c1;
  double tol = 0.01;

  // Initial voltage
  double Vini = 3.68136;
  REQUIRE_THAT(c1.V(), WithinAbs(Vini, tol));
  REQUIRE_THAT(c1.V(), WithinAbs(Vini, tol));

  // Voltage should increase when set to charging
  c1.setCurrent(-1);
  REQUIRE(c1.V() > Vini);

  // Voltage should decrease when set to discharge
  c1.setCurrent(1);
  REQUIRE(c1.V() < Vini);
}


TEST_CASE("test_setStates_SPM", "[CELL_SPM]")
{
  Cell_SPM c1;

  // Set valid new states
  double T = 45_degC, delta = 2e-9, LLI = 1, thickp = 60e-6, thickn = 70e-6,
         ep = 0.4, en = 0.4, ap = 3 * ep / 8.5e-6, an = 3 * en / 1.25e-5,
         CS = 0.1, Dp = 5e-14, Dn = 6e-14, rp = 0.005, rn = 0.005, rcc = 0.002,
         delta_pl = 0.01, SOC = 0.4, I = -1;

  auto &st = c1.getStateObj();

  std::vector<double> sini(st.size()); // Also includes cumulative.

  for (int i = 0; i < State_SPM::nch; i++) {
    sini[i + State_SPM::i_zp] = st.zp(i);
    sini[i + State_SPM::i_zn] = st.zn(i);
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
  c1.setStates(spn, true, true); // This checks states are valid

  // Assertions for all expected states
  REQUIRE_THAT(st.SOC(), WithinAbs(SOC, TOL_EQ));
  REQUIRE_THAT(st.I(), WithinAbs(I, TOL_EQ));
  REQUIRE_THAT(st.T(), WithinAbs(T, TOL_EQ));
  REQUIRE_THAT(st.delta(), WithinAbs(delta, TOL_EQ));
  REQUIRE_THAT(st.CS(), WithinAbs(CS, TOL_EQ));
  REQUIRE_THAT(st.delta_pl(), WithinAbs(delta_pl, TOL_EQ));
  REQUIRE_THAT(st.LLI(), WithinAbs(LLI, TOL_EQ));
  REQUIRE_THAT(st.thickp(), WithinAbs(thickp, TOL_EQ));
  REQUIRE_THAT(st.thickn(), WithinAbs(thickn, TOL_EQ));
  REQUIRE_THAT(st.ep(), WithinAbs(ep, TOL_EQ));
  REQUIRE_THAT(st.en(), WithinAbs(en, TOL_EQ));
  REQUIRE_THAT(st.ap(), WithinAbs(ap, TOL_EQ));
  REQUIRE_THAT(st.an(), WithinAbs(an, TOL_EQ));
  REQUIRE_THAT(st.Dp(), WithinAbs(Dp, TOL_EQ));
  REQUIRE_THAT(st.Dn(), WithinAbs(Dn, TOL_EQ));
  REQUIRE_THAT(st.rDCp(), WithinAbs(rp, TOL_EQ));
  REQUIRE_THAT(st.rDCn(), WithinAbs(rn, TOL_EQ));
  REQUIRE_THAT(st.rDCcc(), WithinAbs(rcc, TOL_EQ));
}

TEST_CASE("test_timeStep_CC_SPM", "[CELL_SPM]")
{
  Cell_SPM c1;
  double soc = c1.SOC();

  // Set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  c1.timeStep_CC(5);
  REQUIRE(c1.V() > V);
  REQUIRE(c1.SOC() > soc);

  // Set to discharge
  V = c1.V();
  soc = c1.SOC();
  c1.setCurrent(1);
  REQUIRE(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  REQUIRE(c1.V() < V);
  REQUIRE(c1.SOC() < soc);
}