/*
 * Cell_test.cpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../../src/slide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <span>

using Catch::Matchers::WithinAbs;
using namespace slide;

TEST_CASE("Cell_Bucket_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;
  REQUIRE_THAT(c1.Cap(), WithinAbs(16, 1e-15));
  REQUIRE_THAT(c1.Vmin(), WithinAbs(2.7, 1e-15));
  REQUIRE_THAT(c1.Vmax(), WithinAbs(4.2, 1e-15));
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5, 1e-15));
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, 1e-15));

  double soc = 1;
  double cap = 5;
  std::string n = "na";

  Cell_Bucket c2(n, cap, soc);
  REQUIRE_THAT(c2.Cap(), WithinAbs(5, 1e-15));
  REQUIRE_THAT(c1.Vmin(), WithinAbs(2.7, 1e-15));
  REQUIRE_THAT(c1.Vmax(), WithinAbs(4.2, 1e-15));
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c2.SOC(), WithinAbs(1, 1e-15));
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, 1e-15));
}


TEST_CASE("getStates_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;
  std::vector<double> s;

  c1.getStates(s);
  REQUIRE_THAT(s[State_Bucket::i_SOC], WithinAbs(0.5, 1e-15));
  REQUIRE_THAT(s[State_Bucket::i_T], WithinAbs(settings::T_ENV, 1e-15));
  REQUIRE_THAT(s[State_Bucket::i_I], WithinAbs(0, 1e-15));
}

TEST_CASE("getV_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;
  std::string n = "na";
  //!< normal cell, should give no errors
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  double v{};
  auto val = c1.checkVoltage(v, false);
  REQUIRE(isStatusSuccessful(val));
  REQUIRE_THAT(v, WithinAbs(3.15, 1e-15));

  //!< cell with SOC out of range
  Cell_Bucket c2(n, 1, 1); //!< make a cell with soC equal to 1
  c2.setCurrent(-1, false, false);
  c2.timeStep_CC(3600); //!< charge further for one hour, now the SOC should be close to 2
  REQUIRE(c2.V() <= 0);
}

TEST_CASE("getParent_test", "[Cell_Bucket]")
{
  Cell_Bucket c;
  REQUIRE(c.getParent() == nullptr);
}

TEST_CASE("setCurrent_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;

  // set I without checking the voltage is valid
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  auto status = c1.setCurrent(1.0, false, false);
  REQUIRE_THAT(c1.I(), WithinAbs(1.0, 1e-15));
  REQUIRE(isStatusSuccessful(status));

  // setCurrent with a valid voltage
  c1.setCurrent(0, true, true);
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  c1.setCurrent(0); // without optional arguments
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
}

TEST_CASE("setT_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, 1e-15));
  c1.setT(0.0_degC);
  REQUIRE_THAT(c1.T(), WithinAbs(0.0_degC, 1e-15));
}

TEST_CASE("setSOC_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;

  // set I without checking the voltage is valid
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5, 1e-15));
  auto status = c1.setSOC(1, false, false);
  REQUIRE_THAT(c1.SOC(), WithinAbs(1.0, 1e-15));
  REQUIRE(isStatusSuccessful(status));

  // setCurrent with a valid voltage
  c1.setSOC(0.5, false, false);
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  c1.setCurrent(0); // without optional arguments
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
}

TEST_CASE("setStates_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;

  //!< set valid new states
  double soc{ 0.7 }, i{ 2 }, t{ 273 };
  std::vector<double> s(6); // 3 + 3 (cumulative)
  s[State_Bucket::i_SOC] = soc;
  s[State_Bucket::i_T] = t;
  s[State_Bucket::i_I] = i;

  int n = 0;
  c1.setStates(s, n, true, true);

  s.clear();
  c1.getStates(s);

  REQUIRE_THAT(s[State_Bucket::i_SOC], WithinAbs(soc, 1e-15)); //!< soc
  REQUIRE_THAT(s[State_Bucket::i_T], WithinAbs(t, 1e-15));     //!< T
  REQUIRE_THAT(s[State_Bucket::i_I], WithinAbs(i, 1e-15));     //!< current
}

TEST_CASE("timeStep_CC_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;
  double I{ -1 }, dt{ 5 };
  double tol = 0.002;

  //!< soc initial = 0.5 and capacity = 10
  //!< so SOC_end = 0.5 + 1*5/3600/capacity =
  c1.setCurrent(I);
  c1.timeStep_CC(dt);
  double target = 0.5 - I * dt / 3600.0 / c1.Cap();
  REQUIRE_THAT(c1.SOC(), WithinAbs(target, tol));

  c1.setCurrent(-I);
  c1.timeStep_CC(dt);
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.50, tol));
}

TEST_CASE("data_test", "[Cell_Bucket]")
{
  slide::Cell_Bucket c1;
  double I{ -1 }, dt{ 5 };
  int nstep = 100;
  c1.setCurrent(I);

  //!< store a data point every time step
  for (int i = 0; i < nstep; i++) {
    c1.timeStep_CC(dt);
    c1.storeData();
  }

  // The writing of data is commented out
  // std::string pref = "CellTest_";
  // c1.writeData(pref);
  REQUIRE(true); // #TODO see if we can add an actual test.
}


TEST_CASE("Copy_test", "[Cell_Bucket]")
{
  Cell_Bucket c1;

  Deep_ptr<StorageUnit> c2{ c1.copy() };
  Cell_Bucket *c22 = dynamic_cast<Cell_Bucket *>(c2.get()); //!< dynamic cast c2 back to a cell
  REQUIRE_THAT(c1.SOC(), WithinAbs(c22->SOC(), 1e-15));
  REQUIRE_THAT(c1.I(), WithinAbs(c22->I(), 1e-15));

  //!< now try with a changed SOC
  double socnew = 0.4;
  c1.setSOC(socnew, false, false);
  c2.reset(c1.copy());
  //!< dynamic cast c2 back to a cell
  c22 = dynamic_cast<Cell_Bucket *>(c2.get()); //!< Dynamic cast from StorageUnit to Cell
  REQUIRE_THAT(c1.SOC(), WithinAbs(c22->SOC(), 1e-15));
  REQUIRE_THAT(c1.I(), WithinAbs(c22->I(), 1e-15));

  c22->setSOC(0.6, false, false);
  REQUIRE_THAT(c1.SOC(), WithinAbs(socnew, 1e-15));
}