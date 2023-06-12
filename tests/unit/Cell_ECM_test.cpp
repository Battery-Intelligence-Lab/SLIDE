/*
 * Cell_ECM_test.cpp
 *
 * Created on: 17 Dec 2019
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../../src/slide.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <span>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using namespace slide;

TEST_CASE("Test ECM constructors", "[ECM]")
{
  Cell_ECM c1;
  REQUIRE_THAT(c1.Cap(), WithinAbs(16, 1e-15));
  REQUIRE_THAT(c1.Vmin(), WithinAbs(2.7, 1e-15));
  REQUIRE_THAT(c1.Vmax(), WithinAbs(4.2, 1e-15));
  REQUIRE_THAT(c1.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.getIr(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5, 1e-15));
  REQUIRE_THAT(c1.T(), WithinAbs(settings::T_ENV, 1e-15));

  double soc = 1;
  double cap = 5;
  Cell_ECM c2(cap, soc);
  REQUIRE_THAT(c2.Cap(), WithinAbs(5, 1e-15));
  REQUIRE_THAT(c2.Vmin(), WithinAbs(2.7, 1e-15));
  REQUIRE_THAT(c2.Vmax(), WithinAbs(4.2, 1e-15));
  REQUIRE_THAT(c2.I(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c2.getIr(), WithinAbs(0, 1e-15));
  REQUIRE_THAT(c2.SOC(), WithinAbs(1, 1e-15));
  REQUIRE_THAT(c2.T(), WithinAbs(settings::T_ENV, 1e-15));
}

TEST_CASE("Test getting ECM states", "[ECM]")
{
  Cell_ECM c1;
  std::vector<double> s;

  c1.getStates(s);
  REQUIRE_THAT(s[State_ECM<1>::i_SOC], WithinAbs(0.5, 1e-15));           //!< soc
  REQUIRE_THAT(s[State_ECM<1>::i_Ir], WithinAbs(0, 1e-15));              //!< Ir
  REQUIRE_THAT(s[State_ECM<1>::i_T], WithinAbs(settings::T_ENV, 1e-15)); //!< T
  REQUIRE_THAT(s[State_ECM<1>::i_I], WithinAbs(0, 1e-15));               //!< current
}

TEST_CASE("Test ECM getV", "[ECM]")
{
  //!< double V(bool print = true); //!< crit is an optional argument
  Cell_ECM c1;

  //!< normal cell, should give no errors
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));
  REQUIRE_THAT(c1.V(), WithinAbs(3.15, 1e-15));

  //!< set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  REQUIRE(V > 3.15);
  c1.timeStep_CC(5);
  REQUIRE(c1.V() > V);

  //!< set to discharge
  V = c1.V();
  c1.setCurrent(1);
  REQUIRE(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  REQUIRE(c1.V() < V);
}

TEST_CASE("Test setting ECM states", "[ECM]")
{
  Cell_ECM c1;

  //!< set valid new states
  double soc{ 0.75 }, ir{ 1 }, i{ 2 }, t{ 273 };
  std::vector<double> s(7); // 4 + 3 (cumulative)
  s[State_ECM<1>::i_SOC] = soc;
  s[State_ECM<1>::i_Ir] = ir;
  s[State_ECM<1>::i_T] = t;
  s[State_ECM<1>::i_I] = i;

  std::span<double> spn{ s };
  c1.setStates(spn, true, true);

  s.clear();
  c1.getStates(s);
  REQUIRE_THAT(s[State_ECM<1>::i_SOC], WithinAbs(soc, 1e-15)); //!< soc
  REQUIRE_THAT(s[State_ECM<1>::i_Ir], WithinAbs(ir, 1e-15));   //!< Ir
  REQUIRE_THAT(s[State_ECM<1>::i_T], WithinAbs(t, 1e-15));     //!< T
  REQUIRE_THAT(s[State_ECM<1>::i_I], WithinAbs(i, 1e-15));     //!< current
}

TEST_CASE("Test ECM validStates", "[ECM]")
{
  Cell_ECM c1;

  //!< set valid new states
  double soc{ 1 }, ir{ 1 }, i{ 2 }, t{ 273 };
  double s[4 + 3]; // 4 + 3 (cumulative)
  s[State_ECM<1>::i_SOC] = soc;
  s[State_ECM<1>::i_Ir] = ir;
  s[State_ECM<1>::i_T] = t;
  s[State_ECM<1>::i_I] = i;
  std::span<double> spn{ s };
  c1.setStates(spn);
  REQUIRE(c1.validStates());
}

TEST_CASE("Test ECM timeStep_CC", "[ECM]")
{
  constexpr double I = -1;
  constexpr double dt = 5;
  constexpr double tol = 0.002;

  Cell_ECM c1;
  //!< soc initial = 0.5 and capacity = 10
  //!< so SOC_end = 0.5 + 1*5/3600 = 0.5014
  c1.setCurrent(I);
  c1.timeStep_CC(dt);
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.5014, tol));

  c1.setCurrent(-I);
  c1.timeStep_CC(dt);
  REQUIRE_THAT(c1.SOC(), WithinAbs(0.50, tol));
}