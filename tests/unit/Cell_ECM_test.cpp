/*
 * Cell_ECM_test.cpp
 *
 * Created on: 17 Dec 2019
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <span>

namespace slide::tests::unit {

bool test_constructor_ECM()
{
  Cell_ECM c1;
  assert(NEAR(c1.Cap(), 16));
  assert(NEAR(c1.Vmin(), 2.7));
  assert(NEAR(c1.Vmax(), 4.2));
  assert(NEAR(c1.I(), 0));
  assert(NEAR(c1.getIr(), 0));
  assert(NEAR(c1.SOC(), 0.5));
  assert(NEAR(c1.T(), settings::T_ENV));

  double soc = 1;
  double cap = 5;
  Cell_ECM c2(cap, soc);
  assert(NEAR(c2.Cap(), 5));
  assert(NEAR(c2.Vmin(), 2.7));
  assert(NEAR(c2.Vmax(), 4.2));
  assert(NEAR(c2.I(), 0));
  assert(NEAR(c2.getIr(), 0));
  assert(NEAR(c2.SOC(), 1));
  assert(NEAR(c2.T(), settings::T_ENV));

  return true;
}

bool test_getStates_ECM()
{
  Cell_ECM c1;
  std::vector<double> s;

  c1.getStates(s);
  assert(NEAR(s[State_ECM<1>::i_SOC], 0.5));           //!< soc
  assert(NEAR(s[State_ECM<1>::i_Ir], 0));              //!< Ir
  assert(NEAR(s[State_ECM<1>::i_T], settings::T_ENV)); //!< T
  assert(NEAR(s[State_ECM<1>::i_I], 0));               //!< current

  std::span<double> spn{ s };
  c1.setStates(spn); //!< #TODO this must throw an error

  return true;
}

bool test_getV_ECM()
{
  //!< double V(bool print = true); //!< crit is an optional argument
  Cell_ECM c1;

  //!< normal cell, should give no errors
  assert(NEAR(c1.V(), 3.15));
  assert(NEAR(c1.V(), 3.15));
  assert(NEAR(c1.V(), 3.15));

  //!< set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  assert(V > 3.15);
  c1.timeStep_CC(5);
  assert(c1.V() > V);

  //!< set to discharge
  V = c1.V();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() < V);

  // //!< cell with SOC out of range #TODO
  // Cell_ECM c2(1, 1); //!< make a cell with soC equal to 1
  // c2.setCurrent(-1, false, false);
  // c2.timeStep_CC(3600); //!< charge further for one hour, now the SOC should be close to 2
  // c2.V(false);

  return true;
}

bool test_setStates_ECM()
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
  assert(NEAR(s[State_ECM<1>::i_SOC], soc)); //!< soc
  assert(NEAR(s[State_ECM<1>::i_Ir], ir));   //!< Ir
  assert(NEAR(s[State_ECM<1>::i_T], t));     //!< T
  assert(NEAR(s[State_ECM<1>::i_I], i));     //!< current

  // //!< set invalid states #TODO
  // soc = 2;
  // t = 0;
  // s[State_ECM::i_SOC] = soc;
  // s[State_ECM::i_T] = t;
  // c1.setStates(spn);

  // //!< set states which violate voltage
  // soc = 1;
  // ir = -5;
  // t = 25_degC;
  // i = -1;
  // s[State_ECM::i_SOC] = soc;
  // s[State_ECM::i_Ir] = ir;
  // s[State_ECM::i_T] = t;
  // s[State_ECM::i_I] = i;

  // c1.setStates(spn);

  return true;
}

bool test_validStates_ECM()
{
  //!< bool validStates(double s[], int nin);
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
  assert(c1.validStates());

  //!< set invalid states
  // s[State_ECM::i_SOC] = 2; //!< soc
  // c1.setStates(spn);

  // assert(!c1.validStates());

  // c1.setStates(spn);
  // s[State_ECM::i_SOC] = 0.5; //!< soc
  // s[State_ECM::i_T] = 0;     //!< T
  // assert(!c1.validStates());

  return true;
}

bool test_timeStep_CC_ECM()
{
  //!< bool timeStep_CC(double dt);
  Cell_ECM c1;
  double I = -1;
  double dt = 5;
  double tol = 0.002;

  //!< soc initial = 0.5 and capacity = 10
  //!< so SOC_end = 0.5 + 1*5/3600 = 0.5014
  c1.setCurrent(I);
  c1.timeStep_CC(dt);
  double err = c1.SOC() - 0.5014;
  assert(err < tol && err > -tol);

  c1.setCurrent(-I);
  c1.timeStep_CC(dt);
  err = c1.SOC() - 0.50;
  assert(err < tol && err > -tol);

  return true;
}

int test_all_Cell_ECM()
{
  /*
   * calls all test-functions
   *
   * IN
   *  	if true, we will also verify that things which have to go wrong, actually do go wrong
   * 					this will result in error messages being printed to the terminal
   * 					all the errors are caught, so the code should not crash
   * 				if false, we only test things which should go well
   */

  //!< if we test the errors, suppress error messages
  if (!TEST(test_constructor_ECM, "test_constructor_ECM")) return 1;
  if (!TEST(test_getStates_ECM, "test_getStates_ECM")) return 2;
  if (!TEST(test_getV_ECM, "test_getV_ECM")) return 3;
  if (!TEST(test_setStates_ECM, "test_setStates_ECM")) return 4;
  if (!TEST(test_validStates_ECM, "test_validStates_ECM")) return 5;
  if (!TEST(test_timeStep_CC_ECM, "test_timeStep_CC_ECM")) return 6;

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Cell_ECM(); }