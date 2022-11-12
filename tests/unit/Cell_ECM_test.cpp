/*
 * Cell_ECM_test.cpp
 *
 *  Created on: 17 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "unit_tests.hpp"
#include "../cells/Cell_ECM/Cell_ECM.hpp"
#include "settings.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <span>

namespace slide::unit_tests {

void test_constructor_ECM()
{
  //!< Cell_ECM();
  //!< Cell_ECM(double capin, double SOCin);

  Cell_ECM c1;
  assert(c1.Cap() == 16);
  assert(c1.Vmin() == 2.7);
  assert(c1.Vmax() == 4.2);
  assert(c1.I() == 0);
  assert(c1.getIr() == 0);
  assert(c1.SOC() == 0.5);
  assert(c1.T() == settings::T_ENV);

  double soc = 1;
  double cap = 5;
  Cell_ECM c2(cap, soc);
  assert(c2.Cap() == 5);
  assert(c2.Vmin() == 2.7);
  assert(c2.Vmax() == 4.2);
  assert(c2.I() == 0);
  assert(c2.getIr() == 0);
  assert(c2.SOC() == 1);
  assert(c2.T() == settings::T_ENV);
}

void test_getStates_ECM(bool testErrors)
{
  //!< void getStates(double s[], int nin, int&nout);
  Cell_ECM c1;
  int nin = settings::CELL_NSTATE_MAX;
  std::vector<double> s;
  int n;

  c1.getStates(s);
  assert(s[0] == 0.5);             //!< soc
  assert(s[1] == 0);               //!< Ir
  assert(s[2] == settings::T_ENV); //!< T
  assert(s[3] == 0);               //!< current

  std::span<double> spn{ s };

  if (testErrors) {
    try {
      //!< cout<<"There must be an error message after this line"<<endl<<flush; 	//!< changed global verbose variable
      //!< c1.getStates(s); //!< this must throw an error
      assert(false);
    } catch (...) {
    }
  }
}

void test_getV_ECM(bool testErrors)
{
  //!< double V(bool print = true); //!< crit is an optional argument
  Cell_ECM c1;

  //!< normal cell, should give no errors
  assert(c1.V(false) == 3.2);
  assert(c1.V(true) == 3.2);
  assert(c1.V() == 3.2);

  //!< set to charging and check the voltage has increased
  c1.setCurrent(-1);
  double V = c1.V();
  assert(V > 3.2);
  c1.timeStep_CC(5);
  assert(c1.V() > V);

  //!< set to discharge
  V = c1.V();
  c1.setCurrent(1);
  assert(c1.V() < V);
  V = c1.V();
  c1.timeStep_CC(5);
  assert(c1.V() < V);

  //!< cell with SOC out of range
  Cell_ECM c2(1, 1); //!< make a cell with soC equal to 1
  c2.setCurrent(-1, false, false);
  c2.timeStep_CC(3600); //!< charge further for one hour, now the SOC should be close to 2
  try {
    //!< cout<<"There should be no error message after this line"<<endl<<flush;	//!< changed global verbose variable
    c2.V(false);
    assert(false);
  } catch (...) {
  }

  if (testErrors) {
    try {
      //!< cout<<"There should be an error message after this line"<<endl<<flush;	//!< changed global verbose variable
      c2.V(true);
      assert(false);
    } catch (...) {
    }
    try {
      //!< cout<<"There should be another error message after this line"<<endl<<flush;	//!< changed global verbose variable
      c2.V();
      assert(false);
    } catch (...) {
    }
  }
}

void test_setStates_ECM(bool testErrors)
{
  //!< double setStates(double s[], int nin, bool checkV = true, bool print = true);
  Cell_ECM c1;
  int nin = settings::CELL_NSTATE_MAX;
  int n;

  //!< set valid new states
  double soc, ir, i, t;
  soc = 0.75;
  ir = 1;
  i = 2;
  t = 273;
  std::vector<double> s{ soc, ir, t, i };
  std::span<double> spn{ s };
  c1.setStates(spn, true, true);

  s.clear();
  c1.getStates(s);
  assert(s[0] == soc); //!< soc
  assert(s[1] == ir);  //!< Ir
  assert(s[2] == t);   //!< T
  assert(s[3] == i);   //!< current

  //!< set invalid states
  if (testErrors) {
    soc = 2;
    t = 0;
    s[0] = soc;
    s[2] = t;
    try {
      //!< cout<<"There must be three error messages, one about invalid SOC and one invalid T, and one illegal state"<<endl<<flush;	//!< changed global verbose variable
      c1.setStates(spn);
      assert(false);
    } catch (...) {
    };

    //!< set states which violate voltage
    soc = 1;
    ir = -5;
    t = 273 + 25;
    i = -1;
    s[0] = soc;
    s[1] = ir;
    s[2] = t;
    s[3] = i;
    try {
      //!< cout<<"There must be one error message about an error when getting the voltage"<<endl<<flush;	//!< changed global verbose variable
      c1.setStates(s, n);
      assert(false);
    } catch (...) {
    };
  }
}

void test_validStates_ECM()
{
  //!< bool validStates(double s[], int nin);
  Cell_ECM c1;
  int nin = settings::CELL_NSTATE_MAX;

  //!< set valid new states
  double soc, ir, i, t;
  soc = 1;
  ir = 1;
  i = 2;
  t = 273;
  double s[nin] = { soc, ir, t, i };

  assert(c1.validStates(s, nin));

  //!< set invalid states
  s[0] = 2; //!< soc
  assert(!c1.validStates(s, nin));
  s[0] = 0.5; //!< soc
  s[2] = 0;   //!< T
  assert(!c1.validStates(s, nin));
}

void test_timeStep_CC_ECM()
{
  //!< void timeStep_CC(double dt);
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
}

void testCell_ECM(bool testErrors)
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

  //!< if we test the errors, suppress error messages
  test_constructor_ECM();
  test_getStates_ECM(testErrors);
  test_getV_ECM(testErrors);
  test_setStates_ECM(testErrors);
  test_validStates_ECM();
  test_timeStep_CC_ECM();
}
} // namespace slide::unit_tests