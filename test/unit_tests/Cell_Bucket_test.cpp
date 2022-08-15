/*
 * Cell_test.cpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "unit_tests.hpp"
#include "../cells/Cell_Bucket/Cell_Bucket.hpp"

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <span>

namespace slide::unit_tests {
void Cell_test()
{
  //!< test the constructors
  Cell_Bucket c1;
  assert(c1.Cap() == 16);
  assert(c1.Vmin() == 2.7);
  assert(c1.Vmax() == 4.2);
  assert(c1.I() == 0);
  assert(c1.SOC() == 0.5);
  assert(c1.T() == settings::T_ENV);

  double soc = 1;
  double cap = 5;
  std::string n = "na";

  Cell_Bucket c2(n, cap, soc);
  assert(c2.Cap() == 5);
  assert(c1.Vmin() == 2.7);
  assert(c1.Vmax() == 4.2);
  assert(c1.I() == 0);
  assert(c2.SOC() == 1);
  assert(c1.T() == settings::T_ENV);
}

void getStates_test(bool fault)
{
  Cell_Bucket c1;
  std::vector<double> s;

  c1.getStates(s);
  assert(s[0] == 0.5);             //!< soc
  assert(s[1] == settings::T_ENV); //!< T
  assert(s[2] == 0);               //!< current

  if (fault) {
    try {
      //!< cout<<"There must be an error message after this line"<<endl<<flush; 	//!< changed global verbose variable
      //!< Something with error.
    } catch (...) {
    }
  }
}

void getV_test(bool fault)
{
  Cell_Bucket c1;
  std::string n = "na";

  //!< normal cell, should give no errors
  assert(c1.V(false) == 3.2);
  assert(c1.V(false) == 3.2);
  assert(c1.V(true) == 3.2);
  assert(c1.V() == 3.2);
  double v;
  int val = c1.checkVoltage(v, false);
  assert(val == 0);
  assert(v == 3.2);

  //!< cell with SOC out of range
  Cell_Bucket c2(n, 1, 1); //!< make a cell with soC equal to 1
  c2.setCurrent(-1, false, false);
  c2.timeStep_CC(3600); //!< charge further for one hour, now the SOC should be close to 2
  if (fault) {
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

void getParent_test()
{
  Cell_Bucket c;

  assert(c.getParent() == nullptr);
}

void setI_test()
{
  Cell_Bucket c1;

  //!< set I without checking the voltage is valid
  assert(c1.I() == 0);
  double V = c1.setCurrent(1.0, false, false);
  assert(c1.I() == 1.0);
  assert(V == 0);

  //!< setCurrent with a valid voltage
  V = c1.setCurrent(0, true, true);
  assert(c1.I() == 0);
  assert(V == 3.2);
  V = c1.setCurrent(0); //!< without optional arguments
  assert(c1.I() == 0);
  assert(V == 3.2);
}
void setT_test()
{
  Cell_Bucket c1;
  assert(c1.T() == settings::T_ENV);
  c1.setT(273);
  assert(c1.T() == 273);
}
void setSOC_test()
{
  Cell_Bucket c1;

  //!< set I without checking the voltage is valid
  assert(c1.SOC() == 0.5);
  auto V = c1.setSOC(1, false, false);
  assert(c1.I() == 1.0);
  assert(V == 0);

  //!< setCurrent with a valid voltage
  V = c1.setCurrent(0, true, true);
  assert(c1.I() == 0);
  assert(V == 3.2);
  V = c1.setCurrent(0); //!< without optional arguments
  assert(c1.I() == 0);
  assert(V == 3.2);
}
void setStates_test(bool fault)
{
  Cell_Bucket c1;
  int nin = settings::CELL_NSTATE_MAX;
  int n;

  //!< set valid new states
  double soc, i, t;
  soc = 0.7;
  i = 2;
  t = 273;

  std::vector<double> s{ soc, t, i };
  std::span<const double> spn{ s };
  c1.setStates(spn, true, true);

  s.clear();
  c1.getStates(s);
  assert(s[0] == soc); //!< soc
  assert(s[1] == t);   //!< T
  assert(s[2] == i);   //!< current

  //!< set invalid states
  if (fault) {
    soc = 2;
    t = 0;
    s[0] = soc;
    s[1] = t;
    try {
      //!< cout<<"There must be three error messages, one about invalid SOC and one invalid T, and one illegal state"<<endl<<flush;	//!< changed global verbose variable
      c1.setStates(spn);
      assert(false);
    } catch (...) {
    };
  }

  //!< set states which violate voltage
  if (fault) {
    soc = 1;
    t = 273 + 25;
    i = -1;
    s[0] = soc;
    s[1] = t;
    s[2] = i;
    try {
      //!< cout<<"There must be one error message about an error when getting the voltage"<<endl<<flush;	//!< changed global verbose variable
      c1.setStates(spn);
      assert(false);
    } catch (...) {
    };
  }
}
void timeStep_CC_test()
{
  Cell_Bucket c1;
  double I = -1;
  double dt = 5;
  double tol = 0.002;

  //!< soc initial = 0.5 and capacity = 10
  //!< so SOC_end = 0.5 + 1*5/3600/capacity =
  c1.setCurrent(I);
  c1.timeStep_CC(dt);
  double err = c1.SOC() - (0.5 - I * dt / 3600.0 / c1.Cap());
  assert(err < tol && err > -tol);

  c1.setCurrent(-I);
  c1.timeStep_CC(dt);
  err = c1.SOC() - 0.50;
  assert(err < tol && err > -tol);
}

void data_test()
{
  slide::Cell_Bucket c1;

  //!< slowly charge while storing data
  double I = -1;
  double dt = 5;
  //!< double tol = 0.002;
  int nstep = 100;
  c1.setCurrent(I);

  //!< store a data point every time step
  for (int i = 0; i < nstep; i++) {
    c1.timeStep_CC(dt, true);
    c1.storeData();
  }

  //!< write the data
  std::string pref = "CellTest_";
  //!< c1.writeData(pref);
}

void Copy_test()
{
  Cell_Bucket c1;

  std::unique_ptr<StorageUnit> c2{ c1.copy() };
  //!< dynamic cast c2 back to a cell
  Cell_Bucket *c22 = dynamic_cast<Cell_Bucket *>(c2.get()); //!< Dynamic cast from smart pointer of StorageUnit to regular pointer of Cell
  assert(c1.SOC() == c22->SOC());
  assert(c1.I() == c22->I());

  //!< now try with a changed SOC
  double socnew = 0.4;
  c1.setSOC(socnew, false, false);
  c2.reset(c1.copy());
  //!< dynamic cast c2 back to a cell
  c22 = dynamic_cast<Cell_Bucket *>(c2.get()); //!< Dynamic cast from StorageUnit to Cell
  assert(c1.SOC() == c22->SOC());
  assert(c1.I() == c22->I());

  //!< now change c22 and double check c1 hasn't changed
  c22->setSOC(0.6, false, false);
  assert(c1.SOC() == socnew);
}

void test_Cell_Bucket(bool testErrors)
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
  Cell_test();
  getStates_test(testErrors);
  getV_test(testErrors);
  getParent_test();
  setI_test();
  setStates_test(testErrors);
  timeStep_CC_test();

  data_test();
  Copy_test();

  std::cout << "test_Cell_Bucket is completed successfully!\n";
}
} // namespace slide::unit_tests
