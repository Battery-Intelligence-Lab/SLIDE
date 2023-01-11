/*
 * Cell_test.cpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <span>

namespace slide::tests::unit {
bool Cell_Bucket_test()
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

  return true;
}

bool getStates_test()
{
  Cell_Bucket c1;
  std::vector<double> s;

  c1.getStates(s);
  assert(s[0] == 0.5);             //!< soc
  assert(s[1] == settings::T_ENV); //!< T
  assert(s[2] == 0);               //!< current
  return true;
}

bool getV_test()
{
  Cell_Bucket c1;
  std::string n = "na";

  //!< normal cell, should give no errors
  assert(c1.V(false) == 3.2);
  assert(c1.V(false) == 3.2);
  assert(c1.V(true) == 3.2);
  assert(c1.V() == 3.2);
  double v;
  auto val = c1.checkVoltage(v, false);
  assert(val == Status::Success);
  assert(v == 3.2);

  //!< cell with SOC out of range
  Cell_Bucket c2(n, 1, 1); //!< make a cell with soC equal to 1
  c2.setCurrent(-1, false, false);
  c2.timeStep_CC(3600); //!< charge further for one hour, now the SOC should be close to 2
  assert(c2.V(true) > 0);
  return true;
}

bool getParent_test()
{
  Cell_Bucket c;
  assert(c.getParent() == nullptr);
  return true;
}

bool setCurrent_test()
{
  Cell_Bucket c1;

  //!< set I without checking the voltage is valid
  assert(c1.I() == 0);
  auto status = c1.setCurrent(1.0, false, false);
  assert(c1.I() == 1.0);
  assert(status != Status::Success);

  //!< setCurrent with a valid voltage
  c1.setCurrent(0, true, true);
  assert(c1.I() == 0);
  assert(c1.V() == 3.2);
  c1.setCurrent(0); //!< without optional arguments
  assert(c1.I() == 0);
  assert(c1.V() == 3.2);
  return true;
}

bool setT_test()
{
  Cell_Bucket c1;
  assert(c1.T() == settings::T_ENV);
  c1.setT(0.0_degC);
  assert(c1.T() == 0.0_degC);
  return true;
}
bool setSOC_test()
{
  Cell_Bucket c1;

  //!< set I without checking the voltage is valid
  assert(c1.SOC() == 0.5);
  auto status = c1.setSOC(1, false, false);
  assert(c1.I() == 1.0);
  assert(status != Status::Success);

  //!< setCurrent with a valid voltage
  c1.setCurrent(0, true, true);
  assert(c1.I() == 0);
  assert(c1.V() == 3.2);
  c1.setCurrent(0); //!< without optional arguments
  assert(c1.I() == 0);
  assert(c1.V() == 3.2);
  return true;
}
bool setStates_test()
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
  soc = 2;
  t = 0;
  s[0] = soc;
  s[1] = t;
  c1.setStates(spn);

  //!< set states which violate voltage
  soc = 1;
  t = 25.0_degC;
  i = -1;
  s[0] = soc;
  s[1] = t;
  s[2] = i;
  c1.setStates(spn);
  return true;
}
bool timeStep_CC_test()
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
  return true;
}

bool data_test()
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
    c1.timeStep_CC(dt);
    c1.storeData();
  }

  //!< write the data
  std::string pref = "CellTest_";
  //!< c1.writeData(pref);
  return true;
}

bool Copy_test()
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
  return true;
}


int test_all_Cell_Bucket()
{
  if (!TEST(Cell_Bucket_test, "Cell_test")) return 1;
  if (!TEST(getStates_test, "getStates_test")) return 2;
  if (!TEST(getV_test, "getV_test")) return 3;
  if (!TEST(getParent_test, "getParent_test")) return 4;
  if (!TEST(setCurrent_test, "setCurrent_test")) return 5;
  if (!TEST(setStates_test, "setStates_test")) return 6;
  if (!TEST(timeStep_CC_test, "timeStep_CC_test")) return 7;
  if (!TEST(data_test, "data_test")) return 8;
  if (!TEST(Copy_test, "Copy_test")) return 9;

  return 0;
}

} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Cell_Bucket(); }
