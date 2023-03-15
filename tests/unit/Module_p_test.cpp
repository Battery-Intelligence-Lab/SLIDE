/*
 * Module_base_p_test.cpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cassert>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

namespace slide::tests::unit {

bool test_Constructor_p()
{
  auto mp = make<Module_p>();
  assert(mp->getNSUs() == 0);
  //  assert(mp->T() == settings::T_ENV); #TODO it returns cool->T() which is nullptr.

  auto cp1 = make<Cell_Bucket>();
  auto cp2 = make<Cell_Bucket>();
  assert(cp1->getID() == "Cell_ECM<0>");
  assert(cp1->getFullID() == "Cell_ECM<0>"); //!< has no parent yet

  Deep_ptr<StorageUnit> cs[] = { std::move(cp1), std::move(cp2) };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp2 = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp2->setSUs(cs, checkCells, true);

  assert(mp2->getNSUs() == std::size(cs));
  assert(mp2->T() == T);

  return true;
}

bool test_BasicGetters_p()
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());
  assert(cp1->getID() == "Cell_ECM<0>");
  assert(cp1->getFullID() == "Cell_ECM<0>"); //!< has no parent yet

  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  assert(mp->Cap() == std::size(cs) * cp1->Cap());
  assert(mp->Vmin() == cp1->Vmin());
  assert(mp->Vmax() == cp1->Vmax());
  assert(mp->VMIN() == cp1->VMIN());
  assert(mp->VMAX() == cp1->VMAX());
  assert(mp->I() == 0);
  assert(mp->V() == cp1->V());

  return true;
}
bool test_setI_p()
{
  //!< double Module_base_s::setCurrent(double Inew, bool checkV, bool print)
  double tol = 0.005;
  double Inew, V;

  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());
  assert(cp1->getID() == "Cell_ECM<0>");
  assert(cp1->getFullID() == "Cell_ECM<0>"); //!< has no parent yet


  double v1 = cp1->V();
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  assert(mp->I() == 0);
  assert(mp->V() == cp1->V());

  //!< discharge
  Inew = 1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  assert(NEAR(mp->I(), Inew, tol));
  assert(mp->V() < v1); //!< voltage must decrease
  //!< do not check individual cells, that is done in getCells

  //!< charge
  Inew = -1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  assert(NEAR(mp->I(), Inew, tol));
  assert(mp->V() > v1); //!< voltage must increase

  //!< rest with different SOC values
  Inew = 0;
  cp2->setSOC(0.4);                                              //!< c2 has lower OCV -> should charge
  mp->setCurrent(Inew, true);                                    //!< the large change in OCV causes a large voltage change, which cannot be fixed by setCurrent
  assert(NEAR(cp1->V(), cp2->V(), settings::MODULE_P_V_ABSTOL)); //!< cell voltages are equal
  assert(NEAR(cp1->I(), cp1->I(), settings::MODULE_P_V_ABSTOL)); //!< cell currents are opposite
  assert(cp1->I() > 0);
  assert(cp2->I() < 0);

  //!< test things which should break
  Inew = 10000;                         //!< very large current, should give too low voltage
  Status status = mp->setCurrent(Inew); //!< should fail because the current equation cannot be solved
  if (isStatusSuccessful(status)) return false;


  Inew = -10000;                 //!< very large current, should give too low voltage
  status = mp->setCurrent(Inew); //!< should fail because the current equation cannot be solved
  if (isStatusSuccessful(status)) return false;

  return true;
}

bool test_validStates_p()
{
  //!< bool Module_base_s::validStates(double s[], int nin, bool print)
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());


  assert(cp1->getID() == "Cell_ECM<0>");
  assert(cp1->getFullID() == "Cell_ECM<0>"); //!< has no parent yet
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< valid states (current states)
  std::vector<double> s;
  mp->getStates(s);
  assert(mp->validStates());

  // if () { #TODO failing tests

  //   //!< wrong length
  //   int nc = settings::CELL_NSTATE_MAX;
  //   double sc[nc];
  //   int noutc;
  //   cp1->getStates(sc, nc, noutc);
  //   assert(!mp->validStates(sc, noutc, false));

  //   //!< an SOC which is too large
  //   s[0] = 2; //!< this is the SOC of cell 1
  //   s[nout - 1] = 273 + 1;
  //   assert(!mp->validStates(s, nout, false));

  //   //!< different voltage values (by changing the SOCs)
  //   s[0] = 0.4; //!< this is the SOC of cell 1 (soc of cell 2 is 0.5)
  //   assert(!mp->validStates(s, nout, false));
  // }

  return true;
}

bool test_timeStep_CC_p()
{
  //!< bool Module_base_s::timeStep_CC(double dt)

  double T = settings::T_ENV;
  bool checkCells = false;


  Deep_ptr<StorageUnit> cs[] = {
    Deep_ptr<StorageUnit>(new Cell_Bucket()),
    Deep_ptr<StorageUnit>(new Cell_Bucket())
  };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());


  std::string n = "na";
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  const double v1{ cp1->V() }, soc1{ cp1->SOC() };

  //!< time step with 0 current
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(mp->V() == cp1->V());

  //!< discharge
  double Inew = 1 * std::size(cs);
  double V, err;
  double tol = settings::MODULE_P_I_ABSTOL;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  V = mp->V();
  assert(V < v1);
  err = mp->I() - Inew;
  assert(err < tol && err > -tol); //!< note: this is not going to be exact because we solve a nonlinear equation to obain I

  //!< check individual cells
  auto &cs2 = mp->getSUs();

  for (auto &su : mp->getSUs()) {
    err = su->I() - Inew / mp->getNSUs(); //!< we know the current has to split equally between both cells
    assert(err < tol && err > -tol);
    assert(su->V() < v1);
    auto cell1 = dynamic_cast<Cell_Bucket *>(su.get()); //!< Dynamic cast from StorageUnit to Cell
    assert(cell1->SOC() < soc1);
    err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
    assert(err < tol && err > -tol);
  }

  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  assert(mp->V() > V);
  err = mp->I() - Inew;
  assert(err < tol && err > -tol);

  for (auto &su : mp->getSUs()) {
    err = su->I() - Inew / mp->getNSUs(); //!< we know the current has to split equally between both cells
    assert(err < tol && err > -tol);
    assert(su->V() > V);
    auto cell1 = dynamic_cast<Cell_Bucket *>(su.get()); //!< Dynamic cast from StorageUnit to Cell
    err = cell1->SOC() - (0.5);
    assert(err < tol && err > -tol);
  }

  return true;
}


template <typename Cell_t>
bool test_Modules_p()
{
  //!< test parallel modules with ECM cells
  double tol = settings::MODULE_P_I_ABSTOL;
  Cell_t *cell1;

  //!< setCurrent
  Deep_ptr<StorageUnit> cs[] = {
    Deep_ptr<StorageUnit>(new Cell_t()),
    Deep_ptr<StorageUnit>(new Cell_t())
  };

  double v1 = cs[0]->V();

  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  assert(mp->I() == 0);
  assert(mp->V() == v1);
  //!< discharge
  double Inew = 1.0 * std::size(cs);
  double V;
  mp->setCurrent(Inew, true);
  V = mp->V();
  assert(NEAR(mp->I(), Inew, tol));
  assert(V < v1); //!< voltage must decrease
  //!< do not check individual cells, that is done in getCells
  //!< charge
  Inew = -1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  V = mp->V();
  assert(NEAR(mp->I(), Inew, tol));
  assert(mp->V() > v1); //!< voltage must increase
  //!< rest with different SOC values
  //!< -> must be different lithium fractions since SOC does not affect the concentration
  //!< too complicated so skip this test here

  //!< validCells
  Deep_ptr<StorageUnit> cs2[] = {
    Deep_ptr<StorageUnit>(new Cell_t()),
    Deep_ptr<StorageUnit>(new Cell_t())
  };
  mp->setSUs(cs2);
  //!< valid cells with the current cells

  //!< CC timestep
  Deep_ptr<StorageUnit> cs3[] = {
    Deep_ptr<StorageUnit>(new Cell_t()),
    Deep_ptr<StorageUnit>(new Cell_t())
  };

  auto cp1 = dynamic_cast<Cell_t *>(cs3[0].get());
  auto cp2 = dynamic_cast<Cell_t *>(cs3[1].get());

  double soc1 = cp1->SOC();
  v1 = cp1->V();
  mp->setSUs(cs3, checkCells, true);
  //!< time step with 0 current
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(mp->V() == cp1->V());
  //!< discharge
  Inew = 1 * std::size(cs3);
  double err;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  V = mp->V();
  assert(V < v1);
  err = mp->I() - Inew;
  assert(err < tol && err > -tol); //!< note: this is not going to be exact because we solve a nonlinear equation to obain I
  //!< check individual cells

  for (auto &su : mp->getSUs()) {
    err = su->I() - Inew / mp->getNSUs(); //!< we know the current has to split equally between both cells
    assert(err < tol && err > -tol);
    assert(su->V() < v1);
    cell1 = dynamic_cast<Cell_t *>(su.get()); //!< Dynamic cast from StorageUnit to Cell
    err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
    assert(cell1->SOC() < soc1);
    assert(err < tol && err > -tol);
  }
  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  auto &cs5 = mp->getSUs();
  assert(mp->V() > V);
  err = mp->I() - Inew;
  assert(err < tol && err > -tol);
  for (int i = 0; i < std::size(cs5); i++) {
    err = cs5[i]->I() - Inew / std::size(cs5); //!< we know the current has to split equally between both cells
    assert(err < tol && err > -tol);
    assert(cs5[i]->V() > V);
    cell1 = dynamic_cast<Cell_t *>(cs5[i].get()); //!< Dynamic cast from StorageUnit to Cell
    err = cell1->SOC() - (0.5);
    assert(err < tol && err > -tol);
  }

  return true;
}


bool test_contactR()
{
  /*
   * Make a module with 3 cells and a contact resistance
   */

  double Rc = 0.01;
  double tol = 0.0001;
  Deep_ptr<StorageUnit> cs[] = {
    make<Cell_Bucket>(),
    make<Cell_Bucket>(),
    make<Cell_Bucket>()
  };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());
  auto cp3 = dynamic_cast<Cell_Bucket *>(cs[2].get());

  double Rcs[] = { Rc, Rc, Rc };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  mp->setRcontact(Rcs);

  //!< total resistance:
  //!< 		-Rp+-Rp-+-Rp-|
  //!< 		   |    |    |
  //!< 		  Rs    Rs   Rs
  //!< where Rp = contact resistance (value Rc) and Rs = cell resistance = 0.01;
  double Rs = cp1->getRtot();
  double Rq = Rc + (Rs * (Rs + Rc)) / (Rc + 2 * Rs); //!< resistance of last two branches
  double Rtot = Rc + Rs * Rq / (Rs + Rq);
  assert(NEAR(mp->getRtot(), Rtot, tol));

  //!< setCurrent
  //!< check voltages at each node from the branch going 'down' and the branch going 'right'
  //!< 	R1*I1 = Rp2*(I2 + I3) + R2*I2
  //!< 	R2*I2 = (Rp3 + R3)*I3
  //!< 	where Ri = resistance of cell i
  //!< 		  Rpi = contact resistance in parallel at cell i
  //!< 	since all cells have the same OCV
  double I = 20;
  mp->setCurrent(I, true, true);
  double I1 = cp1->I();
  double I2 = cp2->I();
  double I3 = cp3->I();

  //!< assert the currents of cells further from the terminal are smaller
  assert(std::abs(I1) > std::abs(I2));
  assert(std::abs(I2) > std::abs(I3));

  //!< Check the voltage at every node
  double Rcell = cp1->getRtot();                //!< all cell resistances are the same
  double V11 = Rcell * I1;                      //!< voltage at the node connecting the first cell, going down [ignoring OCV]
  double V12 = Rcs[1] * (I2 + I3) + Rcell * I2; //!< voltage at the node connecting the first cell, going right
  //!< double V13 = Rcs[1] * (I2 + I3) + Rcs[2]*I3 + Rcell*I3;
  double V22 = Rcell * I2;            //!< voltage at node of 2nd cell going down
  double V23 = (Rcs[2] + Rcell) * I3; //!< voltage at node of 2nc cell going right
  assert(NEAR(V11, V12, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V22, V23, settings::MODULE_P_V_ABSTOL));

  //!< check the total voltage
  double V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
  double V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
  double V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
  assert(NEAR(V1, V2, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V1, V3, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V2, V3, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(mp->V(), V2, settings::MODULE_P_V_ABSTOL));

  // #TODO needs to have getVall but it is protected.
  // assert(NEAR(V1, mp->Vi(0), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V2, mp->Vi(1), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V3, mp->Vi(2), tol)); //!< these numbers should be exactly the same

  //!< set charging current
  I = -20;
  mp->setCurrent(I, true, true);
  I1 = cp1->I();
  I2 = cp2->I();
  I3 = cp3->I();
  assert(std::abs(I1) > std::abs(I2));
  assert(std::abs(I2) > std::abs(I3));
  Rcell = cp1->getRtot();                //!< all cell resistances are the same
  V11 = Rcell * I1;                      //!< voltage at the node connecting the first cell, going down [ignoring OCV]
  V12 = Rcs[1] * (I2 + I3) + Rcell * I2; //!< voltage at the node connecting the first cell, going right
  V22 = Rcell * I2;                      //!< voltage at node of 2nd cell going down
  V23 = (Rcs[2] + Rcell) * I3;           //!< voltage at node of 2nc cell going right
  assert(NEAR(V11, V12, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V22, V23, settings::MODULE_P_V_ABSTOL));

  V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
  V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
  V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
  assert(NEAR(V1, V2, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V1, V3, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(V2, V3, settings::MODULE_P_V_ABSTOL));
  assert(NEAR(mp->V(), V2, settings::MODULE_P_V_ABSTOL));

  // #TODO needs to have getVall but it is protected.
  // assert(NEAR(V1, mp->Vi(0), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V2, mp->Vi(1), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V3, mp->Vi(2), tol)); //!< these numbers should be exactly the same
  return true;
}

bool test_Hierarchichal_p()
{
  //!< test parallel modules made out of other parallel modules
  double tol = settings::MODULE_P_I_ABSTOL;
  std::string ids[] = { "H1", "H2", "H3" };
  Deep_ptr<StorageUnit> SU1[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  Deep_ptr<StorageUnit> SU2[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  Deep_ptr<StorageUnit> SU3[] = { make<Cell_Bucket>(), make<Cell_Bucket>(), make<Cell_Bucket>() };
  auto cp1 = dynamic_cast<Cell_Bucket *>(SU1[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(SU1[1].get());
  auto cp3 = dynamic_cast<Cell_Bucket *>(SU2[0].get());
  auto cp4 = dynamic_cast<Cell_Bucket *>(SU2[1].get());
  auto cp5 = dynamic_cast<Cell_Bucket *>(SU3[0].get());
  auto cp6 = dynamic_cast<Cell_Bucket *>(SU3[1].get());
  auto cp7 = dynamic_cast<Cell_Bucket *>(SU3[2].get());

  double T = settings::T_ENV;
  bool checkCells = false;
  Deep_ptr<StorageUnit> MU[] = {
    make<Module_p>(ids[0], T, true, false, std::size(SU1), 1, 2),
    make<Module_p>(ids[1], T, true, false, std::size(SU2), 1, 2),
    make<Module_p>(ids[2], T, true, false, std::size(SU3), 1, 2)
  };

  auto mp1 = dynamic_cast<Module_p *>(MU[0].get()); //!< pass through cool systems
  auto mp2 = dynamic_cast<Module_p *>(MU[1].get());
  auto mp3 = dynamic_cast<Module_p *>(MU[2].get());

  mp1->setSUs(SU1, checkCells);
  mp2->setSUs(SU2, checkCells);
  mp3->setSUs(SU3, checkCells);

  //!< make the hierarichical module
  std::string n4 = "4";
  checkCells = true;
  auto mp = make<Module_p>(n4, T, true, true, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  assert(NEAR(Vini, mp1->V(), tol));
  assert(NEAR(Vini, cp5->V(), tol));
  assert(mp->getFullID() == "4");
  assert(mp1->getFullID() == "4_H1");
  assert(cp1->getFullID() == "4_H1_Cell_ECM<0>");
  assert(cp4->getFullID() == "4_H2_Cell_ECM<0>");
  assert(cp5->getFullID() == "4_H3_Cell_ECM<0>");

  //!< set a CC current
  double Inew = -14;    //!< should give about 2A per cell
  mp->setCurrent(Inew); // #TODO cannot set current !
  // assert(NEAR(mp->I(), Inew, tol));
  assert(NEAR(mp1->V(), mp2->V(), settings::MODULE_P_V_ABSTOL));
  assert(NEAR(mp3->V(), mp2->V(), settings::MODULE_P_V_ABSTOL));

  //!< time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(std::abs(cp1->SOC() - (0.5 - 2 * dt / 3600.0 / cp1->Cap())) < tol); //!< the SOC must have increased (check just 1 cell out of all 7)
  assert(mp->V() > Vini);
  assert(NEAR(mp2->V(), mp3->V(), tol)); //!< submodules must have same voltage

  return true;
}

bool test_Hierarchical_cross_p()
{
  //!< test parallel module made out of series modules
  //!< note: series modules must have same number of cells to get the same voltage
  double tol = settings::MODULE_P_I_ABSTOL;
  std::string ids[] = { "H1", "H2", "H3" };
  Deep_ptr<StorageUnit> SU1[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  Deep_ptr<StorageUnit> SU2[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  Deep_ptr<StorageUnit> SU3[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };

  double cap1 = SU1[0]->Cap();
  double v5 = SU3[0]->V();
  double T = settings::T_ENV;
  bool checkCells = false;

  Deep_ptr<StorageUnit> MU[] = {
    make<Module_s>(ids[0], T, true, false, std::size(SU1), 1, 2),
    make<Module_s>(ids[1], T, true, false, std::size(SU2), 1, 2),
    make<Module_s>(ids[2], T, true, false, std::size(SU3), 1, 2)
  };


  auto mp1 = dynamic_cast<Module_s *>(MU[0].get()); //!< pass through cool systems
  auto mp2 = dynamic_cast<Module_s *>(MU[1].get());
  auto mp3 = dynamic_cast<Module_s *>(MU[2].get());

  mp1->setSUs(SU1, checkCells);
  mp2->setSUs(SU2, checkCells);
  mp3->setSUs(SU3, checkCells);

  //!< make the hierarichical module
  std::string n4 = "4";
  checkCells = true;
  auto mp = make<Module_p>(n4, T, true, true, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  assert(NEAR(Vini, mp1->V(), tol));
  assert(std::abs(Vini - v5 * 2) < tol); //!< one module has 2 cells so voltage should split in 2

  //!< set a CC current
  double Inew = -6; //!< should give about 2A per cell
  mp->setCurrent(Inew);
  assert(NEAR(mp->I(), Inew, tol));

  // assert(std::abs(cp3->I() + 2) < tol);
  // assert(std::abs(cp6->I() + 2) < tol);
  assert(std::abs(mp1->I() + 2) < tol);  //!< m1 has two cells
  assert(std::abs(mp3->I() + 2) < tol);  //!< m3 has two cells
  assert(NEAR(mp1->V(), mp3->V(), tol)); //!< check voltage is equal

  //!< time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  //  assert(std::abs(cp1->SOC() - (0.5 - 2 * dt / 3600.0 / cap1)) < tol); //!< the SOC must have increased (check just 1 cell out of all 7)
  assert(mp->V() > Vini);
  assert(NEAR(mp2->V(), mp3->V(), tol));
  //!< submodules must have same voltage
  //!< note: there is no check on sub-modules with different SOC but I assume that works since it works with sub-cells of different SOC

  return true;
}

bool test_copy_p()
{
  //!< 	/*
  //!< 	 * test the copy-function
  //!< 	 */

  //!< make module #TODO copy functions are commented out.
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());

  std::string n = "na";
  double v1 = cp1->V();
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< copy this one and check they are identical
  auto cn = mp;
  // Deep_ptr<StorageUnit> cn(mp->copy());
  Module_p *c22 = dynamic_cast<Module_p *>(cn.get()); //!< Dynamic cast from StorageUnit to Cell
  assert(mp->V() == c22->V());

  auto &corig = mp->getSUs();
  auto &cnew = c22->getSUs();
  for (int i = 0; i < mp->getNSUs(); i++)
    assert(corig[i]->V() == cnew[i]->V());

  //!< change the copied version, and ensure the old one is still the same
  c22->setCurrent(1 * std::size(cs), false, false); //!< discharge
  for (int t = 0; t < 10; t++)
    c22->timeStep_CC(2);

  for (int i = 0; i < mp->getNSUs(); i++) {
    assert(corig[i]->V() == v1);
    assert(cnew[i]->V() < v1);
  }

  return true;
}

bool test_equaliseV_timing(Deep_ptr<Module_p> &mp, Deep_ptr<StorageUnit> c[], int nin)
{
  //!< test timing
  //!< IN
  //!< mp		parallel module
  //!< SUs 		array with smart pointers to the children of mp
  mp->setBlockDegAndTherm(true); //!< ignore thermal and degradation during this function (we mess with individual cells time time keeping for thermal gives errors)

  //!< set a 1C current to the individual cells, then redistribute
  double I = mp->Cap();
  double Ii = I / mp->getNcells(); //!< current per cell
  for (int i = 0; i < nin; i++)
    c[i]->setCurrent(Ii * c[i]->getNcells());

  std::cout << "after setCurrent, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  mp->redistributeCurrent(true);
  std::cout << "while after n1 steps in redistributeCurrent, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  //!< take 10 time steps, then redistribute again
  int N = 10;
  double dt = 2;
  for (auto &su : mp->getSUs()) su->timeStep_CC(dt, N);

  std::cout << "after 10 time steps, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  mp->redistributeCurrent(true);

  std::cout << "while after n2 steps in redistributeCurrent, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  //!< go to low SOC (lowest cell voltage is 3V
  while (mp->getVlow() > 3.0)
    mp->timeStep_CC(dt, 1);

  //!< take 10 time steps, then redistribute again
  N = 10;
  for (auto &su : mp->getSUs()) su->timeStep_CC(dt, N);

  std::cout << "after discharge, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  mp->redistributeCurrent(true);
  std::cout << "while after n3 steps in redistributeCurrent, the cell voltages are ";
  for (auto &su : mp->getSUs())
    std::cout << su->V() << " ";
  std::cout << '\n';

  std::cout << "Number of time steps for mp " << mp->getFullID()
            << " is for setCurrent n1, for timestep n2"
            << ", for timestep at low SOC n3, now starting a CC cycle\n";

  //!< Use a Cycler to do a full CC cycle with redistributeCurrent every time step
  Cycler cyc;
  ThroughputData th;
  double lim = 0.0;
  int ndata = 0;
  double vlim;
  cyc.initialise(mp.get(), mp->getFullID());
  vlim = mp->Vmax() - lim;
  cyc.CC(-I, vlim, TIME_INF, dt, ndata, th); //!< CC charge
  vlim = mp->Vmin() + lim;
  cyc.CC(I, vlim, TIME_INF, dt, ndata, th); //!< CC discharge

  std::cout << "Finished CC cycle.\n";

  return true;
}

bool test_equaliseV()
{
  //!< test timing with
  //!< 		5 identical cells
  //!< 		5 cells with minor differences
  //!< 		5 widely different cells
  //!< 		4 similar and one aged cell

  DEG_ID deg;
  deg.SEI_id.add_model(4); //!< chirstensen SEI growth
  deg.SEI_porosity = 0;    //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0); //!< no surface cracks
  deg.CS_diffusion = 0;   //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  deg.LAM_id.add_model(0); //!< no LAM
  deg.pl_id = 0;           //!< no litihium plating
  double T2 = settings::T_ENV;
  bool checkCells2 = false;

  //!< 5 identical cells
  int ncel1 = 5;
  std::string n1 = "mp_identical";
  Deep_ptr<StorageUnit> cs1[] = {
    make<Cell_SPM>("cell1", deg, 1, 1, 1, 1),
    make<Cell_SPM>("cell2", deg, 1, 1, 1, 1),
    make<Cell_SPM>("cell3", deg, 1, 1, 1, 1),
    make<Cell_SPM>("cell4", deg, 1, 1, 1, 1),
    make<Cell_SPM>("cell5", deg, 1, 1, 1, 1)
  };

  auto mpp1 = make<Module_p>(n1, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp1->setSUs(cs1, checkCells2, true);
  test_equaliseV_timing(mpp1, cs1, ncel1);

  //!< 5 cells with small distribution
  std::default_random_engine gen;
  double std1 = 0;
  double std2 = 0;
  std1 = 0.004;
  std2 = 0.025;
  std::normal_distribution<double> distr_c(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate
  std::string n2 = "mp_variation";

  Deep_ptr<StorageUnit> cs2[] = {
    make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell9", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen))
  };

  auto mpp2 = make<Module_p>(n2, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp2->setSUs(cs2, checkCells2, true);
  test_equaliseV_timing(mpp2, cs2, ncel1);

  //!< 5 cells with large distribution
  std1 = 0.1;
  std2 = 0.15;
  std::normal_distribution<double> distr_c2(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r2(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d2(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate
  std::string n3 = "mp_largeVariation";

  Deep_ptr<StorageUnit> cs3[] = {
    make<Cell_SPM>("cell5", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)),
    make<Cell_SPM>("cell6", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)),
    make<Cell_SPM>("cell7", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)),
    make<Cell_SPM>("cell8", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)),
    make<Cell_SPM>("cell9", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen))
  };

  auto mpp3 = make<Module_p>(n3, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp3->setSUs(cs3, checkCells2, true);
  test_equaliseV_timing(mpp3, cs3, ncel1);

  //!< 4 similar and one very different
  std::string n4 = "mp_4and1";
  Deep_ptr<StorageUnit> cs4[] = {
    make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell9", deg, 0.5, 2.0, 1.1, 1.1) //!< one with half the capacity and double the resistance
  };

  auto mpp4 = make<Module_p>(n4, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp4->setSUs(cs4, checkCells2, true);
  test_equaliseV_timing(mpp4, cs4, ncel1);

  return true;
}

int test_all_Module_p()
{
  /*
   * Test the functions from the parallel module
   * note that we already test the function from the base module in the unit test for the series-connected module so there is no point to repeat them
   */

  if (!TEST(test_Constructor_p, "test_Constructor_p")) return 1;
  if (!TEST(test_BasicGetters_p, "test_BasicGetters_p")) return 2;
  if (!TEST(test_setI_p, "test_setI_p")) return 3;
  if (!TEST(test_validStates_p, "test_validStates_p")) return 5;
  if (!TEST(test_timeStep_CC_p, "test_timeStep_CC_p")) return 6;
  if (!TEST(test_contactR, "test_contactR")) return 7;

  // //!< Combinations
  if (!TEST(test_Modules_p<Cell_ECM<1>>, "test_Modules_p<Cell_ECM<1>>")) return 8; //!< parallel from ECM cells
  if (!TEST(test_Modules_p<Cell_SPM>, "test_Modules_p<Cell_SPM>")) return 9;       //!< parallel from SPM cells
  if (!TEST(test_Hierarchichal_p, "test_Hierarchichal_p")) return 10;              //!< parallel from parallel
  if (!TEST(test_Hierarchical_cross_p, "test_Hierarchical_cross_p")) return 11;    //!< parallel from series
  if (!TEST(test_copy_p, "test_copy_p")) return 12;

  return 0;
}

} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Module_p(); }