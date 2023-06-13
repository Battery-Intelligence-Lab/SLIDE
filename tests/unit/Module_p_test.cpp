/*
 * Module_base_p_test.cpp
 *
 *  Created on: 18 Dec 2019
 *  Updated on: 13 Jun 2023 (Catch2)
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../../src/slide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <cassert>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::Equals;
using namespace slide;

constexpr double TOL_EQ = 1e-15;

TEST_CASE("test_Constructor_p")
{
  auto mp = make<Module_p>();
  REQUIRE(mp->getNSUs() == 0);
  //  assert(mp->T() == settings::T_ENV); #TODO it returns cool->T() which is nullptr.

  auto cp1 = make<Cell_Bucket>();
  auto cp2 = make<Cell_Bucket>();
  REQUIRE_THAT(cp1->getID(), Equals("Cell_ECM_0_"));
  REQUIRE_THAT(cp1->getFullID(), Equals("Cell_ECM_0_")); // Has no parent yet

  Deep_ptr<StorageUnit> cs[] = { std::move(cp1), std::move(cp2) };
  std::string n = "na";
  const double T = settings::T_ENV;
  bool checkCells = false;
  auto mp2 = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp2->setSUs(cs, checkCells, true);

  REQUIRE(mp2->getNSUs() == std::size(cs));
  REQUIRE_THAT(mp2->T(), WithinAbs(T, TOL_EQ));
}

TEST_CASE("test_BasicGetters_p")
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());
  REQUIRE_THAT(cp1->getID(), Equals("Cell_ECM_0_", Catch::CaseSensitive::Yes));
  REQUIRE_THAT(cp1->getFullID(), Equals("Cell_ECM_0_", Catch::CaseSensitive::Yes)); // Has no parent yet

  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  REQUIRE(mp->Cap() == std::size(cs) * cp1->Cap());
  REQUIRE(mp->Vmin() == cp1->Vmin());
  REQUIRE(mp->Vmax() == cp1->Vmax());
  REQUIRE(mp->VMIN() == cp1->VMIN());
  REQUIRE(mp->VMAX() == cp1->VMAX());
  REQUIRE(mp->I() == 0);
  REQUIRE(mp->V() == cp1->V());
}

TEST_CASE("test_setI_p")
{
  constexpr double tol = 0.005;
  double Inew;

  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());

  // Check for correct ID and full ID
  REQUIRE_THAT(cp1->getID(), Equals("Cell_ECM_0_"));
  REQUIRE_THAT(cp1->getFullID(), Equals("Cell_ECM_0_"));

  double v1 = cp1->V();
  std::string n = "na";
  constexpr double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  // Initial conditions
  REQUIRE_THAT(mp->I(), WithinAbs(0, TOL_EQ));
  REQUIRE_THAT(mp->V(), WithinAbs(cp1->V(), TOL_EQ));

  // Test for discharge condition
  Inew = 1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));
  REQUIRE(mp->V() < v1); //!< voltage must decrease
  //!< do not check individual cells, that is done in getCells #TODO check if done

  // Test for charge condition
  Inew = -1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));
  REQUIRE(mp->V() > v1); //!< voltage must increase

  // Rest with different SOC values
  Inew = 0;
  cp2->setSOC(0.4); // c2 has lower OCV -> should charge
  mp->setCurrent(Inew, true);
  REQUIRE_THAT(cp1->V(), WithinAbs(cp2->V(), settings::MODULE_P_V_ABSTOL));
  REQUIRE(cp1->I() > 0); //!< cell currents are opposite
  REQUIRE(cp2->I() < 0);

  // Test large currents that should fail
  Inew = 10000;                         // Large current causing low voltage
  Status status = mp->setCurrent(Inew); // Expected to fail due to unsolvable current equation
  REQUIRE_FALSE(isStatusSuccessful(status));

  Inew = -10000;                 // Similarly large current causing low voltage
  status = mp->setCurrent(Inew); // Expected to fail as well
  REQUIRE_FALSE(isStatusSuccessful(status));
}

TEST_CASE("test_validStates_p")
{
  Deep_ptr<StorageUnit> cs[] = { make<Cell_Bucket>(), make<Cell_Bucket>() };
  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());
  auto cp2 = dynamic_cast<Cell_Bucket *>(cs[1].get());

  // Validate ID and Full ID
  REQUIRE_THAT(cp1->getID(), Equals("Cell_ECM_0_"));
  REQUIRE_THAT(cp1->getFullID(), Equals("Cell_ECM_0_"));

  std::string n = "na";
  constexpr double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  // Test valid states
  std::vector<double> s;
  mp->getStates(s);
  REQUIRE(mp->validStates());
}

TEST_CASE("test_timeStep_CC_p")
{
  constexpr double T = settings::T_ENV;
  constexpr bool checkCells = false;
  constexpr double tol = 1e-4;

  Deep_ptr<StorageUnit> cs[] = {
    Deep_ptr<StorageUnit>(new Cell_Bucket()),
    Deep_ptr<StorageUnit>(new Cell_Bucket())
  };

  auto cp1 = dynamic_cast<Cell_Bucket *>(cs[0].get());

  std::string n = "na";
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  const double v1{ cp1->V() }, soc1{ cp1->SOC() };

  // time step with 0 current
  constexpr double dt = 5;
  mp->timeStep_CC(dt);
  REQUIRE_THAT(mp->V(), WithinAbs(cp1->V(), TOL_EQ));

  // discharge
  double Inew = 1 * std::size(cs);
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  double V = mp->V();
  REQUIRE(V < v1);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));

  for (auto &su : mp->getSUs()) {
    REQUIRE_THAT(su->I(), WithinAbs(Inew / mp->getNSUs(), tol));
    REQUIRE(su->V() < v1);
    auto cell1 = dynamic_cast<Cell_Bucket *>(su.get()); // Dynamic cast from StorageUnit to Cell
    REQUIRE(cell1->SOC() < soc1);
    REQUIRE_THAT(cell1->SOC(), WithinAbs((0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap()), tol));
  }

  // charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  REQUIRE(mp->V() > V);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));

  for (auto &su : mp->getSUs()) {
    REQUIRE_THAT(su->I(), WithinAbs(Inew / mp->getNSUs(), tol));
    REQUIRE(su->V() > V);                               // #TODO why V not v1?
    auto cell1 = dynamic_cast<Cell_Bucket *>(su.get()); // Dynamic cast from StorageUnit to Cell
    REQUIRE_THAT(cell1->SOC(), WithinAbs(0.5, tol));
  }
}

TEMPLATE_TEST_CASE("test_Modules_p", "[Module_p]", Cell_ECM<1>, Cell_SPM)
{
  constexpr double tol = 1e-4; // We set a higher tolerance for SOC but current tolerance can be less.

  // Set current
  Deep_ptr<StorageUnit> cs[] = { make<TestType>(), make<TestType>() };
  const double v1 = cs[0]->V();
  std::string n = "na";
  constexpr double T = settings::T_ENV;
  constexpr bool checkCells = false;

  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  REQUIRE_THAT(mp->I(), WithinAbs(0.0, TOL_EQ));
  REQUIRE_THAT(mp->V(), WithinAbs(v1, TOL_EQ));

  // Discharge
  double Inew = 1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  double V = mp->V();
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));
  REQUIRE(V < v1); // Voltage must decrease

  // Charge
  Inew = -1.0 * std::size(cs);
  mp->setCurrent(Inew, true);
  V = mp->V();
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));
  REQUIRE(mp->V() > v1); // Voltage must increase

  // CC timestep
  Deep_ptr<StorageUnit> cs3[] = { make<TestType>(), make<TestType>() };

  auto cp1 = dynamic_cast<TestType *>(cs3[0].get());

  const double soc1 = cp1->SOC();
  const double v1_cc = cp1->V();
  mp->setSUs(cs3, checkCells, true);

  // Time step with 0 current
  const double dt = 5;
  mp->timeStep_CC(dt);
  REQUIRE_THAT(mp->V(), WithinAbs(cp1->V(), TOL_EQ));

  // Discharge
  Inew = 1 * std::size(cs3);
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  V = mp->V();
  REQUIRE(V < v1_cc);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));


  for (auto &su : mp->getSUs()) {
    double expected_I = Inew / mp->getNSUs(); // We know the current has to split equally between both cells
    REQUIRE_THAT(su->I(), WithinAbs(expected_I, TOL_EQ));
    REQUIRE(su->V() < v1_cc);
    auto cell1 = dynamic_cast<TestType *>(su.get()); // Dynamic cast from StorageUnit to Cell
    REQUIRE(cell1->SOC() < soc1);
    double expected_SOC = 0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap();
    REQUIRE_THAT(cell1->SOC(), WithinAbs(expected_SOC, TOL_EQ));
  }

  // Charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  auto &cs5 = mp->getSUs();
  REQUIRE(mp->V() > V);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));

  for (size_t i = 0; i < std::size(cs5); i++) {
    double expected_I = Inew / std::size(cs5); // We know the current has to split equally between both cells
    REQUIRE_THAT(cs5[i]->I(), WithinAbs(expected_I, TOL_EQ));
    REQUIRE(cs5[i]->V() > V);
    auto cell1 = dynamic_cast<TestType *>(cs5[i].get()); // Dynamic cast from StorageUnit to Cell
    REQUIRE_THAT(cell1->SOC(), WithinAbs(0.5, tol));     // #TODO due to Cell_SPM it throws
  }
}

TEST_CASE("test_contactR", "[Module_p]")
{
  // Make a module with 3 cells and a contact resistance
  double Rc = 0.01;
  constexpr double tol = 0.0001;
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

  // Calculating total resistance
  double Rs = cp1->getRtot();
  double Rq = Rc + (Rs * (Rs + Rc)) / (Rc + 2 * Rs); //!< resistance of last two branches
  double Rtot = Rc + Rs * Rq / (Rs + Rq);
  REQUIRE_THAT(mp->getRtot(), WithinAbs(Rtot, tol));

  //!< setCurrent
  //!< check voltages at each node from the branch going 'down' and the branch going 'right'
  //!< 	R1*I1 = Rp2*(I2 + I3) + R2*I2
  //!< 	R2*I2 = (Rp3 + R3)*I3
  //!< 	where Ri = resistance of cell i
  //!< 		  Rpi = contact resistance in parallel at cell i
  //!< 	since all cells have the same OCV
  double I = 20;
  mp->setCurrent(I, true, true);
  double I1{ cp1->I() }, I2{ cp2->I() }, I3{ cp3->I() };

  // Asserting the currents of cells further from the terminal are smaller
  REQUIRE(std::abs(I1) > std::abs(I2));
  REQUIRE(std::abs(I2) > std::abs(I3));

  // Checking the voltage at every node
  double Rcell = cp1->getRtot();                //!< all cell resistances are the same
  double V11 = Rcell * I1;                      //!< voltage at the node connecting the first cell, going down [ignoring OCV]
  double V12 = Rcs[1] * (I2 + I3) + Rcell * I2; //!< voltage at the node connecting the first cell, going right
  //!< double V13 = Rcs[1] * (I2 + I3) + Rcs[2]*I3 + Rcell*I3;
  double V22 = Rcell * I2;
  double V23 = (Rcs[2] + Rcell) * I3;
  REQUIRE_THAT(V11, WithinAbs(V12, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V22, WithinAbs(V23, settings::MODULE_P_V_ABSTOL));

  // Checking the total voltage
  double V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
  double V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
  double V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
  REQUIRE_THAT(V1, WithinAbs(V2, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V1, WithinAbs(V3, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V2, WithinAbs(V3, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(mp->V(), WithinAbs(V2, settings::MODULE_P_V_ABSTOL));

  // #TODO needs to have getVall but it is protected.
  // assert(NEAR(V1, mp->Vi(0), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V2, mp->Vi(1), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V3, mp->Vi(2), tol)); //!< these numbers should be exactly the same

  // Setting charging current
  I = -20;
  mp->setCurrent(I, true, true);
  I1 = cp1->I(), I2 = cp2->I(), I3 = cp3->I();

  // Asserting the currents of cells further from the terminal are smaller
  REQUIRE(std::abs(I1) > std::abs(I2));
  REQUIRE(std::abs(I2) > std::abs(I3));

  // Checking the voltage at every node
  Rcell = cp1->getRtot();
  V11 = Rcell * I1;
  V12 = Rcs[1] * (I2 + I3) + Rcell * I2;
  V22 = Rcell * I2;
  V23 = (Rcs[2] + Rcell) * I3;
  REQUIRE_THAT(V11, WithinAbs(V12, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V22, WithinAbs(V23, settings::MODULE_P_V_ABSTOL));

  // Checking the total voltage
  V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
  V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
  V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
  REQUIRE_THAT(V1, WithinAbs(V2, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V1, WithinAbs(V3, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(V2, WithinAbs(V3, settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(mp->V(), WithinAbs(V2, settings::MODULE_P_V_ABSTOL));

  // #TODO needs to have getVall but it is protected.
  // assert(NEAR(V1, mp->Vi(0), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V2, mp->Vi(1), tol)); //!< these numbers should be exactly the same
  // assert(NEAR(V3, mp->Vi(2), tol)); //!< these numbers should be exactly the same
}

TEST_CASE("test_Hierarchichal_p", "[module]")
{
  //!< test parallel modules made out of other parallel modules
  double tol = 1e-3;
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

  //!< make the hierarchical module
  std::string n4 = "4";
  checkCells = true;
  auto mp = make<Module_p>(n4, T, true, true, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  REQUIRE_THAT(Vini, WithinAbs(mp1->V(), tol));
  REQUIRE_THAT(Vini, WithinAbs(cp5->V(), tol));
  REQUIRE_THAT(mp->getFullID(), Equals("4"));
  REQUIRE_THAT(mp1->getFullID(), Equals("4_H1"));
  REQUIRE_THAT(cp1->getFullID(), Equals("4_H1_Cell_ECM_0_"));
  REQUIRE_THAT(cp4->getFullID(), Equals("4_H2_Cell_ECM_0_"));
  REQUIRE_THAT(cp5->getFullID(), Equals("4_H3_Cell_ECM_0_"));

  //!< set a CC current
  double Inew = -14;    //!< should give about 2A per cell
  mp->setCurrent(Inew); // #TODO cannot set current !
                        // REQUIRE_THAT(mp->I(), Catch::WithinAbs(Inew, tol));

  REQUIRE_THAT(mp1->V(), WithinAbs(mp2->V(), settings::MODULE_P_V_ABSTOL));
  REQUIRE_THAT(mp3->V(), WithinAbs(mp2->V(), settings::MODULE_P_V_ABSTOL));

  //!< time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  REQUIRE_THAT(cp1->SOC(), WithinAbs((0.5 - 2 * dt / 3600.0 / cp1->Cap()), tol)); //!< the SOC must have increased (check just 1 cell out of all 7)
  REQUIRE(mp->V() > Vini);
  REQUIRE_THAT(mp2->V(), mp3->V(), tol);                                          //!< submodules must have same voltage
}

TEST_CASE("Hierarchical_cross_p", "[Module_p]")
{
  // Test parallel module made out of series modules
  // Note: series modules must have same number of cells to get the same voltage

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

  auto mp1 = dynamic_cast<Module_s *>(MU[0].get()); // Pass through cool systems
  auto mp2 = dynamic_cast<Module_s *>(MU[1].get());
  auto mp3 = dynamic_cast<Module_s *>(MU[2].get());

  mp1->setSUs(SU1, checkCells);
  mp2->setSUs(SU2, checkCells);
  mp3->setSUs(SU3, checkCells);

  // Make the hierarichical module
  std::string n4 = "4";
  checkCells = true;
  auto mp = make<Module_p>(n4, T, true, true, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  REQUIRE_THAT(Vini, WithinAbs(mp1->V(), tol));
  REQUIRE_THAT(Vini, WithinAbs(v5 * 2, tol)); // One module has 2 cells so voltage should split in 2

  // Set a CC current
  double Inew = -6; // Should give about 2A per cell
  mp->setCurrent(Inew);
  REQUIRE_THAT(mp->I(), WithinAbs(Inew, tol));
  REQUIRE_THAT(mp1->I(), WithinAbs(-2, tol));       // m1 has two cells
  REQUIRE_THAT(mp3->I(), WithinAbs(-2, tol));       // m3 has two cells
  REQUIRE_THAT(mp1->V(), WithinAbs(mp3->V(), tol)); // Check voltage is equal

  // Time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  REQUIRE(mp->V() > Vini);
  REQUIRE_THAT(mp2->V(), WithinAbs(mp3->V(), tol)); // Submodules must have same voltage
                                                    // Note: there is no check on sub-modules with different SOC but I assume that works since it works with sub-cells of different SOC
}

TEMPLATE_TEST_CASE("test_copy_p", "[Module_p]", Cell_Bucket, Cell_ECM<1>, Cell_SPM)
{
  // make module
  Deep_ptr<StorageUnit> cs[] = { make<TestType>(), make<TestType>() };

  auto cp1 = dynamic_cast<TestType *>(cs[0].get());

  std::string n = "na";
  const double v1 = cp1->V();
  constexpr double T = settings::T_ENV;
  constexpr bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  // copy this one and check they are identical
  auto cn = mp;
  auto c22 = dynamic_cast<Module_p *>(cn.get()); // Dynamic cast from StorageUnit to Module
  REQUIRE_THAT(mp->V(), WithinAbs(c22->V(), TOL_EQ));

  auto &corig = mp->getSUs();
  auto &cnew = c22->getSUs();
  for (size_t i = 0; i < mp->getNSUs(); ++i)
    REQUIRE_THAT(corig[i]->V(), WithinAbs(cnew[i]->V(), TOL_EQ));

  // change the copied version, and ensure the old one is still the same
  c22->setCurrent(1 * std::size(cs), false, false); // discharge
  for (int t = 0; t < 10; ++t)
    c22->timeStep_CC(2);

  for (size_t i = 0; i < mp->getNSUs(); ++i) {
    REQUIRE_THAT(corig[i]->V(), WithinAbs(v1, TOL_EQ));
    REQUIRE(cnew[i]->V() < v1);
  }
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
  ThroughputData th;
  double lim = 0.0;
  int ndata = 0;
  double vlim;
  Cycler cyc(mp.get(), mp->getFullID());
  vlim = mp->Vmax() - lim;
  cyc.CC(-I, vlim, TIME_INF, dt, ndata, th); //!< CC charge
  vlim = mp->Vmin() + lim;
  cyc.CC(I, vlim, TIME_INF, dt, ndata, th);  //!< CC discharge

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

  deg.CS_id.add_model(0);  //!< no surface cracks
  deg.CS_diffusion = 0;    //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

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

// int test_all_Module_p()
// {
//   /*
//    * Test the functions from the parallel module
//    * note that we already test the function from the base module in the unit test for the series-connected module so there is no point to repeat them
//    */

//   if (!TEST(test_contactR, "test_contactR")) return 7;

//   // //!< Combinations
//   if (!TEST(test_Modules_p<Cell_ECM<1>>, "test_Modules_p<Cell_ECM<1>>")) return 8; //!< parallel from ECM cells
//   if (!TEST(test_Modules_p<Cell_SPM>, "test_Modules_p<Cell_SPM>")) return 9;       //!< parallel from SPM cells
//   if (!TEST(test_Hierarchichal_p, "test_Hierarchichal_p")) return 10;              //!< parallel from parallel
//   if (!TEST(test_Hierarchical_cross_p, "test_Hierarchical_cross_p")) return 11;    //!< parallel from series

//   return 0;
// }