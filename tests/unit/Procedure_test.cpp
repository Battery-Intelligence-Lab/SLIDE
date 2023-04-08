/*
 * CyclerProcedure_test.cpp
 *
 *  Created on: 3 Mar 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <fstream>

namespace slide::tests::unit {
bool test_Procedure_cycleAge(double Rc, bool spread, int cool)
{
  /*
   * Test the ageing procedure with a few different SUs
   */

  //!< Settings
  std::default_random_engine gen;
  double std1{ 1e-10 }, std2{ 1e-10 }, std3{ 1e-10 };
  if (spread) {
    std1 = 0.004;
    std2 = 0.025;
    std3 = 0.1;
  }

  std::normal_distribution<double> distr_c(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d(1.0, std3); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate

  //!< Make the procedure with standard settings
  double Vbal = 3.5;
  int ndata = 20;
  bool balance = true;
  bool unittest = true;

  std::cout << "Procedure_test start cells\n";

  //!< test with SPM cell
  DEG_ID deg;
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  auto cp1 = make<Cell_SPM>("proctest_cell", deg, 1, 1, 1, 1);
  auto p = Procedure(balance, Vbal, ndata, unittest);

  p.cycleAge(cp1.get(), true); //!< WITH CV -> 1C CCCV cycling
  //!< just check this runs without producing an error warning

  std::cout << "Procedure_test start series module\n";

  //!< test with series module
  std::string n = "proctest_s_module_" + std::to_string(Rc);
  Deep_ptr<StorageUnit> cs[] = {
    make<Cell_SPM>("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen))
  };

  double Rcs[] = { Rc, Rc, Rc, Rc, Rc };
  constexpr double T = 25_degC;
  constexpr bool checkCells = false;

  auto mp = make<Module_s>(n, T, true, false, std::size(cs), cool, 1); //!< no multithreading

  mp->setSUs(cs, checkCells, true);
  mp->setRcontact(Rcs);
  auto p2 = Procedure(balance, Vbal, ndata, unittest);
  p2.cycleAge(mp.get(), false); //!< this should write a file called s_module_capacities.csv
                                //!< 	check that all cells age more or less the same
                                //!< 	the capacity of the string should be the capacity of the smallest cell
  std::cout << "Procedure_test start battery of parallel module\n";

  //!< test with Battery from parallel module
  constexpr size_t ncel2 = 9;
  std::string n22 = "proctest_p_module_" + std::to_string(Rc) + "_batt";
  std::string n2 = "mp1";

  Deep_ptr<StorageUnit> cs2[ncel2];

  for (size_t i = 0; i < ncel2; i++) {
    std::string name = "cell" + std::to_string(i);
    cs2[i] = make<Cell_SPM>(name, deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  }

  double Rcs2[] = { Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc };
  constexpr double T2 = 273 + 25;
  constexpr bool checkCells2 = false;

  auto mpp = make<Module_p>(n2, T2, true, false, std::size(cs2), cool, 2); //!< no multithreading, nt_Vcheck time steps between checking SU voltage

  mpp->setSUs(cs2, checkCells2, true);
  mpp->setRcontact(Rcs2);

  auto bat = make<Battery>(n22);
  bat->setModule(std::move(mpp));

  auto p3 = Procedure(balance, Vbal, ndata, unittest);

  p3.cycleAge(bat.get(), false); //!< this should write a file called p_module_capacities.csv
  //!< check the document with the capacities
  //!< 	if contact resistance is not zero, then cells with higher numbers (right columns)
  //!< 		should have more remaining capacities than cells with low numbers (left columns)
  //!< 	if Rc == 0, all cells should have more or less the same capacity

  return true;
}

bool test_Procedure_CoolSystem()
{
  /*
   * test the cool system design with proper cycle ageing
   *
   * NOTE: these tests use CycleAge with all coolcontrol settings, so they can take a while to run
   *
   */

  std::cout << "Procedure_test test coolsystems.\n";

  //!< Parameters from Cell_SPM which are needed to calculate the heat balance
  constexpr double rho = 1626;
  constexpr double Cp = 750;
  constexpr double L = 1.6850e-4;               //!< thickness of one layer
  constexpr double width = 0.1;                 //!< width of the pouch
  constexpr double height = 0.2;                //!< height of the pouch
  constexpr int nlayers = 31;                   //!< number of layers in the pouch
  constexpr double Acell = width * height;      //!< geometric surface area of the pouch
  constexpr double elec_surf = Acell * nlayers; //!< total 'unrolled' surface area of the electrodes

  //!< settings of the cycle age
  double Ccha = 1;
  double Cdis = 1;
  int Ncycle = 100;
  int ncheck = 10;
  int nbal = 10;

  //!< General settings
  double T = settings::T_ENV;
  bool checkCells = true;
  double Vbal = 3.5;
  int ndata = 20;
  bool balance = true;
  bool unitest = true; //!< controls printing for cycle age, should be true unless for debugging
  auto p = Procedure(balance, Vbal, ndata, unitest);

  //!< Loop for each setting of the cool controller
  for (int coolControl = 1; coolControl < 6; coolControl++) {
    //!< ****************************************************************************************************************************************************
    //!< Make a simple module with one SPM cell
    constexpr size_t ncel = 1;
    Deep_ptr<StorageUnit> cs[] = { make<Cell_SPM>() };

    auto cp0 = dynamic_cast<Cell_SPM *>(cs[0].get());

    std::string n = "testCoolSystem_cell";
    auto mp = make<Module_s>(n, T, true, false, ncel, coolControl, 1);

    double Tini[1] = { cs[0]->T() };
    mp->setSUs(cs, checkCells, true);

    //!< Use procedure age to cycle it
    //!< cout<<"Procedure_test start coolsystem test with a single cell for cool control setting "<<coolControl<<endl;
    p.cycleAge(mp.get(), Ncycle, ncheck, nbal, true, Ccha, Cdis, mp->Vmax(), mp->Vmin()); //!< include CV phase

    //!< check the energy balance of the outer module
    double Qgen = cp0->thermal_getTotalHeat();         //!< total heat generated by cells
    double Qcool = mp->getCoolSystem()->getHeatEvac(); //!< total heat extracted by the coolsystem from the cells
    double Tnew[1] = { cp0->T() };
    double Qheat = 0; //!< total energy in heating up the cells
    for (int i = 0; i < 1; i++)
      Qheat += (Tnew[i] - Tini[i]) * (rho * Cp * L * elec_surf);
    assert(std::abs(Qgen - Qcool - Qheat) <= 1e-10 * std::abs(Qgen)); // Equality to include zero condition.

    //!< **********************************************************************************************************************************************************
    //!< Make a simple module with SPM cells
    constexpr size_t ncel2 = 4;
    Deep_ptr<StorageUnit> cs2[ncel2];
    Cell_SPM *cs2_ptr[ncel2];
    for (size_t i = 0; i < ncel2; i++) {
      auto temp = make<Cell_SPM>();
      cs2_ptr[i] = temp.get();
      cs2[i] = std::move(temp);
    }

    std::string n2 = "testCoolSystem_module";
    auto mp2 = make<Module_s>(n2, T, true, false, std::size(cs2), coolControl, 1);
    mp2->setSUs(cs2, checkCells, true);
    double Tini2[4] = { cs2_ptr[0]->T(), cs2_ptr[1]->T(), cs2_ptr[2]->T(), cs2_ptr[3]->T() };

    //!< Use procedure age to cycle it
    //!< cout<<"Procedure_test start coolsystem test with a simple module for cool control setting "<<coolControl<<endl;
    p.cycleAge(mp2.get(), Ncycle, ncheck, nbal, false, Ccha, Cdis, mp2->Vmax(), mp2->Vmin()); //!< exclude CV phase, CC only

    //!< check the energy balance of the outer module
    double Qgen2 = cs2_ptr[0]->thermal_getTotalHeat() + cs2_ptr[1]->thermal_getTotalHeat() + cs2_ptr[2]->thermal_getTotalHeat() + cs2_ptr[3]->thermal_getTotalHeat(); //!< total heat generated by cells
    double Qcool2 = mp2->getCoolSystem()->getHeatEvac();                                                                                                              //!< total heat extracted by the coolsystem from the cells
    double Tnew2[] = { cs2_ptr[0]->T(), cs2_ptr[1]->T(), cs2_ptr[2]->T(), cs2_ptr[3]->T() };
    double Qheat2 = 0; //!< total energy in heating up the cells
    for (int i = 0; i < std::size(Tnew2); i++)
      Qheat2 += (Tnew2[i] - Tini2[i]) * (rho * Cp * L * elec_surf);

    assert(std::abs(Qgen2 - Qcool2 - Qheat2) <= std::abs(Qgen2) * 1e-10);

    //!< ******************************************************************************************************************************************************
    //!< make the hierarchical module
    constexpr size_t ncel11 = 2;
    constexpr size_t ncel22 = 2;
    //!< constexpr size_t ncel33 = 3; -> gives problems with voltage limits in cycleAge so skip for now.
    //!<  this unit test focuses on thermal model which does not depend on different sized modules
    constexpr size_t ncel33 = 2;
    std::string n11 = "H1";
    std::string n22 = "H2";
    std::string n33 = "H3";

    Deep_ptr<StorageUnit> SU1[] = { make<Cell_SPM>(), make<Cell_SPM>() };
    Deep_ptr<StorageUnit> SU2[] = { make<Cell_SPM>(), make<Cell_SPM>() };
    Deep_ptr<StorageUnit> SU3[] = { make<Cell_SPM>(), make<Cell_SPM>() };

    auto cp11 = dynamic_cast<Cell_SPM *>(SU1[0].get());
    auto cp22 = dynamic_cast<Cell_SPM *>(SU1[1].get());
    auto cp33 = dynamic_cast<Cell_SPM *>(SU2[0].get());
    auto cp44 = dynamic_cast<Cell_SPM *>(SU2[1].get());
    auto cp55 = dynamic_cast<Cell_SPM *>(SU3[0].get());
    auto cp66 = dynamic_cast<Cell_SPM *>(SU3[1].get());

    Deep_ptr<StorageUnit> MU[] = { make<Module_s>(n11, T, true, false, ncel11, coolControl, 2),
                                   make<Module_s>(n22, T, true, false, ncel22, coolControl, 2),
                                   make<Module_s>(n33, T, true, false, ncel33, coolControl, 2) };


    auto mp11 = dynamic_cast<Module_s *>(MU[0].get());
    auto mp22 = dynamic_cast<Module_s *>(MU[1].get());
    auto mp33 = dynamic_cast<Module_s *>(MU[2].get());

    mp11->setSUs(SU1, checkCells);
    mp22->setSUs(SU2, checkCells);
    mp33->setSUs(SU3, checkCells);

    constexpr size_t nm = 3;
    std::string n44 = "testCoolSystem_complexModule";

    auto mp44 = make<Module_s>(n44, T, true, true, 6, coolControl, 1);
    mp44->setSUs(MU, checkCells, true);
    double Tini22[6] = { cp11->T(), cp22->T(), cp33->T(), cp44->T(), cp55->T(), cp66->T() };

    //!< Use procedure age to cycle it
    //!< cout<<"Procedure_test start coolsystem test with a complex module for cool control setting "<<coolControl<<endl;
    p.cycleAge(mp44.get(), Ncycle, ncheck, nbal, false, Ccha, Cdis, mp->Vmax(), mp->Vmin()); //!< exclude CV phase, CC only

    double Qgen3, Qcool3, Qheat3;
    //!< check balance of module mp11
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp11->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) <= std::abs(Qgen3) * 1e-10);
    //!< check balance of module mp22
    Qgen3 = cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp22->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) <= std::abs(Qgen3) * 1e-10);
    //!< check balance of module mp33
    //!< Qgen3 = cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat(); //!< total heat generated by cells
    Qgen3 = cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat(); //!< total heat generated by cells
    Qcool3 = mp33->getCoolSystem()->getHeatEvac();                       //!< total heat extracted by the coolsystem from the cells
    //!< Qheat3 = -((Tini22[4] - cp55->T()) * (rho*Cp*L*elec_surf) + (Tini22[5] - cp66->T()) * (rho*Cp*L*elec_surf)+ (Tini22[6] - cp77->T()) * (rho*Cp*L*elec_surf));		//!< total energy in heating up the cells
    Qheat3 = -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) <= std::abs(Qgen3) * 1e-10);
    //!< check balance of the top level module
    Qgen3 = mp11->getCoolSystem()->getHeatEvac() + mp22->getCoolSystem()->getHeatEvac() + mp33->getCoolSystem()->getHeatEvac();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) <= std::abs(Qgen3) * 1e-10);

    //!< check balance of total system
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();
    Qgen3 += cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();
    //!< Qgen3 += cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();
    Qgen3 += cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf));
    //!< Qheat3 += -((Tini22[4] - cp55->T()) * (rho*Cp*L*elec_surf) + (Tini22[5] - cp66->T()) * (rho*Cp*L*elec_surf)+ (Tini22[6] - cp77->T()) * (rho*Cp*L*elec_surf));
    Qheat3 += -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();

    //!< Comparison of cool system performance in the different control strategies: print out the following statement
    //!< cout<<"Total heat balance of coolsystem complex module entire "<<coolControl<<" is Qgen = "<<Qgen3<<", Qheat = "<<Qheat3<<", Qcool = "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl<<flush;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) <= std::abs(Qgen3) * 1e-10);
    assert(Qheat3 >= 0);
  }

  return true;
}

bool test_Procedure_cycleAge_stress()
{
  /*
   * Stress test for CycleAge with parallel modules with large variations
   */

  //!< 		5 widely different cells
  //!< 		4 similar and one aged cell

  DEG_ID deg;
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  double T2 = settings::T_ENV;
  bool checkCells2 = false;
  double Vbal = 3.5;
  int ndata = 20;
  bool balance = true;
  bool unittest = true;

  //!< 5 cells with large distribution
  std::cout << "Start procedure test with a wide distribution in capacity\n";
  constexpr size_t ncel1 = 5;
  double std1 = 0.1;
  double std2 = 0.15;
  std::default_random_engine gen;
  std::normal_distribution<double> distr_c2(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r2(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d2(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate
  std::string n3 = "proctest_mp_largeVariation";

  Deep_ptr<StorageUnit> cs3[ncel1];

  cs3[0] = make<Cell_SPM>("cell5", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen));
  cs3[1] = make<Cell_SPM>("cell6", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen));
  cs3[2] = make<Cell_SPM>("cell7", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen));
  cs3[3] = make<Cell_SPM>("cell8", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen));
  cs3[4] = make<Cell_SPM>("cell9", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen));

  auto mpp3 = make<Module_p>(n3, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage

  mpp3->setSUs(cs3, checkCells2, true);
  auto p1 = Procedure(balance, Vbal, ndata, unittest);
  p1.cycleAge(mpp3.get(), false);

  //!< 4 similar and one very different
  std::cout << "Start procedure test with 4 similar and one very different cell.\n";
  std::string n4 = "proctest_mp_4and1";
  std1 = 0.004;
  std2 = 0.025;
  std::normal_distribution<double> distr_c(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate

  Deep_ptr<StorageUnit> cs4[ncel1];
  cs4[0] = make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[1] = make<Cell_SPM>("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[2] = make<Cell_SPM>("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[3] = make<Cell_SPM>("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[4] = make<Cell_SPM>("cell9", deg, 0.5, 2.0, 1.1, 1.1); //!< one with half the capacity and double the resistance and 10% more degradation

  auto mpp4 = make<Module_p>(n4, T2, true, false, ncel1, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp4->setSUs(cs4, checkCells2, true);
  auto p2 = Procedure(balance, Vbal, ndata, unittest);
  p2.cycleAge(mpp4.get(), false); //!< this should write a file called s_module_capacities.csv

  return true;
}

bool test_degradationModel(bool capsread, bool Rspread, bool degspread, DEG_ID deg, int cool)
{
  /*
   * Test the cycle-ageing procedure with different degradation models to see which produces different knee points.
   * This function cycles a series and parallel module with the given settings
   *
   * IN
   * parameters on cell-to-cell variations and degradation models
   */

  //!< name of files
  Cycler cyc1;
  std::string degid = deg.print();
  std::string ncap = "noCapSpread";
  std::string nR = "noRSpread";
  std::string ndeg = "noDegSpread";

  //!< Settings
  std::default_random_engine gen;
  double std_cap = 0;
  double std_R = 0;
  double std_d = 0;
  if (capsread) {
    std_cap = 0.004;
    ncap = "capSpread";
  }
  if (Rspread) {
    std_R = 0.025;
    nR = "RSpread";
  }
  if (degspread) {
    std_d = 0.05;
    ndeg = "degSpread";
  }
  std::string pref = degid + '_' + ncap + '_' + nR + '_' + ndeg;

  //!< Control printing
  std::cout << "Start simulations for prefix " << pref << '\n';
  bool unitTest = true; //!< if true, less stuff is printed

  std::normal_distribution<double> distr_c(1.0, std_cap); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r(1.0, std_R);   //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d(1.0, std_d);   //!< normal distribution with mean 1 and std 2.5% for cell degradation rate

  //!< Make the procedure with standard settings
  double Vbal = 3.5;
  int ndata = 20;
  bool balance = true;
  double Ccha = 1;
  double Cdis = 1;
  double Vmax;
  double Vmin;
  int Ncycle = 5000; //!< 1 CC cycle takes about 0.5 min
  int ncheck = 500;  //!< do a checkup ever 100 cycles
  int nbal = 10;     //!< balance every 10 cycles

  //!< test with series module
  constexpr size_t ncel = 5;
  std::string n = pref + "_proctest_sMod";

  Deep_ptr<StorageUnit> cs[ncel];

  cs[0] = make<Cell_SPM>("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs[1] = make<Cell_SPM>("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs[2] = make<Cell_SPM>("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs[3] = make<Cell_SPM>("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs[4] = make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));

  double T = 25_degC;
  bool checkCells = false;

  auto mp = make<Module_s>(n, T, true, false, ncel, cool, true); //!< no multithreading

  mp->setSUs(cs, checkCells, true);
  Vmax = mp->Vmax();
  Vmin = mp->Vmin();
  auto p2 = Procedure(balance, Vbal, ndata, unitTest);
  try {
    p2.cycleAge(mp.get(), Ncycle, ncheck, nbal, true, Ccha, Cdis, Vmax, Vmin);
  } catch (...) {
    std::cout << "Something happened.\n";
  }
  //!< 	check that all cells age more or less the same
  //!< 	the capacity of the string should be the capacity of the smallest cell

  std::cout << "Procedure_test start parallel module\n";

  //!< test with parallel module
  std::string n2 = pref + "_proctest_pMod";

  Deep_ptr<StorageUnit> cs2[] = {
    make<Cell_SPM>("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)),
    make<Cell_SPM>("cell9", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen))
  };

  double T2 = 25_degC;
  bool checkCells2 = false;

  auto mpp = make<Module_p>(n2, T2, true, false, std::size(cs2), cool, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage

  mpp->setSUs(cs2, checkCells2, true);
  Vmax = mpp->Vmax();
  Vmin = mpp->Vmin();

  // auto p3 = Procedure(balance, Vbal, ndata, unitTest);

  // try {
  //   p3.cycleAge(mpp.get(), Ncycle, ncheck, nbal, false, Ccha, Cdis, Vmax, Vmin);
  // } catch (...) {
  // };

  return true;
}

bool test_allDegradationModels(int cool)
{
  /*
   * Test the cycle-ageing procedure with different degradation models to see which produces different knee points.
   * This function calls the other one with a range of settings
   */
  bool caps{ true }, Rs{ true }, degs{ true };
  DEG_ID deg{};

  //!< kinetic SEI
  deg.SEI_id.add_model(1);
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< kinetic SEI from paper
  deg.SEI_id[0] = 4;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< kinetic SEI + porosity
  deg.SEI_id[0] = 1;
  deg.SEI_porosity = 1;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< kinetic SEI from paper + porosity
  deg.SEI_id[0] = 4;
  deg.SEI_porosity = 1;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< reset SEI
  deg.SEI_id[0] = 0;
  deg.SEI_porosity = 0;

  //!< DAI + Laresgoiti LAM
  deg.LAM_id.add_model(1);
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< Delacourt LAM [linear with Ah]
  deg.LAM_id[0] = 2;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< Kindermann LAM
  deg.LAM_id[0] = 3;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< Kindermann LAM + SEI
  deg.SEI_id[0] = 1;
  deg.LAM_id[0] = 3;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< Kindermann LAM + SEI + porosity
  deg.SEI_id[0] = 1;
  deg.LAM_id[0] = 3;
  deg.SEI_porosity = 1;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< SEI + Ekstrom cracks
  deg.SEI_id[0] = 1;
  deg.CS_id.add_model(5);
  deg.LAM_id[0] = 0;
  deg.SEI_porosity = 0;
  test_degradationModel(caps, Rs, degs, deg, cool);

  //!< SEI + Laresgoiti cracks (pure laresgoiti, i.e. with the simplified stress model)
  deg.SEI_id[0] = 1;
  deg.CS_id[0] = 1;
  deg.LAM_id[0] = 0;
  deg.SEI_porosity = 0;
  test_degradationModel(caps, Rs, degs, deg, cool);

  return true;
}

int test_all_Procedure()
{
  int cool = 1;
  //!< Test normal procedures, with and without contact resistance and CV phases
  test_Procedure_cycleAge(0, true, cool);
  //!< test with two different values for contact resistance
  test_Procedure_cycleAge(0, false, cool); //!< no contact resistance, no cell-to-cell variation
  test_Procedure_cycleAge(0.001 / 5.0, false, cool);
  test_Procedure_cycleAge(0.001 / 5.0, true, cool); //!< 0.2 mOhm contact resistance, with cell-to-cell variation

  //!< test with large variation of cells in P module
  test_Procedure_cycleAge_stress();

  //!< Test the cooling system
  test_Procedure_CoolSystem();

  //!< Test various degradation models
  test_allDegradationModels(cool); //!< test them all

  //!< test a specific model
  DEG_ID deg;
  deg.SEI_id.add_model(1); //!< kinetic SEI and porosity, with Dai/Laresgoiti LAM
  deg.SEI_porosity = 1;
  deg.LAM_id.add_model(1);
  test_degradationModel(true, true, true, deg, cool);

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Procedure(); }
