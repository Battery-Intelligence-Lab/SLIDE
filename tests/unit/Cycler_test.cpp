/*
 * Cycler_test.cpp
 *
 *  Created on: 19 Dec 2019
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
bool test_Cycler_SU(std::unique_ptr<StorageUnit> su, bool testCV)
{
  /*
   * define all the tests with one storage unit
   * This function can then be called with a cell, s or p module, etc
   *
   * You always stay 0.2V away from the minimum and maximum voltage limits
   * 	complex module has 3 series modules, one of which has 3 cells and the other 2 cells
   * 	so if you go to the max voltage of the series string, then the two modules with 2 cells will exceed the voltage limits
   * 	while the module with three cells is not yet at the maximum voltage and neither is the total voltage
   */

  Cycler cyc;
  double tol = settings::MODULE_P_V_ABSTOL; //!< complex modules have larger under- and overshoot due to the larger numbers involved
  double lim = 0.0;
  double Ah, Wh, dtime, V1;

  cyc.initialise(su, "Cycler_test");
  double vlim, tlim, Ilim;
  double dt = 2;
  int ndata = 0;
  double I = su->Cap(); //!< use a 1C rate

  //!< CC charge
  //!< cout<<"\t CC charge"<<endl<<flush;
  vlim = su->Vmax() - lim;
  tlim = 99999999;
  auto succ = cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh, dtime); //!< CC charge must do a slight overshoot
  assert(succ == Status::ReachedVoltageLimit);

  assert(su->V() - vlim < tol);
  assert(su->V() >= vlim);

  //!< rest for 1h to relax
  //!< cout<<"\t resting 1h after CC charge for SU "<<su->getFullID()<<endl;
  tlim = 3600;
  succ = cyc.rest(tlim, dt, ndata, Ah, Wh);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  //!< rest for 10h and check the voltage does not change
  //!< cout<<"\t resting 10h after CC charge for SU "<<su->getFullID()<<endl;
  tlim = 10 * 3600;
  V1 = su->V();
  succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));

  //!< CC discharge
  //!< cout<<"\t CC discharge"<<endl<<flush;
  vlim = su->Vmin() + lim;
  tlim = 99999999;
  succ = cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh, dtime); //!< CC discharge must do a slight overshoot
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() <= vlim);

  //!< rest for 1h to relax
  //!< cout<<"\t resting 1h after CC discharge for SU "<<su->getFullID()<<endl;
  tlim = 3600;
  succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  //!< rest for 10h and check the voltage does not change
  //!< cout<<"\t resting 10h after CC discharge for SU "<<su->getFullID()<<endl;
  tlim = 10 * 3600;
  V1 = su->V();
  succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  assert(NEAR(V1, su->V(), tol));

  //!< CC cycle for a specified time
  vlim = su->Vmax() + 100;
  tlim = 1000;
  succ = cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
  assert(succ == Status::ReachedTimeLimit);
  assert(std::abs(Ah + I * tlim / 3600) < tol);
  vlim = su->Vmin() - 100;
  tlim = 1000;
  succ = cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
  //!< note: this brings the complex module from SPM cells close to their minimum voltage (because not all sub-modules take the same current the charge-discharge is not symmetric)
  //!< 	so if cells are cooled below 20 degrees, this test fails since we reach the minimum voltage just before the time limit
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(Ah, I * tlim / 3600, tol));

  if (testCV) {
    //!< CCCV charge
    //!< cout<<"\t CCCV charge"<<endl<<flush;
    vlim = su->Vmax() - lim;
    tlim = 99999999;
    //!< cout<<"\t \t starting CC phase"<<endl<<flush;
    succ = cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh, dtime); //!< CC charge must do a slight overshoot
    //!< cout<<"\t \t terminating CC phase"<<endl<<flush;
    assert(succ == Status::ReachedVoltageLimit);
    assert(su->V() - vlim < tol);
    assert(su->V() >= vlim);
    Ilim = 0.1;
    //!< cout<<"\t \t starting CV phase"<<endl<<flush;
    succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime); //!< CV charge must do a slight overshoot
    //!< cout<<"\t \t terminating CV phase with "<<succ<<"\t"<<su->V()<<"\t"<<su->I()<<endl<<flush;
    assert(succ == Status::ReachedVoltageLimit);
    assert(NEAR(su->V(), vlim, tol));
    assert(-su->I() <= Ilim);

    //!< rest for 1h to relax
    //!< cout<<"\t resting 1h after CCCV charge for SU "<<su->getFullID()<<endl;
    tlim = 3600;
    succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    //!< rest for 10h and check the voltage does not change
    //!< cout<<"\t resting 10h after CCCV charge for SU "<<su->getFullID()<<endl;
    tlim = 10 * 3600;
    V1 = su->V();
    succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    assert(NEAR(V1, su->V(), tol));

    //!< CCCV discharge
    //!< cout<<"\t CCCV discharge"<<endl<<flush;
    vlim = su->Vmin() + lim;
    tlim = 99999999;
    succ = cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh, dtime); //!< CC charge must do a slight overshoot
    assert(succ == Status::ReachedVoltageLimit);
    assert(su->V() - vlim < tol);
    assert(su->V() <= vlim);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime); //!< CV discharge must do a slight overshoot
    //!< cout<<"\t \t terminating CV phase with "<<succ<<"\t"<<su->V()<<"\t"<<su->I()<<" and voltage error "<<abs(su->V()-vlim) <<endl<<flush;
    assert(succ == Status::ReachedCurrentLimit);
    assert(NEAR(su->V(), vlim, tol));
    assert(su->I() <= Ilim);

    //!< rest for 1h to relax
    //!< cout<<"\t resting 1h after CCCV discharge for SU "<<su->getFullID()<<endl;
    tlim = 3600;
    succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    //!< rest for 10h and check the voltage does not change
    //!< cout<<"\t resting 10h after CCCV discharge for SU "<<su->getFullID()<<endl;
    tlim = 10 * 3600;
    V1 = su->V();
    succ = cyc.rest(tlim, dt, ndata, Ah, Ah);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    assert(NEAR(V1, su->V(), tol));

    //!< getCapacity
    //!< cout<<"\t getCapacity 1"<<endl<<flush;
    double cap = cyc.testCapacity(Ah, dtime);
    assert(std::abs(cap - su->Cap()) / su->Cap() < 0.15); //!< check we are close to the nominal capacity, i.e. error < 15%
    su->setBlockDegAndTherm(true);
    double cap2 = cyc.testCapacity(Ah, dtime);
    assert(NEAR(cap, cap2, tol));                                     //!< get the capacity while ignoring degradation
    cyc.CC(-cap / 10.0, su->Vmin() + 1.4, 3600, 2, 0, Ah, Wh, dtime); //!< charge to get in the valid voltage region
    cyc.setDiagnostic(true);
    double cap3 = cyc.testCapacity(Ah, dtime);
    //!< cout<<"\t capacity measurements: "<<cap<<"\t"<<cap2<<"\t"<<cap3<<"\t for nominal capacity "<<su->Cap()<<endl; //
    //!< assert(std::abs(cap-cap3) < tol);								//!< get the capacity while stopping when one cell reached the voltage limit
    //!<  this should fail for a complex hierarchial module (see testCyclerSPM)
  }
}

bool test_CyclerCell()
{
  /*
   * Test a cycler with cells, and modules made out of cells
   */
  double T = settings::T_ENV;
  bool checkCells = false;
  bool checkCV = false; //!< Cells can't do a CV phase since they have no diffusion
                        //!< the resistive effect is much to quick for the PI controller in Cycler

  //!< test cell
  //!< cout<<"Test cycler made of one Cell"<<endl;


  auto cp1 = std::make_unique<Cell>();
  test_Cycler_SU(cp1, checkCV);

  //!< test series of cell
  //!< cout<<"Test cycler made of a series module of Cells"<<endl;
  constexpr int ncel = 2;
  auto cp2 = std::make_unique<Cell>();
  auto cp3 = std::make_unique<Cell>();
  std::unique_ptr<StorageUnit> cs[] = { cp2, cp3 };
  std::string n = "module_series_cell";
  auto ms = std::make_unique<Module_s>(n, T, true, false, ncel, 1, 1);
  ms->setSUs(cs, checkCells, true);
  test_Cycler_SU(ms, checkCV);

  //!< test parallel of cells
  //!< cout<<"Test cycler made of a parallel module of Cells"<<endl;
  auto cp4 = std::make_unique<Cell>();
  auto cp5 = std::make_unique<Cell>();
  std::unique_ptr<StorageUnit> cs2[] = { cp4, cp5 };
  n = "module_paralell_cell";
  auto mp = std::make_unique<Module_p>(n, T, true, false, ncel, 1, 1);
  mp->setSUs(cs2, checkCells, true);
  test_Cycler_SU(mp, checkCV);

  //!< test complex module
  //!< cout<<"Cycler test complex module of Cells"<<endl;
  constexpr int ncel1 = 2;
  constexpr int ncel2 = 2;
  int ncel3 = 3;
  std::string n1 = "1";
  std::string n2 = "2";
  std::string n3 = "3";
  auto cp10 = std::make_unique<Cell>();
  auto cp20 = std::make_unique<Cell>();
  auto cp30 = std::make_unique<Cell>();
  auto cp40 = std::make_unique<Cell>();
  auto cp50 = std::make_unique<Cell>();
  auto cp60 = std::make_unique<Cell>();
  auto cp70 = std::make_unique<Cell>();
  std::unique_ptr<StorageUnit> SU1[] = { cp10, cp20 };
  std::unique_ptr<StorageUnit> SU2[] = { cp30, cp40 };
  std::unique_ptr<StorageUnit> SU3[] = { cp50, cp60, cp70 };
  auto mp1 = std::make_unique<Module_p>(n1, T, true, false, ncel1, 1, 2);
  auto mp2 = std::make_unique<Module_p>(n2, T, true, false, ncel2, 1, 2);
  auto mp3 = std::make_unique<Module_p>(n3, T, true, false, ncel3, 1, 2);
  mp1->setSUs(SU1, ncel1, checkCells);
  mp2->setSUs(SU2, ncel2, checkCells);
  mp3->setSUs(SU3, ncel3, checkCells);
  int nm = 3;
  std::string n4 = "cells_complex";
  std::unique_ptr<StorageUnit> MU[] = { mp1, mp2, mp3 };
  checkCells = true;
  auto msp = std::make_unique<Module_s>(n4, T, true, false, 7, 1, true);
  msp->setSUs(MU, checkCells, true); //!< three module_p in series
  test_Cycler_SU(msp, checkCV);
}
bool test_CyclerECM()
{
  /*
   * Test a cycler with ECM cells, and modules made out of ECM cells
   */
  double T = settings::T_ENV;
  bool checkCells = false;
  bool checkCV = false; //!< Cells can't do a CV phase since they have no diffusion
                        //!< the resistive effect is much to quick for the PI controller in Cycler
  int ncel;
  std::string n;

  //!< test cell
  //!< cout<<"Test cycler made of one ECM Cell"<<endl;
  auto cp1 = std::make_unique<Cell_ECM>();
  test_Cycler_SU(cp1, checkCV);

  //!< test series of cell
  //!< cout<<"Test cycler made of a series module of Cells"<<endl;
  ncel = 2;
  auto cp2 = std::make_unique<Cell_ECM>();
  auto cp3 = std::make_unique<Cell_ECM>();
  std::unique_ptr<StorageUnit> cs[] = { cp2, cp3 };
  n = "module_series_ECMcell";
  auto ms = std::make_unique<Module_s>(n, T, true, false, ncel, 1, 1);
  ms->setSUs(cs, checkCells, true);
  test_Cycler_SU(ms, checkCV);

  //!< test parallel of cells
  //!< cout<<"Test cycler made of a parallel module of Cells"<<endl;
  ncel = 2;
  auto cp4 = std::make_unique<Cell_ECM>();
  auto cp5 = std::make_unique<Cell_ECM>();
  std::unique_ptr<StorageUnit> cs2[] = { cp4, cp5 };
  n = "module_parallel_ECMcell";
  auto mp = std::make_unique<Module_p>(n, T, true, false, ncel, 1, 1);
  mp->setSUs(cs2, checkCells, true);
  test_Cycler_SU(mp, checkCV);

  //!< test complex module
  //!< cout<<"Cycler test complex module of Cells"<<endl;
  constexpr int ncel1 = 2;
  constexpr int ncel2 = 2;
  int ncel3 = 3;
  std::string n1 = "1";
  std::string n2 = "2";
  std::string n3 = "3";
  auto cp10 = std::make_unique<Cell_ECM>();
  auto cp20 = std::make_unique<Cell_ECM>();
  auto cp30 = std::make_unique<Cell_ECM>();
  auto cp40 = std::make_unique<Cell_ECM>();
  auto cp50 = std::make_unique<Cell_ECM>();
  auto cp60 = std::make_unique<Cell_ECM>();
  auto cp70 = std::make_unique<Cell_ECM>();
  std::unique_ptr<StorageUnit> SU1[] = { cp10, cp20 };
  std::unique_ptr<StorageUnit> SU2[] = { cp30, cp40 };
  std::unique_ptr<StorageUnit> SU3[] = { cp50, cp60, cp70 };
  auto mp1 = std::make_unique<Module_p>(n1, T, true, false, ncel1, 1, 2);
  auto mp2 = std::make_unique<Module_p>(n2, T, true, false, ncel2, 1, 2);
  auto mp3 = std::make_unique<Module_p>(n3, T, true, false, ncel3, 1, 2);
  mp1->setSUs(SU1, ncel1, checkCells);
  mp2->setSUs(SU2, ncel2, checkCells);
  mp3->setSUs(SU3, ncel3, checkCells);
  int nm = 3;
  std::string n4 = "ECMcells_complex";
  std::unique_ptr<StorageUnit> MU[] = { mp1, mp2, mp3 };
  checkCells = true;
  auto msp = std::make_unique<Module_s>(n4, T, true, false, 7, 1, 1);
  msp->setSUs(MU, checkCells, true); //!< three module_p in series
  test_Cycler_SU(msp, checkCV);
}

bool test_CyclerSPM()
{
  /*
   * Test a cycler with SPM cells, s and p modules out of SPM cells and complex modules
   */
  double T = settings::T_ENV;
  bool checkCells = false;
  bool checkCV = true;
  int ncel;
  std::string n;

  //!< test cell
  //!< cout<<"Test cycler made of one SPM Cell"<<endl;
  auto cp1 = std::make_unique<Cell_SPM>();
  test_Cycler_SU(cp1, checkCV);

  //!< test series of cell
  //!< cout<<"Test cycler made of a series module of SPM Cells"<<endl;
  ncel = 2;
  auto cp2 = std::make_unique<Cell_SPM>();
  auto cp3 = std::make_unique<Cell_SPM>();
  std::unique_ptr<StorageUnit> cs[] = { cp2, cp3 };
  n = "module_series_SPMcell";
  auto ms = std::make_unique<Module_s>(n, T, true, false, ncel, 1, 1);
  ms->setSUs(cs, checkCells, true);
  test_Cycler_SU(ms, checkCV);

  //!< test parallel of cells
  //!< cout<<"Test cycler made of a parallel module of SPM Cells"<<endl;
  auto cp4 = std::make_unique<Cell_SPM>();
  auto cp5 = std::make_unique<Cell_SPM>();
  ncel = 2;
  std::unique_ptr<StorageUnit> cs2[] = { cp4, cp5 };
  n = "module_parallel_SPMcell";
  auto mp = std::make_unique<Module_p>(n, T, true, false, ncel, 1, 1);
  mp->setSUs(cs2, checkCells, true);
  test_Cycler_SU(mp, checkCV);

  //!< test complex module
  //!< cout<<"Cycler test complex module of SPM Cells"<<endl;
  constexpr int ncel1 = 2;
  constexpr int ncel2 = 2;
  int ncel3 = 3;
  std::string n1 = "1";
  std::string n2 = "2";
  std::string n3 = "3";
  auto cp10 = std::make_unique<Cell_SPM>();
  auto cp20 = std::make_unique<Cell_SPM>();
  auto cp30 = std::make_unique<Cell_SPM>();
  auto cp40 = std::make_unique<Cell_SPM>();
  auto cp50 = std::make_unique<Cell_SPM>();
  auto cp60 = std::make_unique<Cell_SPM>();
  auto cp70 = std::make_unique<Cell_SPM>();
  std::unique_ptr<StorageUnit> SU1[] = { cp10, cp20 };
  std::unique_ptr<StorageUnit> SU2[] = { cp30, cp40 };
  std::unique_ptr<StorageUnit> SU3[] = { cp50, cp60, cp70 };
  auto mp1 = std::make_unique<Module_p>(n1, T, true, false, ncel1, 1, 2); //!< pass through coolsystem
  auto mp2 = std::make_unique<Module_p>(n2, T, true, false, ncel2, 1, 2);
  auto mp3 = std::make_unique<Module_p>(n3, T, true, false, ncel3, 1, 2);
  mp1->setSUs(SU1, ncel1, checkCells);
  mp2->setSUs(SU2, ncel2, checkCells);
  mp3->setSUs(SU3, ncel3, checkCells);
  int nm = 3;
  std::string n4 = "SPMcells_complex";
  std::unique_ptr<StorageUnit> MU[] = { mp1, mp2, mp3 };
  checkCells = true;
  auto msp = std::make_unique<Module_s>(n4, T, true, false, 7, 1, 1); //!< top level coolsystem
  msp->setSUs(MU, checkCells, true);                                  //!< three module_p in series
  test_Cycler_SU(msp, checkCV);
  //!< note the capacity of this module will be larger than the parallel module
  //!< 	getCapacity (dis)charges to the voltage limit of the total module (so in this case 3*cell voltage)
  //!< 	The 2 modules with 2 cells will reach their voltage limit after 2*cell capacity
  //!< 	but the 3rd module has three cells, so it will not have reached its voltage limit yet
  //!< 	so the overall module won't have reached its voltage limit and you keep (dis)charging
  //!< 	and then the overall module reaches its voltage limit after > 2 * cell capacity, the first two will be over(dis)charged
  //!< so the capacity check with and without diagnostic will give a different result
}

bool test_CyclerVariations(double Rc)
{
  /*
   * Test modules with cell to cell variations and contact resistance (using SPM cells)
   *
   * IN
   * Rc 	value of the contact resistance
   */

  //!< cout<<"start test with cell-to-cell variations with contact resistance "<<Rc<<endl;

  //!< cell-to-cell variations
  default_random_engine gen;
  double std1 = 0.004;
  double std2 = 0.025;
  normal_distribution<double> distr_c(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  normal_distribution<double> distr_r(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  normal_distribution<double> distr_d(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate

  //!< degradation
  DEG_ID deg;
  deg.SEI_n = 1;        //!< there is 1 SEI model
  deg.SEI_id[0] = 4;    //!< chirstensen SEI growth
  deg.SEI_porosity = 0; //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)
  deg.CS_n = 1;         //!< there is 1 model (that there are no cracks)
  deg.CS_id[0] = 0;     //!< no surface cracks
  deg.CS_diffusion = 0; //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)
  deg.LAM_n = 1;        //!< there are 1 LAM model
  deg.LAM_id[0] = 0;    //!< no LAM
  deg.pl_id = 0;        //!< no litihium plating

  //!< Make a parallel module with 9 cells, and contact resistance Rc
  int ncel2 = 9;
  std::string n2 = "Variations_p_module_" + std::to_string(Rc);
  std::unique_ptr<Cell_SPM> cp7(new Cell_SPM("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp8(new Cell_SPM("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp9(new Cell_SPM("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp10(new Cell_SPM("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp11(new Cell_SPM("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp12(new Cell_SPM("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp13(new Cell_SPM("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp14(new Cell_SPM("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<Cell_SPM> cp15(new Cell_SPM("cell9", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
  std::unique_ptr<StorageUnit> cs2[] = { cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15 };
  double Rcs2[ncel2] = { Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc };
  double T2 = settings::T_ENV;
  bool checkCells2 = false;
  auto mpp = std::make_unique<Module_p>(n2, T2, true, false, ncel2, 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp->setSUs(cs2, checkCells2, true);
  mpp->setRcontact(Rcs2, ncel2);

  //!< call the test function
  test_Cycler_SU(mpp, false); //!< don't do CV phases
}

bool test_Cycler_writeData(int control)
{
  /*
   * Write the cell data during a CCCV cycle.
   * The results are written to a csv file, and have to be plotted by the Matlab script 'readTestData.m'
   * Then users have to manually verify the results look ok.
   */

  auto su = std::make_unique<Cell_SPM>();
  Cycler cyc;
  double tol = settings::MODULE_P_V_ABSTOL; //!< complex modules have larger under- and overshoot due to the larger numbers involved
  double lim = 0.0;
  double Ah, Wh;

  //!< ******************************************************* test a single cell doing a single CCCV cycle *********************************************
  cyc.initialise(su, "Cycler_test_cell");
  double vlim, tlim, Ilim;
  double dt = 2;
  int ndata = 2; //!< store data every 2 seconds (or every dt)
  int succ;
  double I = su->Cap();

  //!< CCCV charge
  vlim = su->Vmax() - lim;
  tlim = 99999999;
  succ = cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() >= vlim);
  Ilim = 0.1;
  succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime);
  assert(succ == Status::ReachedCurrentLimit);
  assert(NEAR(su->V(), vlim, tol));
  assert(-su->I() <= Ilim);

  //!< CCCV discharge
  vlim = su->Vmin() + lim;
  tlim = 99999999;
  succ = cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() <= vlim);
  Ilim = 0.1;
  succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime);
  assert(succ == Status::ReachedCurrentLimit);
  assert(NEAR(su->V(), vlim, tol));
  assert(su->I() <= Ilim);

  //!< write the data
  if (DATASTORE_CELL == 0)
    std::cerr << "Warning in testCycler, we want to write data to test the cell's behaviour but the global variable CYCLER_STORE_CELL is 0 so no data will be stored\n";
  su->writeData("test_writeData");

  //!< ******************************************************* test three cells in a module doing lots of cycles *********************************************
  int ncel = 3;
  bool checkCells = false;
  DEG_ID deg;
  deg.SEI_n = 1;        //!< there is 1 SEI model
  deg.SEI_id[0] = 4;    //!< chirstensen SEI growth
  deg.SEI_porosity = 0; //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)
  deg.CS_n = 1;         //!< there is 1 model (that there are no cracks)
  deg.CS_id[0] = 0;     //!< no surface cracks
  deg.CS_diffusion = 0; //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)
  deg.LAM_n = 1;        //!< there are 1 LAM model
  deg.LAM_id[0] = 0;    //!< no LAM
  deg.pl_id = 0;        //!< no litihium plating
  auto cp2 = std::make_unique<Cell_SPM>("cell1", deg, 1, 1, 1, 1);
  auto cp3 = std::make_unique<Cell_SPM>("cell2", deg, 1, 1, 1, 1); //!< check that the middle cell heats up more
  auto cp4 = std::make_unique<Cell_SPM>("cell3", deg, 1, 1, 1, 1);
  std::unique_ptr<StorageUnit> cs[] = { cp2, cp3, cp4 };
  std::string n = "mod1";
  auto ms = std::make_unique<Module_s>(n, settings::T_ENV, true, false, ncel, control, 1);
  ms->setSUs(cs, checkCells, true);
  cyc.initialise(ms, "test_writeData_sModule");

  //!< do 100 cycles
  for (int i = 0; i < 100; i++) {
    //!< CCCV charge
    vlim = ms->Vmax() - lim;
    tlim = 99999999;
    succ = cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime);

    //!< CCCV discharge
    vlim = ms->Vmin() + lim;
    tlim = 99999999;
    succ = cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh, dtime);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, tlim, dt, ndata, Ah, Wh, dtime);
  }

  //!< Some data might be written during the 100 cycles, push the rest out too (note the prefix must be the same or this last batch will end up in a different file)
  ms->writeData("test_writeData_sModule");
}

bool test_Cycler_CoolSystem()
{
  /*
   * test the cool system design with proper cycle ageing
   */

  //!< Parameters from Cell_SPM which are needed to calculate the heat balance
  double rho = 1626;
  double Cp = 750;
  double L = 1.6850 * pow(10, -4);    //!< thickness of one layer
  double width = 0.1;                 //!< width of the pouch
  double height = 0.2;                //!< height of the pouch
  int nlayers = 31;                   //!< number of layers in the pouch
  double Acell = width * height;      //!< geometric surface area of the pouch
  double elec_surf = Acell * nlayers; //!< total 'unrolled' surface area of the electrodes

  //!< General settings
  double T = settings::T_ENV;
  bool checkCells = true;
  double Icha, Idis;
  double dt = 2;
  int N = 10;
  Cycler cyc;
  double lim = 0.0;
  double Ah, Wh;
  double vlim, tlim;
  int ndata = 0;

  //!< Loop for each setting of the cool controller
  for (int coolControl = 1; coolControl < 6; coolControl++) {

    //!< ****************************************************************************************************************************************************
    //!< Make a simple module with one SPM cell
    int ncel = 1;
    auto cp0 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> cs[] = { cp0 };
    std::string n = "testCoolSystem";
    auto mp = std::make_unique<Module_s>(n, T, true, false, ncel, coolControl, 1);
    mp->setSUs(cs, checkCells, true);
    double Tini[1] = { cp0->T() };
    cyc.initialise(mp, "Cycler_cooltest_oneCell");

    //!< do a few 1C cycles
    Icha = -cp0->Cap();
    Idis = cp0->Cap();
    for (int i = 0; i < N; i++) {
      //!< charge
      vlim = mp->Vmax() - lim;
      tlim = 99999999;
      cyc.CC(Icha, vlim, tlim, dt, ndata, Ah, Wh);

      //!< CC discharge
      vlim = mp->Vmin() + lim;
      cyc.CC(Idis, vlim, tlim, dt, ndata, Ah, Wh);
    }

    //!< check the energy balance of the outer module
    double Qgen = cp0->thermal_getTotalHeat();         //!< total heat generated by cells
    double Qcool = mp->getCoolSystem()->getHeatEvac(); //!< total heat extracted by the coolsystem from the cells
    double Tnew[1] = { cp0->T() };
    double Qheat = 0; //!< total energy in heating up the cells
    for (int i = 0; i < 1; i++)
      Qheat += (Tnew[i] - Tini[i]) * (rho * Cp * L * elec_surf);
    //!< cout<<"Total heat balance of coolsystem single cell "<<coolControl<<" is "<<Qgen<<", "<<Qheat<<", "<<Qcool<<" and error "<<abs(Qgen - Qcool - Qheat)<<endl;
    assert(std::abs(Qgen - Qcool - Qheat) / std::abs(Qgen) < 1e-10);

    //!< **********************************************************************************************************************************************************
    //!< Make a simple module with SPM cells
    int ncel2 = 4;
    auto cp1 = std::make_unique<Cell_SPM>();
    auto cp2 = std::make_unique<Cell_SPM>();
    auto cp3 = std::make_unique<Cell_SPM>();
    auto cp4 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> cs2[] = { cp1, cp2, cp3, cp4 };
    std::string n2 = "testCoolSystem";
    auto mp2 = std::make_unique<Module_s>(n2, T, true, false, ncel2, coolControl, 1);
    mp2->setSUs(cs2, checkCells, true);
    double Tini2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };
    cyc.initialise(mp2, "Cycler_cooltest_simpleModule");

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp1->Cap();
    Idis = cp1->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      vlim = mp2->Vmax() - lim;
      tlim = 99999999;
      cyc.CC(Icha, vlim, tlim, dt, ndata, Ah, Wh);

      //!< CC discharge
      vlim = mp2->Vmin() + lim;
      cyc.CC(Idis, vlim, tlim, dt, ndata, Ah, Wh);
    }

    //!< check the energy balance of the outer module
    double Qgen2 = cp1->thermal_getTotalHeat() + cp2->thermal_getTotalHeat() + cp3->thermal_getTotalHeat() + cp4->thermal_getTotalHeat(); //!< total heat generated by cells
    double Qcool2 = mp2->getCoolSystem()->getHeatEvac();                                                                                  //!< total heat extracted by the coolsystem from the cells
    double Tnew2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };
    double Qheat2 = 0; //!< total energy in heating up the cells
    for (int i = 0; i < 4; i++)
      Qheat2 += (Tnew2[i] - Tini2[i]) * (rho * Cp * L * elec_surf);
    //!< cout<<"Total heat balance of coolsystem simle module "<<coolControl<<" is "<<Qgen2<<", "<<Qheat2<<", "<<Qcool2<<" and error "<<abs(Qgen2 - Qcool2 - Qheat2)<<endl;
    assert(std::abs(Qgen2 - Qcool2 - Qheat2) / std::abs(Qgen2) < 1e-10);

    //!< ******************************************************************************************************************************************************
    //!< make the hierarchical module
    int ncel11 = 2;
    int ncel22 = 2;
    int ncel33 = 3;
    std::string n11 = "H1";
    std::string n22 = "H2";
    std::string n33 = "H3";
    auto cp11 = std::make_unique<Cell_SPM>();
    auto cp22 = std::make_unique<Cell_SPM>();
    auto cp33 = std::make_unique<Cell_SPM>();
    auto cp44 = std::make_unique<Cell_SPM>();
    auto cp55 = std::make_unique<Cell_SPM>();
    auto cp66 = std::make_unique<Cell_SPM>();
    auto cp77 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> SU1[] = { cp11, cp22 };
    std::unique_ptr<StorageUnit> SU2[] = { cp33, cp44 };
    std::unique_ptr<StorageUnit> SU3[] = { cp55, cp66, cp77 };
    auto mp11 = std::make_unique<Module_s>(n11, T, true, false, ncel11, coolControl, 2);
    auto mp22 = std::make_unique<Module_s>(n22, T, true, false, ncel22, coolControl, 2);
    auto mp33 = std::make_unique<Module_s>(n33, T, true, false, ncel33, coolControl, 2);
    mp11->setSUs(SU1, ncel11, checkCells);
    mp22->setSUs(SU2, ncel22, checkCells);
    mp33->setSUs(SU3, ncel33, checkCells);
    int nm = 3;
    std::string n44 = "H4";
    std::unique_ptr<StorageUnit> MU[] = { mp11, mp22, mp33 };
    auto mp44 = std::make_unique<Module_s>(n44, T, true, true, 7, coolControl, 1);
    mp44->setSUs(MU, checkCells, true);
    double Tini22[7] = { cp11->T(), cp22->T(), cp33->T(), cp44->T(), cp55->T(), cp66->T(), cp77->T() };
    cyc.initialise(mp44, "Cycler_cooltest_complexModule");

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp11->Cap();
    Idis = cp11->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      vlim = mp44->Vmax() - lim;
      tlim = 99999999;
      cyc.CC(Icha, vlim, tlim, dt, ndata, Ah, Wh);

      //!< CC discharge
      vlim = mp44->Vmin() + lim;
      cyc.CC(Idis, vlim, tlim, dt, ndata, Ah, Wh);
    }

    double Qgen3, Qcool3, Qheat3;
    //!< check balance of module mp11
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp11->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 1 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);
    //!< check balance of module mp22
    Qgen3 = cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp22->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 2 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);
    //!< check balance of module mp33
    Qgen3 = cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();                                                                             //!< total heat generated by cells
    Qcool3 = mp33->getCoolSystem()->getHeatEvac();                                                                                                                                  //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf) + (Tini22[6] - cp77->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 3 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);

    //!< check balance of the top level module
    Qgen3 = mp11->getCoolSystem()->getHeatEvac() + mp22->getCoolSystem()->getHeatEvac() + mp33->getCoolSystem()->getHeatEvac();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();
    //!< cout<<"Total heat balance of coolsystem complex module top "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);

    //!< check balance of total system
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();
    Qgen3 += cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();
    Qgen3 += cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf) + (Tini22[6] - cp77->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);
    assert(Qheat3 > 0);

    //!< Comparison of cool system performance in the different control strategies: print out the following statement
    //!< cout<<"Total heat balance of coolsystem complex module entire "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
  }
}

int test_all_Cycler()
{
  auto test_CyclerVariations_0 = []() { return test_CyclerVariations(0.0); };
  auto test_CyclerVariations_high = []() { return test_CyclerVariations(0.001 / 5.0); };
  auto test_Cycler_writeData_1 = []() { return test_Cycler_writeData(1); };

  if (!TEST(test_CyclerCell, "test_CyclerCell")) return 1;
  if (!TEST(test_CyclerECM, "test_CyclerECM")) return 2;
  if (!TEST(test_CyclerSPM, "test_CyclerSPM")) return 3;
  if (!TEST(test_CyclerVariations_0, "test_CyclerVariations_0")) return 4;
  if (!TEST(test_CyclerVariations_high, "test_CyclerVariations_high")) return 5;
  if (!TEST(test_Cycler_writeData_1, "test_Cycler_writeData_1")) return 6;
  if (!TEST(test_Cycler_CoolSystem, "test_Cycler_CoolSystem")) return 7;

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Cycler(); }