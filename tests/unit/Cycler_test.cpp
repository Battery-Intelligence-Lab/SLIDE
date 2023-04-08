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

bool test_CyclerSU(StorageUnit *su, bool testCV)
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
  double lim{ 0.0 }, V1{};
  ThroughputData th{};

  cyc.initialise(su, "Cycler_test");
  double vlim, Ilim, dt{ 1 };
  int ndata = 0;
  double I = su->Cap(); //!< use a 1C rate

  //!< CC charge
  //!< cout<<"\t CC charge"<<endl<<flush;
  vlim = su->Vmax() - lim;
  auto succ = cyc.CC(-I, vlim, TIME_INF, dt, ndata, th); //!< CC charge must do a slight overshoot
  assert(succ == Status::ReachedVoltageLimit);
  assert(NEAR(su->V(), vlim, tol));

  //!< rest for 1h to relax
  //!< cout<<"\t resting 1h after CC charge for SU "<<su->getFullID()<<endl;
  double tlim = 3600;
  succ = cyc.rest(tlim, dt, ndata, th);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  //!< rest for 10h and check the voltage does not change
  //!< cout<<"\t resting 10h after CC charge for SU "<<su->getFullID()<<endl;
  tlim = 10 * 3600;
  V1 = su->V();
  succ = cyc.rest(tlim, dt, ndata, th);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  assert(NEAR(V1, su->V(), tol));


  //!< CC discharge
  //!< cout<<"\t CC discharge"<<endl<<flush;
  vlim = su->Vmin() + lim;
  succ = cyc.CC(I, vlim, TIME_INF, dt, ndata, th); //!< CC discharge must do a slight overshoot
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() <= vlim);

  //!< rest for 1h to relax
  //!< cout<<"\t resting 1h after CC discharge for SU "<<su->getFullID()<<endl;
  tlim = 3600;
  succ = cyc.rest(tlim, dt, ndata, th);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  //!< rest for 10h and check the voltage does not change
  //!< cout<<"\t resting 10h after CC discharge for SU "<<su->getFullID()<<endl;
  tlim = 10 * 3600;
  V1 = su->V();
  succ = cyc.rest(tlim, dt, ndata, th);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(su->I(), 0.0, 1e-10));
  assert(NEAR(V1, su->V(), tol));

  //!< CC cycle for a specified time
  vlim = su->Vmax() + 100;
  tlim = 1000;
  th = {};
  succ = cyc.CC(-I, vlim, tlim, dt, ndata, th);
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(th.Ah(), std::abs(I * tlim / 3600), tol));
  vlim = su->Vmin() - 100;
  tlim = 1000;
  th = {};
  succ = cyc.CC(I, vlim, tlim, dt, ndata, th);
  //!< note: this brings the complex module from SPM cells close to their minimum voltage (because not all sub-modules take the same current the charge-discharge is not symmetric)
  //!< 	so if cells are cooled below 20 degrees, this test fails since we reach the minimum voltage just before the time limit
  assert(succ == Status::ReachedTimeLimit);
  assert(NEAR(th.Ah(), std::abs(I * tlim / 3600), tol));

  if (testCV) {
    //!< CCCV charge
    //!< cout<<"\t CCCV charge"<<endl<<flush;
    vlim = su->Vmax() - lim;
    //!< cout<<"\t \t starting CC phase"<<endl<<flush;
    succ = cyc.CC(-I, vlim, TIME_INF, dt, ndata, th); //!< CC charge must do a slight overshoot
    //!< cout<<"\t \t terminating CC phase"<<endl<<flush;
    assert(succ == Status::ReachedVoltageLimit);
    assert(su->V() - vlim < tol);
    assert(su->V() >= vlim);
    Ilim = 0.1;
    //!< cout<<"\t \t starting CV phase"<<endl<<flush;
    succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th); //!< CV charge must do a slight overshoot
    //!< cout<<"\t \t terminating CV phase with "<<succ<<"\t"<<su->V()<<"\t"<<su->I()<<endl<<flush;
    assert(succ == Status::ReachedVoltageLimit);
    assert(NEAR(su->V(), vlim, tol));
    assert(-su->I() <= Ilim);

    //!< rest for 1h to relax
    //!< cout<<"\t resting 1h after CCCV charge for SU "<<su->getFullID()<<endl;
    tlim = 3600;
    succ = cyc.rest(tlim, dt, ndata, th);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    //!< rest for 10h and check the voltage does not change
    //!< cout<<"\t resting 10h after CCCV charge for SU "<<su->getFullID()<<endl;
    tlim = 10 * 3600;
    V1 = su->V();
    succ = cyc.rest(tlim, dt, ndata, th);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    assert(NEAR(V1, su->V(), tol));

    //!< CCCV discharge
    //!< cout<<"\t CCCV discharge"<<endl<<flush;
    vlim = su->Vmin() + lim;
    succ = cyc.CC(I, vlim, TIME_INF, dt, ndata, th); //!< CC charge must do a slight overshoot
    assert(succ == Status::ReachedVoltageLimit);
    assert(su->V() - vlim < tol);
    assert(su->V() <= vlim);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th); //!< CV discharge must do a slight overshoot
    //!< cout<<"\t \t terminating CV phase with "<<succ<<"\t"<<su->V()<<"\t"<<su->I()<<" and voltage error "<<abs(su->V()-vlim) <<endl<<flush;
    assert(succ == Status::ReachedCurrentLimit);
    assert(NEAR(su->V(), vlim, tol));
    assert(su->I() <= Ilim);

    //!< rest for 1h to relax
    //!< cout<<"\t resting 1h after CCCV discharge for SU "<<su->getFullID()<<endl;
    tlim = 3600;
    succ = cyc.rest(tlim, dt, ndata, th);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    //!< rest for 10h and check the voltage does not change
    //!< cout<<"\t resting 10h after CCCV discharge for SU "<<su->getFullID()<<endl;
    tlim = 10 * 3600;
    V1 = su->V();
    succ = cyc.rest(tlim, dt, ndata, th);
    assert(succ == Status::ReachedTimeLimit);
    assert(NEAR(su->I(), 0.0, 1e-10));
    assert(NEAR(V1, su->V(), tol));

    //!< getCapacity
    //!< cout<<"\t getCapacity 1"<<endl<<flush;
    double Ah, dtime;
    double cap = cyc.testCapacity(Ah, dtime);
    assert(std::abs(cap - su->Cap()) / su->Cap() < 0.15); //!< check we are close to the nominal capacity, i.e. error < 15%
    su->setBlockDegAndTherm(true);
    double cap2 = cyc.testCapacity(Ah, dtime);
    assert(NEAR(cap, cap2, tol));                          //!< get the capacity while ignoring degradation
    cyc.CC(-cap / 10.0, su->Vmin() + 1.4, 3600, 2, 0, th); //!< charge to get in the valid voltage region
    cyc.setDiagnostic(true);
    //  double cap3 = cyc.testCapacity(Ah, dtime);
    //!< cout<<"\t capacity measurements: "<<cap<<"\t"<<cap2<<"\t"<<cap3<<"\t for nominal capacity "<<su->Cap()<<endl; //
    //!< assert(std::abs(cap-cap3) < tol);								//!< get the capacity while stopping when one cell reached the voltage limit
    //!<  this should fail for a complex hierarchial module (see testCyclerSPM)
  }

  return true;
}


template <typename Cell_t>
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
  auto cp1 = make<Cell_t>();
  test_CyclerSU(cp1.get(), checkCV);

  //!< test series of cell
  //!< cout<<"Test cycler made of a series module of Cells"<<endl;
  Deep_ptr<StorageUnit> cs[] = { make<Cell_t>(), make<Cell_t>() };
  std::string n = "module_series_cell";
  auto ms = make<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  ms->setSUs(cs, checkCells, true);
  test_CyclerSU(ms.get(), checkCV);

  //!< test parallel of cells
  //!< cout<<"Test cycler made of a parallel module of Cells"<<endl;
  Deep_ptr<StorageUnit> cs2[] = { make<Cell_t>(), make<Cell_t>() };
  n = "module_paralell_cell";
  auto mp = make<Module_p>(n, T, true, false, std::size(cs2), 1, 1);
  mp->setSUs(cs2, checkCells, true);
  test_CyclerSU(mp.get(), checkCV);

  //!< test complex module
  //!< cout<<"Cycler test complex module of Cells"<<endl;
  std::string ids[] = { "1", "2", "3" };
  Deep_ptr<StorageUnit> SU1[] = { make<Cell_t>(), make<Cell_t>() };
  Deep_ptr<StorageUnit> SU2[] = { make<Cell_t>(), make<Cell_t>() };
  Deep_ptr<StorageUnit> SU3[] = { make<Cell_t>(), make<Cell_t>(), make<Cell_t>() };


  auto mu1 = make<Module_p>(ids[0], T, true, false, std::size(SU1), 1, 2);
  auto mu2 = make<Module_p>(ids[1], T, true, false, std::size(SU2), 1, 2);
  auto mu3 = make<Module_p>(ids[2], T, true, false, std::size(SU3), 1, 2);

  mu1->setSUs(SU1, checkCells);
  mu2->setSUs(SU2, checkCells);
  mu3->setSUs(SU3, checkCells);

  Deep_ptr<StorageUnit> MU[] = {
    std::move(mu1),
    std::move(mu2),
    std::move(mu3),
  };

  checkCells = true;
  auto msp = make<Module_s>("cells_complex", T, true, false, 7, 1, true);
  msp->setSUs(MU, checkCells, true); //!< three module_p in series

  test_CyclerSU(msp.get(), checkCV);
  //!< note the capacity of this module will be larger than the parallel module
  //!< 	getCapacity (dis)charges to the voltage limit of the total module (so in this case 3*cell voltage)
  //!< 	The 2 modules with 2 cells will reach their voltage limit after 2*cell capacity
  //!< 	but the 3rd module has three cells, so it will not have reached its voltage limit yet
  //!< 	so the overall module won't have reached its voltage limit and you keep (dis)charging
  //!< 	and then the overall module reaches its voltage limit after > 2 * cell capacity, the first two will be over(dis)charged
  //!< so the capacity check with and without diagnostic will give a different result

  return true;
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
  std::default_random_engine gen;
  double std1 = 0.004;
  double std2 = 0.025;
  std::normal_distribution<double> distr_c(1.0, std1); //!< normal distribution with mean 1 and std 0.4% for cell capacity
  std::normal_distribution<double> distr_r(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell resistance
  std::normal_distribution<double> distr_d(1.0, std2); //!< normal distribution with mean 1 and std 2.5% for cell degradation rate

  //!< degradation

  DEG_ID deg;

  deg.SEI_id.add_model(4); //!< chirstensen SEI growth
  deg.SEI_porosity = 0;    //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0); //!< no surface cracks
  deg.CS_diffusion = 0;   //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  deg.LAM_id.add_model(0); //!< no LAM
  deg.pl_id = 0;           //!< no litihium plating

  //!< Make a parallel module with 9 cells, and contact resistance Rc
  std::string n2 = "Variations_p_module_" + std::to_string(Rc);

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

  double Rcs2[] = { Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc, Rc };
  double T2 = settings::T_ENV;
  bool checkCells2 = false;
  auto mpp = make<Module_p>(n2, T2, true, false, std::size(cs2), 1, 1); //!< no multithreading, nt_Vcheck time steps between checking SU voltage
  mpp->setSUs(cs2, checkCells2, true);
  mpp->setRcontact(Rcs2);

  //!< call the test function
  test_CyclerSU(mpp.get(), false); //!< don't do CV phases

  return true;
}

bool test_Cycler_writeData(int control)
{
  /*
   * Write the cell data during a CCCV cycle.
   * The results are written to a csv file, and have to be plotted by the Matlab script 'readTestData.m'
   * Then users have to manually verify the results look ok.
   */

  auto su = make<Cell_SPM>();
  Cycler cyc;
  double tol = settings::MODULE_P_V_ABSTOL; //!< complex modules have larger under- and overshoot due to the larger numbers involved
  double lim = 0.0;
  ThroughputData th;
  double Ah{}, Wh{}, dtime{};

  //!< ******************************************************* test a single cell doing a single CCCV cycle *********************************************
  cyc.initialise(su.get(), "Cycler_test_cell");
  double vlim, Ilim;
  double dt = 2;
  int ndata = 2; //!< store data every 2 seconds (or every dt)
  Status succ;
  double I = su->Cap();

  //!< CCCV charge
  vlim = su->Vmax() - lim;
  succ = cyc.CC(-I, vlim, TIME_INF, dt, ndata, th);
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() >= vlim);
  Ilim = 0.1;
  succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th);
  assert(succ == Status::ReachedCurrentLimit);
  assert(NEAR(su->V(), vlim, tol));
  assert(-su->I() <= Ilim);

  //!< CCCV discharge
  vlim = su->Vmin() + lim;
  succ = cyc.CC(I, vlim, TIME_INF, dt, ndata, th);
  assert(succ == Status::ReachedVoltageLimit);
  assert(su->V() - vlim < tol);
  assert(su->V() <= vlim);
  Ilim = 0.1;
  succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th);
  assert(succ == Status::ReachedCurrentLimit);
  assert(NEAR(su->V(), vlim, tol));
  assert(su->I() <= Ilim);

  //!< write the data
  su->writeData("test_writeData");

  //!< ******************************************************* test three cells in a module doing lots of cycles *********************************************
  bool checkCells = false;
  DEG_ID deg;

  deg.SEI_id.add_model(4); //!< chirstensen SEI growth
  deg.SEI_porosity = 0;    //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0); //!< no surface cracks
  deg.CS_diffusion = 0;   //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  deg.LAM_id.add_model(0); //!< no LAM
  deg.pl_id = 0;           //!< no litihium plating
  Deep_ptr<StorageUnit> cs[] = { make<Cell_SPM>("cell1", deg, 1, 1, 1, 1),
                                 make<Cell_SPM>("cell2", deg, 1, 1, 1, 1), //!< check that the middle cell heats up more
                                 make<Cell_SPM>("cell3", deg, 1, 1, 1, 1) };
  std::string n = "mod1";
  auto ms = make<Module_s>(n, settings::T_ENV, true, false, std::size(cs), control, 1);
  ms->setSUs(cs, checkCells, true);
  cyc.initialise(ms.get(), "test_writeData_sModule");

  //!< do 100 cycles
  for (int i = 0; i < 100; i++) {
    //!< CCCV charge
    vlim = ms->Vmax() - lim;
    succ = cyc.CC(-I, vlim, TIME_INF, dt, ndata, th);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th);

    //!< CCCV discharge
    vlim = ms->Vmin() + lim;
    succ = cyc.CC(I, vlim, TIME_INF, dt, ndata, th);
    Ilim = 0.1;
    succ = cyc.CV(vlim, Ilim, TIME_INF, dt, ndata, th);
  }

  //!< Some data might be written during the 100 cycles, push the rest out too (note the prefix must be the same or this last batch will end up in a different file)
  ms->writeData("test_writeData_sModule");

  return true;
}

bool test_Cycler_CoolSystem()
{
  /*
   * test the cool system design with proper cycle ageing
   */

  //!< Parameters from Cell_SPM which are needed to calculate the heat balance
  double rho = 1626;
  double Cp = 750;
  double L = 1.6850e-4;               //!< thickness of one layer
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
  ThroughputData th{};
  double vlim;
  int ndata = 0;

  //!< Loop for each setting of the cool controller
  for (int coolControl = 1; coolControl < 6; coolControl++) {

    //!< ****************************************************************************************************************************************************
    //!< Make a simple module with one SPM cell
    Deep_ptr<StorageUnit> cs[] = { make<Cell_SPM>() };
    auto cp0 = dynamic_cast<Cell_SPM *>(cs[0].get());

    std::string n = "testCoolSystem";
    auto mp = make<Module_s>(n, T, true, false, std::size(cs), coolControl, 1);
    mp->setSUs(cs, checkCells, true);
    double Tini[1] = { cp0->T() };
    cyc.initialise(mp.get(), "Cycler_cooltest_oneCell");

    //!< do a few 1C cycles
    Icha = -cp0->Cap();
    Idis = cp0->Cap();
    for (int i = 0; i < N; i++) {
      //!< charge
      vlim = mp->Vmax() - lim;
      cyc.CC(Icha, vlim, TIME_INF, dt, ndata, th);

      //!< CC discharge
      vlim = mp->Vmin() + lim;
      cyc.CC(Idis, vlim, TIME_INF, dt, ndata, th);
    }

    if (settings::T_MODEL == 2) {
      //!< check the energy balance of the outer module
      double Qgen = cp0->thermal_getTotalHeat();         //!< total heat generated by cells
      double Qcool = mp->getCoolSystem()->getHeatEvac(); //!< total heat extracted by the coolsystem from the cells
      double Tnew[1] = { cp0->T() };
      double Qheat = 0; //!< total energy in heating up the cells
      for (int i = 0; i < 1; i++)
        Qheat += (Tnew[i] - Tini[i]) * (rho * Cp * L * elec_surf);
      //!< cout<<"Total heat balance of coolsystem single cell "<<coolControl<<" is "<<Qgen<<", "<<Qheat<<", "<<Qcool<<" and error "<<abs(Qgen - Qcool - Qheat)<<endl;
      assert(NEAR(Qgen, Qcool + Qheat, 1e-9));
    }

    //!< **********************************************************************************************************************************************************
    //!< Make a simple module with SPM cells
    Deep_ptr<StorageUnit> cs2[] = {
      make<Cell_SPM>(),
      make<Cell_SPM>(),
      make<Cell_SPM>(),
      make<Cell_SPM>()
    };

    auto cp1 = dynamic_cast<Cell_SPM *>(cs2[0].get());
    auto cp2 = dynamic_cast<Cell_SPM *>(cs2[1].get());
    auto cp3 = dynamic_cast<Cell_SPM *>(cs2[2].get());
    auto cp4 = dynamic_cast<Cell_SPM *>(cs2[3].get());

    std::string n2 = "testCoolSystem";
    auto mp2 = make<Module_s>(n2, T, true, false, std::size(cs2), coolControl, 1);
    mp2->setSUs(cs2, checkCells, true);
    double Tini2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };
    cyc.initialise(mp2.get(), "Cycler_cooltest_simpleModule");

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp1->Cap();
    Idis = cp1->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      vlim = mp2->Vmax() - lim;
      cyc.CC(Icha, vlim, TIME_INF, dt, ndata, th);

      //!< CC discharge
      vlim = mp2->Vmin() + lim;
      cyc.CC(Idis, vlim, TIME_INF, dt, ndata, th);
    }

    if (settings::T_MODEL == 2) {
      //!< check the energy balance of the outer module
      double Qgen2 = cp1->thermal_getTotalHeat() + cp2->thermal_getTotalHeat() + cp3->thermal_getTotalHeat() + cp4->thermal_getTotalHeat(); //!< total heat generated by cells
      double Qcool2 = mp2->getCoolSystem()->getHeatEvac();                                                                                  //!< total heat extracted by the coolsystem from the cells
      double Tnew2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };
      double Qheat2 = 0; //!< total energy in heating up the cells
      for (int i = 0; i < 4; i++)
        Qheat2 += (Tnew2[i] - Tini2[i]) * (rho * Cp * L * elec_surf);
      //!< cout<<"Total heat balance of coolsystem simle module "<<coolControl<<" is "<<Qgen2<<", "<<Qheat2<<", "<<Qcool2<<" and error "<<abs(Qgen2 - Qcool2 - Qheat2)<<endl;
      assert(std::abs(Qgen2 - Qcool2 - Qheat2) / std::abs(Qgen2) < 1e-10);
    }

    //!< ******************************************************************************************************************************************************
    //!< make the hierarchical module
    std::string n11 = "H1";
    std::string n22 = "H2";
    std::string n33 = "H3";
    std::array<Deep_ptr<StorageUnit>, 2> SU1{ make<Cell_SPM>(), make<Cell_SPM>() };
    std::array<Deep_ptr<StorageUnit>, 2> SU2{ make<Cell_SPM>(), make<Cell_SPM>() };
    std::array<Deep_ptr<StorageUnit>, 3> SU3{ make<Cell_SPM>(),
                                              make<Cell_SPM>(),
                                              make<Cell_SPM>() };

    auto cp11 = dynamic_cast<Cell_SPM *>(SU1[0].get());
    auto cp22 = dynamic_cast<Cell_SPM *>(SU1[1].get());
    auto cp33 = dynamic_cast<Cell_SPM *>(SU2[0].get());
    auto cp44 = dynamic_cast<Cell_SPM *>(SU2[1].get());
    auto cp55 = dynamic_cast<Cell_SPM *>(SU3[0].get());
    auto cp66 = dynamic_cast<Cell_SPM *>(SU3[1].get());
    auto cp77 = dynamic_cast<Cell_SPM *>(SU3[2].get());

    Deep_ptr<StorageUnit> MU[] = {
      make<Module_s>(n11, T, true, false, SU1.size(), coolControl, 2),
      make<Module_s>(n22, T, true, false, SU2.size(), coolControl, 2),
      make<Module_s>(n33, T, true, false, SU3.size(), coolControl, 2)
    };

    auto mp11 = dynamic_cast<Module_s *>(MU[0].get());
    auto mp22 = dynamic_cast<Module_s *>(MU[1].get());
    auto mp33 = dynamic_cast<Module_s *>(MU[2].get());
    mp11->setSUs(SU1, checkCells);
    mp22->setSUs(SU2, checkCells);
    mp33->setSUs(SU3, checkCells);

    std::string n44 = "H4";
    auto mp44 = make<Module_s>(n44, T, true, true, 7, coolControl, 1);
    mp44->setSUs(MU, checkCells, true);
    double Tini22[7] = { cp11->T(), cp22->T(), cp33->T(), cp44->T(), cp55->T(), cp66->T(), cp77->T() };
    cyc.initialise(mp44.get(), "Cycler_cooltest_complexModule");

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp11->Cap();
    Idis = cp11->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      vlim = mp44->Vmax() - lim;
      cyc.CC(Icha, vlim, TIME_INF, dt, ndata, th);

      //!< CC discharge
      vlim = mp44->Vmin() + lim;
      cyc.CC(Idis, vlim, TIME_INF, dt, ndata, th);
    }

    double Qgen3, Qcool3, Qheat3;
    //!< check balance of module mp11
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp11->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 1 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(NEAR(Qgen3, Qcool3 + Qheat3, 1e-9));
    //!< check balance of module mp22
    Qgen3 = cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp22->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 2 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(NEAR(Qgen3, Qcool3 + Qheat3, 1e-9));
    //!< check balance of module mp33
    Qgen3 = cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();                                                                             //!< total heat generated by cells
    Qcool3 = mp33->getCoolSystem()->getHeatEvac();                                                                                                                                  //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf) + (Tini22[6] - cp77->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    //!< cout<<"Total heat balance of coolsystem complex module 3 "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(NEAR(Qgen3, Qcool3 + Qheat3, 1e-9));

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

  return true;
}

int test_all_Cycler()
{
  auto test_CyclerVariations_0 = []() { return test_CyclerVariations(0.0); };
  auto test_CyclerVariations_high = []() { return test_CyclerVariations(0.001 / 5.0); };
  auto test_Cycler_writeData_1 = []() { return test_Cycler_writeData(1); };

  if (!TEST(test_CyclerCell<Cell_Bucket>, "test_Cycler_Cell_Bucket")) return 1;
  if (!TEST(test_CyclerCell<Cell_ECM<1>>, "test_Cycler_Cell_ECM<1>")) return 2;
  if (!TEST(test_CyclerCell<Cell_SPM>, "test_Cycler_Cell_SPM")) return 3;
  if (!TEST(test_CyclerVariations_0, "test_CyclerVariations_0")) return 4;
  if (!TEST(test_CyclerVariations_high, "test_CyclerVariations_high")) return 5;
  if (!TEST(test_Cycler_writeData_1, "test_Cycler_writeData_1")) return 6; // #TODO passes but takes long!

  if (settings::T_MODEL == 2) // #TODO only valid if T_MODEL==2
    if (!TEST(test_Cycler_CoolSystem, "test_Cycler_CoolSystem")) return 7;

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Cycler(); }