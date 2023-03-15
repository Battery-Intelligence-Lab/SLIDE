/*
 * paperCode.cpp
 *
 * Created on: 9 Jul 2020
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "slide.hpp"
#include "paperCode.hpp"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <fstream>

namespace slide::paperCode::paper2022 {
void Vequalisation_Rdc(double Rdc)
{
  /*
   * Function to do the simulation showing how the voltage is equalised in a parallel module of 5 cells
   * 4 of which are relatively similar and one is very different
   *
   * IN
   * Rdc 	contact resistance
   */

  std::string ID = "paper_Vequalisation_" + std::to_string(Rdc);

  //!< Make the module with 0 contact resistance
  std::string name = "pmod";
  unsigned seed = 2;
  std::default_random_engine gen(seed);
  std::normal_distribution<double> distr_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distr_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distr_d(1.0, 0.100); //!< normal distribution with mean 1 and std 10%
  slide::DEG_ID deg;
  deg.SEI_id.add_model(4); //!< add model 4.
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0); //!< add model 0.
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  Deep_ptr<StorageUnit> cs4[5];
  cs4[0] = make<Cell_SPM>("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[1] = make<Cell_SPM>("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[2] = make<Cell_SPM>("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[3] = make<Cell_SPM>("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[4] = make<Cell_SPM>("cell5", deg, 0.5, 2.0, 1.1, 1.1); //!< one with half the capacity and double the resistance and 10% more degradation

  double Rcs[5] = { Rdc, Rdc, Rdc, Rdc, Rdc };

  Module_p mpp4(name, settings::T_ENV, true, false, 5, 1, 1); //!< print, no multithread, coolcontrol 1, HVAC coolsystem

  mpp4.setSUs(cs4, true, true);
  mpp4.setRcontact(Rcs);

  //!< Make the cycler
  Cycler cyc;
  cyc.initialise(&mpp4, ID); //!< #TODO we may make this shared pointer. Otherwise it will cause problems. Or unique.
  double vlim;
  double dt = 2;
  int ndata = 2; //!< store data every 2 seconds (or every dt)

  //!< do a CC charge-discharge
  ThroughputData th{};

  double I = mpp4.Cap();
  vlim = mpp4.Vmax();
  cyc.CC(-I, vlim, TIME_INF, dt, ndata, th);

  //!< CC discharge
  vlim = mpp4.Vmin();
  cyc.CC(I, vlim, TIME_INF, dt, ndata, th);

  //!< write the data
  mpp4.writeData(ID);
}

void Vequalisation()
{
  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeTimeData)
    std::cerr << "Warning in paperCode::Vequalisation, we want to write data to "
                 "test the cell's behaviour but the global variable CYCLER_STORE_CELL "
                 "is not 2 so no data will be stored.\n";

  if constexpr (settings::T_MODEL > 1)
    std::cerr << "Warning in paperCode::Vequalisation. We want to use a large contact "
                 "R which will lead to overheating if you use a full thermal model. "
                 "Set settings::T_MODEL to 0 or 1 to enable this simulation.\n";

  //!< Simulate with no contact R and with a value of 1 mOhm
  Vequalisation_Rdc(0);
  Vequalisation_Rdc(0.001);
}

void thermalModel()
{
  using namespace slide;
  //!< calculate the cell temperatures in a small module
  //!< use exactly the same code as for V equalisation, but with Rcontact = 0 and the global setting correctly

  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeTimeData)
    std::cerr << "Warning in paperCode::thermalModel, we want to write data to "
                 "test the cell's behaviour but the global variable CYCLER_STORE_CELL "
                 "is not 2 so no data will be stored.\n";

  if constexpr (settings::T_MODEL != 2)
    std::cerr << "Warning in paperCode::thermalModel. We want to calculate the thermal "
                 "model but the global settings are wrong. Set settings::T_MODEL "
                 "to 2 to enable this simulation.\n";

  std::string ID = "paper_thermalModel";

  //!< Make the module with 0 contact resistance
  std::string name = "pmod";
  unsigned seed = 2;
  std::default_random_engine gen(seed);
  std::normal_distribution<double> distr_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distr_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distr_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  slide::DEG_ID deg;

  deg.SEI_id.add_model(4); //!< add model 4.
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0); //!< add model 0.
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  Deep_ptr<StorageUnit> cs4[5];
  cs4[0] = make<Cell_SPM>("cell1", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[1] = make<Cell_SPM>("cell2", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[2] = make<Cell_SPM>("cell3", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[3] = make<Cell_SPM>("cell4", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen));
  cs4[4] = make<Cell_SPM>("cell5", deg, 0.5, 2.0, 1.1, 1.1); //!< one with half the capacity and double the resistance and 10% more degradation

  auto mpp4 = make<Module_p>(name, settings::T_ENV, true, false, 5, 1, 1); //!< print, no multithread, coolcontrol 1, HVAC coolsystem
  mpp4->setSUs(cs4, true, true);

  //!< Make the cycler
  Cycler cyc;
  cyc.initialise(mpp4.get(), ID);
  double vlim;
  double dt = 2;
  int ndata = 2; //!< store data every 2 seconds (or every dt)


  for (int i = 0; i < 5; i++) {
    //!< do a CC charge-discharge
    ThroughputData th{};
    double I = mpp4->Cap();
    vlim = mpp4->Vmax();
    cyc.CC(-I, vlim, TIME_INF, dt, ndata, th);

    //!< CC discharge
    vlim = mpp4->Vmin();
    cyc.CC(I, vlim, TIME_INF, dt, ndata, th);
  }

  //!< write the data
  mpp4->writeData(ID);
}

void degrade(auto &&su)
{
  /*
   * Test some cycle ageing for an arbitrarily large battery
   *
   */

  //!< standard settings for degradation simulation
  double Ccha = 1;
  double Cdis = 1;
  bool testCV = false;
  double Vmax = su.Vmax();
  double Vmin = su.Vmin();
  int Ncycle = 10000;
  int ncheck = 250; //!< do a checkup ever 250 cycles
  int nbal = 10;    //!< balance every 10 cycles

  //!< Make the procedure
  bool balance = true;
  double Vbal = 3.5; //!< rebalance to 3.1V
  int ndata = 10;    //!< how often to store data (if 0, nothing is stored, not even usage stats)

  Procedure p(balance, Vbal, ndata, false);

  //!< Simulate the cycle ageing
  p.cycleAge(&su, Ncycle, ncheck, nbal, testCV, Ccha, Cdis, Vmax, Vmin);
}

void degradation_1cell()
{
  /*
   * Simulate degradation for a large battery where one cell's behaviour is simulated and then multiplied by the number of series and parallel cells
   */

  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeHistogramData)
    std::cerr << "Warning in paperCode::degradation_1cell, we want to do a "
                 "degradation simulation but the global variable CYCLER_STORE_CELL "
                 "is not 1 so no data or too much data will be stored.\n";

  if constexpr (settings::T_MODEL > 0)
    std::cerr << "Warning in paperCode::degradation_1cell. We want to simulate one "
                 "cell without thermal model. Set settings::T_MODEL to 0 to "
                 "enable this simulation.\n";

  constexpr int ser = 15 * 20;
  constexpr int par = 9 * 7;

  //!< Make one cell with standard parameters//!<Degradation settings
  slide::DEG_ID deg; //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM

  deg.SEI_id.add_model(4); //!< add model 4.
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0); //!< add model 0.
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  auto mc = make<Cell_SPM>("cell1", deg, 1, 1, 1, 1); //!< cell-to-cell variation parameters are 0

  //!< Wrap the cell in a series-module
  auto modulei = make<Module_s>("module1", settings::T_ENV, true, false, 1, 1, 0); //!< single-threaded, conventional coolsystem
                                                                                   //!<(std::string IDi, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype)
  Deep_ptr<StorageUnit> RinB[1] = { std::move(mc) };
  modulei->setSUs(RinB, true, true);

  //!< Wrap the module in a Battery and set number of series and parallel
  auto bat = make<Battery>("single_cell_multiplied_by_300s_63p");
  bat->setSeriesandParallel(ser, par);
  bat->setModule(std::move(modulei));

  //!< call the degradation function
  degrade(*bat); //!< #TODO
}

void degradation_electricalModel()
{
  /*
   * Simulate a battery with identical cells but contact resistances
   */

  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeHistogramData)
    std::cerr << "Warning in paperCode::degradation_electricalModel, we want to "
                 "do a degradation simulation but the global variable CYCLER_STORE_CELL "
                 "is not 1 so no data or too much data will be stored.\n";

  if constexpr (settings::T_MODEL > 0)
    std::cerr << "Warning in paperCode::degradation_electricalModel. We want to simulate "
                 "N cells without thermal model. Set settings::T_MODEL to 0 to "
                 "enable this simulation.\n";

  //!< Get the battery from makeBattery
  int coolControl = 1; //!< just for naming since there is no thermal model
  bool capspread = false;
  bool Rcellspread = false;
  bool degratespread = false;
  bool contactR = true;
  std::string IDaddition = "_noThermalModel";
  //!< todo	std::shared_ptr<StorageUnit> su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR,coolControl,IDaddition,1);

  //!< degrade the battery
  //!< todo	degrade(su);

  //!< Sensitivty with 10 times higher Rcontact
  IDaddition = "_noThermalModel_HighRContact";

  auto su2 = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 10);
  degrade(*su2);
}

void degradation_variations()
{
  /*
   * Simulate the additional effect of the cell-to-cell variations
   *
   * R only
   * cap only
   * degrate only
   * R and cap and degrate
   */

  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeHistogramData)
    std::cerr << "Warning in paperCode::degradation_variations, we want to do "
                 "a degradation simulation but the global variable CYCLER_STORE_CELL "
                 "is not 1 so no data or too much data will be stored.\n";

  if constexpr (settings::T_MODEL > 0)
    std::cerr << "Warning in paperCode::degradation_variations. We want to simulate "
                 "cells without thermal model. Set settings::T_MODEL to 0 to "
                 "enable this simulation.\n";

  //!< resistance spread
  int coolControl = 1; //!< just for naming since there is no thermal model
  bool capspread = false;
  bool Rcellspread = true;
  bool degratespread = false;
  bool contactR = true;
  std::string IDaddition = "_noThermalModel";
  auto su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1.0);
  degrade(*su);

  //!< capacity spread
  coolControl = 1; //!< just for naming since there is no thermal model
  capspread = true;
  Rcellspread = false;
  degratespread = false;
  contactR = true;
  IDaddition = "_noThermalModel";
  auto su2 = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);
  degrade(*su2);

  //!< degradation rate spread
  coolControl = 1; //!< just for naming since there is no thermal model
  capspread = false;
  Rcellspread = false;
  degratespread = true;
  contactR = true;
  IDaddition = "_noThermalModel";
  auto su3 = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);
  degrade(*su3);

  //!< all three spread
  coolControl = 1; //!< just for naming since there is no thermal model
  capspread = true;
  Rcellspread = true;
  degratespread = true;
  contactR = true;
  IDaddition = "_noThermalModel";
  auto su4 = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);
  degrade(*su4);
}

void degradation_thermal()
{
  /*
   * Simulate the effect of the thermal model
   *
   * no thermal model = before
   * individual cell model = here
   * coupled thermal model = coolsystem simulation (cool1)
   */

  //!< ensure the global settings are appropriately
  if constexpr (settings::DATASTORE_CELL != settings::cellDataStorageLevel::storeHistogramData)
    std::cerr << "Warning in paperCode::degradation_thermal, we want to do a "
                 "degradation simulation but the global variable CYCLER_STORE_CELL "
                 "is not 1 so no data or too much data will be stored.\n";

  if constexpr (settings::T_MODEL != 1)
    std::cerr << "Warning in paperCode::degradation_thermal. We want to simulate "
                 "cells with a noncoupled thermal model. Set settings::T_MODEL "
                 "to 1 to enable this simulation.\n";

  //!< all three spread
  int coolControl = 1; //!< just for naming since there is no thermal model
  bool capspread = true;
  bool Rcellspread = true;
  bool degratespread = true;
  bool contactR = true;
  std::string IDaddition = "_individualThermalModel";
  auto su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);
  degrade(*su);
}
} // namespace slide::paperCode::paper2022
