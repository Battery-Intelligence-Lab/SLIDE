/*
 * Battery.cpp
 *
 *  Created on: 2 Mar 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "makeBattery.hpp"

#include "../cells/cells.hpp"
#include "../modules/Module_s.hpp"
#include "../modules/Module_p.hpp"
#include "../procedures/Cycler.hpp"
#include "../procedures/Procedure.hpp"
// #include "unit_tests.hpp"
#include "../settings/settings.hpp"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ctime>
#include <typeinfo>

namespace slide {
Deep_ptr<StorageUnit> makeBattery(bool balance, bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl)
{
  /*
   * Function to make a large battery pack.
   * For now, make the system from EPFL
   * 		20s cells in module, 15s modules in string, 9p strings in battery
   * 			ignored the 3p in the bottom module (which is 20s3p
   *
   * Names are hierarchical
   */

  //!< String prefix for the names of the results to indicate the settings of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread) settings_str += "_capSpread";
  if (RcellSpread) settings_str += "_RSpread";
  if (degrateSpread) settings_str += "_degSpread";
  if (contactR) settings_str += "_contactR";
  if (balance) settings_str += "_balance";


  //!< Cell-to-cell variations, numbers based on paper from Trish, Jorge (results of screening over 200 prestine lithium-ion cells
  std::default_random_engine generator;
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf = 1;
  double capf = 1;
  double degf = 1;
  double degflam = 1;
  //!< contact resistances
  double Rc_p = 0; //!< contact R for parallel connection
  double Rc_s = 0; //!< contact R for series connection
  if (contactR) {
    Rc_p = 0.001 / 5.0; //!< use a 1mOhm resistance. Value from paper Taheri et al, 2011,
                        //!< Investigating electrical contact resistance losses in lithium-ion battery assemblies for hybrid and electric vehicles
                        //!< 	j Power sources, 2011
                        //!< pouch cell of 5*the capacity -> 5* larger currents -> divide R by 5 or the resistive Vdrop becomes huge
                        //!< this can be a large resistance since they are in parallel, so Rtot reduces with each branch
    Rc_s = 1e-4;        //!< this must be a much smaller resistance
                        //!< 15s20s means we have 15*20 contactR in series. so if Rc_s = 0.01, then the total R is 3 ohm
                        //!< so the resistive voltage drop means cells immediately go outside their voltage limit
                        //!< todo: unrealistically large values, see EPFL battery
                        //!< Paper schimpe has values of 0.1 mOhm = 10^-4
                        //!< Rumpf has 9 * 10^-4
  }

  //!< Degradation settings
  DEG_ID deg;
  /*//!< SEI ONLY

  deg.SEI_id.add_model(4); //!< chirstensen SEI growth
  deg.SEI_porosity = 0; //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0; //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  deg.LAM_id.add_model(0); //!< no LAM
  deg.pl_id = 0; //!< no litihium plating

  /*	//!< Generalised degradation model
  deg.SEI_id.add_model(1); //!< kinetic SEI
  deg.SEI_id.add_model(2); //!< diffusion SEI
  deg.SEI_porosity = 0; //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0; //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  deg.LAM_id.add_model(2); //!< Delacourt LAM
  deg.LAM_id.add_model(3); //!< Kindermann LAM
  deg.pl_id =  1; //!< kinetic litihium plating
  */

  //!< Kinetic SEI + porosity + Dai LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(0);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t nc = 20; //!< number of cells in modules
  constexpr size_t ns = 15; //!< number of modules in strings
  constexpr size_t np = 9;  //!< number of strings in battery
  constexpr bool checkCells = true;

  std::cout << "Total number of cells: " << nc * ns * np << '\n';
  Deep_ptr<StorageUnit> SinP[np]; //!< array with strings in one battery
  Deep_ptr<StorageUnit> MinS[ns]; //!< array with modules in one string
  Deep_ptr<StorageUnit> CinM[nc]; //!< array with cells in one module

  double Rc1[np], Rc2[ns], Rc3[nc];        //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t ip = 0; ip < np; ip++) {     //!< loop for strings in battery
    for (size_t is = 0; is < ns; is++) {   //!< loop for modules in string
      for (size_t ic = 0; ic < nc; ic++) { //!< loop for cells in modules
        if (capSpread)
          capf = distribution_c(generator);
        if (RcellSpread)
          Rf = distribution_r(generator);
        if (degrateSpread) {
          degf = distribution_d(generator);
          degflam = distribution_d(generator);
        }
        CinM[ic] = make_Deep_ptr<Cell_SPM>("cell" + std::to_string(ic), deg, capf, Rf, degf, degflam);
        Rc3[ic] = Rc_s; //!< in series module, so every cell has a resistance of Rc
      }
      auto mi = make_Deep_ptr<Module_s>("s" + std::to_string(is), T_ENV, true, false, nc, coolControl, 0); //!< print warning messages, single-threaded

      mi->setSUs(CinM, checkCells, true);
      mi->setRcontact(Rc3);
      MinS[is] = std::move(mi);
      Rc2[is] = Rc_s; //!< in series module, so every cell has a resistance of Rc
    }
    auto mj = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, nc * ns, coolControl, 2);
    mj->setSUs(MinS, checkCells, true);
    mj->setRcontact(Rc2);
    SinP[ip] = std::move(mj);
    Rc1[ip] = Rc_p; //!< in parallel module, every horizontal branch has the same value of Rc
  }

  auto mp = make<Module_p>("p", T_ENV, true, true, nc * ns * np, coolControl, 2); //!< multithreaded parallel module
  mp->setSUs(SinP, checkCells, true);
  mp->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_pss");
  bat->setModule(std::move(mp));

  std::cout << "The voltage of the battery is " << bat->V() << " and it has a capacity of " << bat->Cap() << '\n';
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*Ah = 115.699Ah
  return bat;

  /*
  //!< Test some cycle ageing
  bool CV = false;
  double Vbal = 3.5;							//!< rebalance to 3.1V
  int ndata = 100;							//!< how often to store data (if 0, nothing is stored, not even usage stats)
  Deep_ptr<Procedure> p(new Procedure(balance, Vbal, ndata, false));
  //p->useCaseAge(bat, coolControl);
  p->cycleAge(bat, CV);*/
}

Deep_ptr<StorageUnit> makeBattery2(bool balance, bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl)
{
  /*
   * Function to make a large battery pack.
   * This battery has the opposite structure from the other one:
   * 		9p 15s 20s
   *
   * Names are hierarchical
   */

  //!< String prefix for the names of the results to indicate the settings_str of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread) settings_str += "_capSpread";
  if (RcellSpread) settings_str += "_RSpread";
  if (degrateSpread) settings_str += "_degSpread";
  if (contactR) settings_str += "_contactR";
  if (balance) settings_str += "_balance";

  //!< Cell-to-cell variations

  //!< unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 2;

  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf{ 1 }, capf{ 1 }, degf{ 1 }, degflam{ 1 };
  //!< contact resistances
  double Rc_p = 0; //!< contact R for parallel connection
  double Rc_s = 0; //!< contact R for series connection
  if (contactR) {
    Rc_p = 0.001 / 5.0; //!< use a 1mOhm resistance. Value from paper Taheri et al, 2011,
                        //!< Investigating electrical contact resistance losses in lithium-ion battery assemblies for hybrid and electric vehicles
                        //!< 	j Power sources, 2011
                        //!< pouch cell of 5*the capacity -> 5* larger currents -> divide R by 5 or the resistive Vdrop becomes huge
    Rc_s = 1e-4;        //!< this must be a much smaller resistance
                        //!< 15s20s means we have 15*20 contactR in series. so if Rc_s = 0.01, then the total R is 3 ohm
                        //!< so the resistive voltage drop means cells immediately go outside their voltage limit
                        //!< todo: unrealistically large values see make EFFL battery
  }

  //!< Degradation settings
  DEG_ID deg;
  /*//!< SEI ONLY
  deg.SEI_id.N = 1;										//!< there is 1 SEI model
  deg.SEI_id[0] = 4;									//!< chirstensen SEI growth
  deg.SEI_porosity = 0;								//!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)
  deg.CS_id.N = 1;										//!< there is 1 model (that there are no cracks)
  deg.CS_id[0] = 0;									//!< no surface cracks
  deg.CS_diffusion = 0;								//!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)
  deg.LAM_id.N = 1;										//!< there are 1 LAM model
  deg.LAM_id[0] = 0;									//!< no LAM
  deg.pl_id = 0;										//!< no litihium plating*/

  /*	//!< Generalised degradation model
  deg.SEI_id.N = 2;										//!< there is 1 SEI model
  deg.SEI_id[0] = 1;									//!< kinetic SEI
  deg.SEI_id[1] = 2;									//!< diffusion SEI
  deg.SEI_porosity = 0;								//!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)
  deg.CS_id.N = 1;										//!< there is 1 model (that there are no cracks)
  deg.CS_id[0] = 0;									//!< no surface cracks
  deg.CS_diffusion = 0;								//!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)
  deg.LAM_id.N = 2;										//!< there are 1 LAM model
  deg.LAM_id[0] = 2;									//!< Delacourt LAM
  deg.LAM_id[1] = 3;									//!< Kindermann LAM
  deg.pl_id = 1;										//!< kinetic litihium plating*/
  //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t nc = 9;  //!< number of cells in (parallel) modules
  constexpr size_t ns = 15; //!< number of modules in (s) strings
  constexpr size_t np = 20; //!< number of strings in (s) battery
  constexpr bool checkCells = true;
  std::cout << "Total number of cells: " << nc * ns * np << '\n';

  Deep_ptr<StorageUnit> SinP[np];   //!< array with strings in one battery
  Deep_ptr<StorageUnit> MinS[ns];   //!< array with modules in one string
  Deep_ptr<StorageUnit> CinM[nc];   //!< array with cells in one module
  double Rc1[np], Rc2[ns], Rc3[nc]; //!< arrays with the contact resistances for the different levels (3 = lowest level)

  for (size_t ip = 0; ip < np; ip++) //!< loop for strings in battery
  {
    for (size_t is = 0; is < ns; is++) //!< loop for modules in string
    {
      for (size_t ic = 0; ic < nc; ic++) { //!< loop for cells in modules
        if (capSpread)
          capf = distribution_c(generator);
        if (RcellSpread)
          Rf = distribution_r(generator);
        if (degrateSpread) {
          degf = distribution_d(generator);
          degflam = distribution_d(generator);
        }
        CinM[ic] = make<Cell_SPM>("cell" + std::to_string(ic), deg, capf, Rf, degf, degflam);
        Rc3[ic] = Rc_p; //!< in parallel module, all branches have same R
      }
      auto mi = make<Module_p>("p" + std::to_string(is), T_ENV, true, false, nc, coolControl, 0); //!< print warning messages, single-threaded
      mi->setSUs(CinM, checkCells, true);
      mi->setRcontact(Rc3);

      MinS[is] = std::move(mi);
      Rc2[is] = Rc_s; //!< in series module, so every cell has a resistance of Rc
    }
    auto mj = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, nc * ns, coolControl, 2); //!< single-threaded
    mj->setSUs(MinS, checkCells, true);
    mj->setRcontact(Rc2);

    SinP[ip] = std::move(mj);
    Rc1[ip] = Rc_s; //!< in parallel module, every horizontal branch has the same value of Rc
  }
  auto ms = make<Module_s>("s", T_ENV, true, true, nc * ns * np, coolControl, 2); //!< multithreaded series module
  ms->setSUs(SinP, checkCells, true);
  ms->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_ssp");
  bat->setModule(std::move(ms));

  std::cout << "The voltage of the battery is " << bat->V() << " and it has a capacity of " << bat->Cap() << '\n';
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*3Ah = 27Ah

  return bat;
}

Deep_ptr<StorageUnit> makeBattery_EPFL_smaller(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM)
{
  /*
   * Function to make a large battery pack.
   * Make the structure of the EPFL battery
   * 	modules are 20s 6p [note: they have 3p of larger cells]
   * 	15 of those modules in s to form a rack
   * 	9 of these racks in p to form battery
   *
   * 	RM 	resistance multiplier (Rcontact = standard * RM)
   */

  //!< String prefix for the names of the results to indicate the settings of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread)
    settings_str += "_capSpread";
  if (RcellSpread)
    settings_str += "_RSpread";
  if (degrateSpread)
    settings_str += "_degSpread";
  if (contactR)
    settings_str += "_contactR" + std::to_string(RM);
  settings_str += IDadditions;

  //!< Cell-to-cell variations

  //!< unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 2;

  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf = 1;
  double capf = 1;
  double degf = 1;
  double degflam = 1;
  //!< contact resistances
  double Rc_p_cells = 0; //!< contact R for parallel connection or cells
  double Rc_s_cells = 0; //!< contact R for series connections of cells
  double Rc_p = 0;       //!< contact R for parallel connection of racks
  double Rc_s = 0;       //!< contact R for series connection of modules
  if (contactR) {
    Rc_p_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_s_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_p = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
    Rc_s = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
  }

  //!< Degradation settings
  DEG_ID deg;
  //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 0;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t ncp = 7;      //!< number of cells in parallel in a module [lowest level]
  constexpr size_t ncs = 3;      //!< number of cells in series in a module
  constexpr size_t nmodules = 3; //!< number of modules in (s) racks
  constexpr size_t nracks = 2;   //!< number of strings in (p) battery

  constexpr bool checkCells = true;
  Deep_ptr<StorageUnit> RinB[nracks];   //!< array with racks in one battery
  Deep_ptr<StorageUnit> MinR[nmodules]; //!< array with modules in one rack
  Deep_ptr<StorageUnit> CinsM[ncs];     //!< array with cells in series in one module
  Deep_ptr<StorageUnit> CinpM[ncp];     //!< array with cells in parallel in one module

  double Rc1[nracks], Rc2[nmodules], Rc3[ncs], Rc4[ncp]; //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t ip = 0; ip < nracks; ip++) {               //!< loop for strings in battery
    for (size_t is = 0; is < nmodules; is++) {           //!< loop for modules in string
      for (size_t ics = 0; ics < ncs; ics++) {           //!< loop for cells in series in modules
        for (size_t icp = 0; icp < ncp; icp++) {         //!< loop for cells in parallel in modules
          if (capSpread)
            capf = distribution_c(generator);
          if (RcellSpread)
            Rf = distribution_r(generator);
          if (degrateSpread) {
            degf = distribution_d(generator);
            degflam = distribution_d(generator);
          }
          CinpM[icp] = make<Cell_SPM>("cell" + std::to_string(icp), deg, capf, Rf, degf, degflam);
          Rc4[icp] = Rc_p_cells; //!< in parallel module, all branches have same R
        }                        //!< loop to make cells connected in parallel in modules

        //!< put the cells in parallel to the lowest-level module
        auto mip = make<Module_p>("p" + std::to_string(ics), T_ENV, true, false, ncp, coolControl, 2); //!< print warning messages, single-threaded, pass through coolsystem

        mip->setSUs(CinpM, checkCells, true);
        mip->setRcontact(Rc4);
        CinsM[ics] = std::move(mip);
        Rc3[ics] = Rc_s_cells; //!< in series module, so every cell has a resistance of Rc
      }                        //!< loop to make the parallel-connected cells which goes in series to make one module

      //!< assemble the cells in series for the
      auto modulei = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ncp * ncs, coolControl, 0); //!< single-threaded, conventional coolsystem
      modulei->setSUs(CinsM, checkCells, true);
      modulei->setRcontact(Rc3);
      MinR[is] = std::move(modulei);
      Rc2[is] = Rc_s;


    } //!< loop to make the modules for one rack

    //!< assemble the modules in series for a rack
    auto racki = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, ncp * ncs * nmodules, coolControl, 2); //!< single-threaded, pass through coolsystem
    racki->setSUs(MinR, checkCells, true);
    racki->setRcontact(Rc2);
    RinB[ip] = std::move(racki);
    Rc1[ip] = Rc_p;

  } //!< loop to make the racks

  //!< Assemble the racks in the battery compartment (bc)
  auto bc = make<Module_p>("p", T_ENV, true, true, ncp * ncs * nmodules * nracks, coolControl, 2); //!< multithreaded parallel module, pass through coolsystem
  bc->setSUs(RinB, checkCells, true);
  bc->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_EPFL"); //!< battery, gets HVAC coolsystem
  bat->setModule(std::move(bc));

  std::cout << "Total number of cells: " << bat->getNcells() << '\n';

  const auto bat_v = bat->V();
  const auto bat_cap = bat->Cap();
  std::cout << "The voltage of the battery is " << bat_v << " and it has a capacity of " << bat_cap << std::endl;
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*3Ah = 27Ah
  return bat;
}

Deep_ptr<StorageUnit> makeBattery_EPFL(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM)
{
  /*
   * Function to make a large battery pack.
   * Make the structure of the EPFL battery
   * 	modules are 20s 6p [note: they have 3p of larger cells]
   * 	15 of those modules in s to form a rack
   * 	9 of these racks in p to form battery
   *
   * 	RM 	resistance multiplier (Rcontact = standard * RM)
   */

  //!< String prefix for the names of the results to indicate the settings of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread)
    settings_str += "_capSpread";
  if (RcellSpread)
    settings_str += "_RSpread";
  if (degrateSpread)
    settings_str += "_degSpread";
  if (contactR)
    settings_str += "_contactR" + std::to_string(RM);
  settings_str += IDadditions;

  //!< Cell-to-cell variations

  //!< unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 2;

  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf = 1;
  double capf = 1;
  double degf = 1;
  double degflam = 1;
  //!< contact resistances
  double Rc_p_cells = 0; //!< contact R for parallel connection or cells
  double Rc_s_cells = 0; //!< contact R for series connections of cells
  double Rc_p = 0;       //!< contact R for parallel connection of racks
  double Rc_s = 0;       //!< contact R for series connection of modules
  if (contactR) {
    Rc_p_cells = 0.0075e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_s_cells = 0.0075e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_p = 0.25e-3 * RM;         //!< this must be a much smaller resistance
    Rc_s = 0.25e-3 * RM;         //!< this must be a much smaller resistance
  }

  //!< Degradation settings
  DEG_ID deg;
  //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t ncp = 7;       //!< number of cells in parallel in a module [lowest level]
  constexpr size_t ncs = 20;      //!< number of cells in series in a module
  constexpr size_t nmodules = 15; //!< number of modules in (s) racks
  constexpr size_t nracks = 9;    //!< number of strings in (p) battery

  constexpr bool checkCells = true;
  Deep_ptr<StorageUnit> RinB[nracks];   //!< array with racks in one battery
  Deep_ptr<StorageUnit> MinR[nmodules]; //!< array with modules in one rack
  Deep_ptr<StorageUnit> CinsM[ncs];     //!< array with cells in series in one module
  Deep_ptr<StorageUnit> CinpM[ncp];     //!< array with cells in parallel in one module

  double Rc1[nracks], Rc2[nmodules], Rc3[ncs], Rc4[ncp]; //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t ip = 0; ip < nracks; ip++) {               //!< loop for strings in battery
    for (size_t is = 0; is < nmodules; is++) {           //!< loop for modules in string
      for (size_t ics = 0; ics < ncs; ics++) {           //!< loop for cells in series in modules
        for (size_t icp = 0; icp < ncp; icp++) {         //!< loop for cells in parallel in modules
          if (capSpread)
            capf = distribution_c(generator);
          if (RcellSpread)
            Rf = distribution_r(generator);
          if (degrateSpread) {
            degf = distribution_d(generator);
            degflam = distribution_d(generator);
          }
          CinpM[icp] = make<Cell_SPM>("cell" + std::to_string(icp), deg, capf, Rf, degf, degflam);
          Rc4[icp] = Rc_p_cells; //!< in parallel module, all branches have same R
        }                        //!< loop to make cells connected in parallel in modules

        //!< put the cells in parallel to the lowest-level module
        auto mip = make<Module_p>("p" + std::to_string(ics), T_ENV, true, false, ncp, coolControl, 2); //!< print warning messages, single-threaded, pass through coolsystem

        mip->setSUs(CinpM, checkCells, true);
        mip->setRcontact(Rc4);
        CinsM[ics] = std::move(mip);
        Rc3[ics] = Rc_s_cells; //!< in series module, so every cell has a resistance of Rc
      }                        //!< loop to make the parallel-connected cells which goes in series to make one module

      //!< assemble the cells in series for the
      auto modulei = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ncp * ncs, coolControl, 0); //!< single-threaded, conventional coolsystem
      modulei->setSUs(CinsM, checkCells, true);
      modulei->setRcontact(Rc3);
      MinR[is] = std::move(modulei);
      Rc2[is] = Rc_s;

    } //!< loop to make the modules for one rack

    //!< assemble the modules in series for a rack
    auto racki = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, ncp * ncs * nmodules, coolControl, 2); //!< single-threaded, pass through coolsystem
    racki->setSUs(MinR, checkCells, true);
    racki->setRcontact(Rc2);
    RinB[ip] = std::move(racki);
    Rc1[ip] = Rc_p;

  } //!< loop to make the racks

  //!< Assemble the racks in the battery compartment (bc)
  auto bc = make<Module_p>("p", T_ENV, true, true, ncp * ncs * nmodules * nracks, coolControl, 2); //!< multithreaded parallel module, pass through coolsystem
  bc->setSUs(RinB, checkCells, true);
  bc->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_EPFL"); //!< battery, gets HVAC coolsystem
  bat->setModule(std::move(bc));

  std::cout << "Total number of cells: " << bat->getNcells() << '\n';
  std::cout << "The voltage of the battery is " << bat->V() << " and it has a capacity of " << bat->Cap() << std::endl;
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*3Ah = 27Ah
  return bat;
}

Deep_ptr<StorageUnit> makeBattery_Test(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM)
{
  /*
   * Function to make a large battery pack.
   * Make the structure of the EPFL battery
   * 	modules are 20s 6p [note: they have 3p of larger cells]
   * 	15 of those modules in s to form a rack
   * 	9 of these racks in p to form battery
   *
   * 	RM 	resistance multiplier (Rcontact = standard * RM)
   */

  //!< String prefix for the names of the results to indicate the settings of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread)
    settings_str += "_capSpread";
  if (RcellSpread)
    settings_str += "_RSpread";
  if (degrateSpread)
    settings_str += "_degSpread";
  if (contactR)
    settings_str += "_contactR" + std::to_string(RM);
  settings_str += IDadditions;

  //!< Cell-to-cell variations

  //!< unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 2;

  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf = 1;
  double capf = 1;
  double degf = 1;
  double degflam = 1;
  //!< contact resistances
  double Rc_p_cells = 0; //!< contact R for parallel connection or cells
  double Rc_s_cells = 0; //!< contact R for series connections of cells
  double Rc_p = 0;       //!< contact R for parallel connection of racks
  double Rc_s = 0;       //!< contact R for series connection of modules
  if (contactR) {
    Rc_p_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_s_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_p = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
    Rc_s = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
  }

  //!< Degradation settings
  DEG_ID deg;
  //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t ncp = 3;      //!< number of cells in parallel in a module [lowest level]
  constexpr size_t ncs = 4;      //!< number of cells in series in a module
  constexpr size_t nmodules = 2; //!< number of modules in (s) racks
  constexpr size_t nracks = 3;   //!< number of strings in (p) battery
  std::cout << "Total number of cells: " << ncp * ncs * nmodules * nracks << '\n';

  constexpr bool checkCells = true;
  Deep_ptr<StorageUnit> RinB[nracks];   //!< array with racks in one battery
  Deep_ptr<StorageUnit> MinR[nmodules]; //!< array with modules in one rack
  Deep_ptr<StorageUnit> CinsM[ncs];     //!< array with cells in series in one module
  Deep_ptr<StorageUnit> CinpM[ncp];     //!< array with cells in parallel in one module

  double Rc1[nracks], Rc2[nmodules], Rc3[ncs], Rc4[ncp]; //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t ip = 0; ip < nracks; ip++) {               //!< loop for strings in battery
    for (size_t is = 0; is < nmodules; is++) {           //!< loop for modules in string
      for (size_t ics = 0; ics < ncs; ics++) {           //!< loop for cells in series in modules
        for (size_t icp = 0; icp < ncp; icp++) {         //!< loop for cells in parallel in modules
          if (capSpread)
            capf = distribution_c(generator);
          if (RcellSpread)
            Rf = distribution_r(generator);
          if (degrateSpread) {
            degf = distribution_d(generator);
            degflam = distribution_d(generator);
          }
          CinpM[icp] = make<Cell_SPM>("cell" + std::to_string(icp), deg, capf, Rf, degf, degflam);
          Rc4[icp] = Rc_p_cells; //!< in parallel module, all branches have same R
        }                        //!< loop to make cells connected in parallel in modules

        //!< put the cells in parallel to the lowest-level module
        auto mip = make<Module_p>("p" + std::to_string(ics), T_ENV, true, false, ncp, coolControl, 2); //!< print warning messages, single-threaded, pass through coolsystem

        mip->setSUs(CinpM, checkCells, true);
        mip->setRcontact(Rc4);
        CinsM[ics] = std::move(mip);
        Rc3[ics] = Rc_s_cells; //!< in series module, so every cell has a resistance of Rc
      }                        //!< loop to make the parallel-connected cells which goes in series to make one module

      //!< assemble the cells in series for the
      auto modulei = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ncp * ncs, coolControl, 0); //!< single-threaded, conventional coolsystem
      modulei->setSUs(CinsM, checkCells, true);
      modulei->setRcontact(Rc3);
      MinR[is] = std::move(modulei);
      Rc2[is] = Rc_s;

    } //!< loop to make the modules for one rack

    //!< assemble the modules in series for a rack
    auto racki = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, ncp * ncs * nmodules, coolControl, 2); //!< single-threaded, pass through coolsystem
    racki->setSUs(MinR, checkCells, true);
    racki->setRcontact(Rc2);
    RinB[ip] = std::move(racki);
    Rc1[ip] = Rc_p;

  } //!< loop to make the racks

  //!< Assemble the racks in the battery compartment (bc)
  auto bc = make<Module_p>("p", T_ENV, true, true, ncp * ncs * nmodules * nracks, coolControl, 2); //!< multithreaded parallel module, pass through coolsystem
  bc->setSUs(RinB, checkCells, true);
  bc->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_EPFL"); //!< battery, gets HVAC coolsystem
  bat->setModule(std::move(bc));

  std::cout << "The voltage of the battery is " << bat->V() << " and it has a capacity of " << bat->Cap() << '\n';
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*3Ah = 27Ah
  return bat;
}

Deep_ptr<StorageUnit> makeBattery_TestParallel(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM)
{
  /*
   * Function to make a large battery pack.
   * Make the structure of the EPFL battery
   * 	modules are 20s 6p [note: they have 3p of larger cells]
   * 	15 of those modules in s to form a rack
   * 	9 of these racks in p to form battery
   *
   * 	RM 	resistance multiplier (Rcontact = standard * RM)
   */

  //!< String prefix for the names of the results to indicate the settings of the simulation
  using settings::T_ENV;
  std::string settings_str = "cool" + std::to_string(coolControl);
  if (capSpread)
    settings_str += "_capSpread";
  if (RcellSpread)
    settings_str += "_RSpread";
  if (degrateSpread)
    settings_str += "_degSpread";
  if (contactR)
    settings_str += "_contactR" + std::to_string(RM);
  settings_str += IDadditions;

  //!< Cell-to-cell variations

  //!< unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 2;

  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution_c(1.0, 0.004); //!< normal distribution with mean 1 and std 0.4%
  std::normal_distribution<double> distribution_r(1.0, 0.025); //!< normal distribution with mean 1 and std 2.5%
  std::normal_distribution<double> distribution_d(1.0, 0.10);  //!< normal distribution with mean 1 and std 10%
  //!< factors with the relative change in resistance and capacity of the cells
  double Rf = 1;
  double capf = 1;
  double degf = 1;
  double degflam = 1;
  //!< contact resistances
  double Rc_p_cells = 0; //!< contact R for parallel connection or cells
  double Rc_s_cells = 0; //!< contact R for series connections of cells
  double Rc_p = 0;       //!< contact R for parallel connection of racks
  double Rc_s = 0;       //!< contact R for series connection of modules
  if (contactR) {
    Rc_p_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_s_cells = 0.0075 * 1e-3 * RM; //!< use a 0.1 mOhm resistance. Value from paper Schimpe
    Rc_p = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
    Rc_s = 0.25 * 1e-3 * RM;         //!< this must be a much smaller resistance
  }

  //!< Degradation settings
  DEG_ID deg;
  //!< Kinetic SEI + porosity + Dai/Laresgoiti LAM
  deg.SEI_id.add_model(4);
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0);
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  //!< Make the battery
  constexpr size_t ncp = 3;      //!< number of cells in parallel in a module [lowest level]
  constexpr size_t ncs = 4;      //!< number of cells in series in a module
  constexpr size_t nmodules = 2; //!< number of modules in (s) racks
  constexpr size_t nracks = 3;   //!< number of strings in (p) battery
  std::cout << "Total number of cells: " << ncp * ncs * nmodules * nracks << '\n';

  constexpr bool checkCells = true;
  Deep_ptr<StorageUnit> RinB[nracks];   //!< array with racks in one battery
  Deep_ptr<StorageUnit> MinR[nmodules]; //!< array with modules in one rack
  Deep_ptr<StorageUnit> CinsM[ncs];     //!< array with cells in series in one module
  Deep_ptr<StorageUnit> CinpM[ncp];     //!< array with cells in parallel in one module

  double Rc1[nracks], Rc2[nmodules], Rc3[ncs], Rc4[ncp]; //!< arrays with the contact resistances for the different levels (3 = lowest level)

  for (size_t icp = 0; icp < ncp; icp++) { //!< loop for cells in parallel in modules
    if (capSpread)
      capf = distribution_c(generator);
    if (RcellSpread)
      Rf = distribution_r(generator);
    if (degrateSpread) {
      degf = distribution_d(generator);
      degflam = distribution_d(generator);
    }
    CinpM[icp] = make<Cell_SPM>("cell" + std::to_string(icp), deg, capf, Rf, degf, degflam);
    Rc4[icp] = Rc_p_cells; //!< in parallel module, all branches have same R
  }                        //!< loop to make cells connected in parallel in modules

  //!< put the cells in parallel to the lowest-level module
  auto mip = make<Module_p>("p" + std::to_string(0), T_ENV, true, false, ncp, coolControl, 2); //!< print warning messages, single-threaded, pass through coolsystem

  mip->setSUs(CinpM, checkCells, true);
  mip->setRcontact(Rc4);

  //!< for (size_t ip = 0; ip < nracks; ip++)
  //!< { //!< loop for strings in battery
  //!< 	for (size_t is = 0; is < nmodules; is++)
  //!< 	{ //!< loop for modules in string
  //!< 		for (size_t ics = 0; ics < ncs; ics++)
  //!< 		{ //!< loop for cells in series in modules
  //!< 			for (size_t icp = 0; icp < ncp; icp++)
  //!< 			{ //!< loop for cells in parallel in modules
  //!< 				if (capSpread)
  //!< 					capf = distribution_c(generator);
  //!< 				if (RcellSpread)
  //!< 					Rf = distribution_r(generator);
  //!< 				if (degrateSpread)
  //!< 				{
  //!< 					degf = distribution_d(generator);
  //!< 					degflam = distribution_d(generator);
  //!< 				}
  //!< 				CinpM[icp] = make<Cell_SPM>("cell" + std::to_string(icp), deg, capf, Rf, degf, degflam);
  //!< 				Rc4[icp] = Rc_p_cells; //!< in parallel module, all branches have same R
  //!< 			}						   //!< loop to make cells connected in parallel in modules

  //!< 			//!< put the cells in parallel to the lowest-level module
  //!< 			auto mip = make<Module_p>("p" + std::to_string(ics), T_ENV, true, false, ncp, coolControl, 2); //!< print warning messages, single-threaded, pass through coolsystem

  //!< 			mip->setSUs(CinpM, checkCells, true);
  //!< 			mip->setRcontact(Rc4);
  //!< 			CinsM[ics] = std::move(mip);
  //!< 			Rc3[ics] = Rc_s_cells; //!< in series module, so every cell has a resistance of Rc
  //!< 		}						   //!< loop to make the parallel-connected cells which goes in series to make one module

  //!< 		//!< assemble the cells in series for the
  //!< 		auto modulei = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ncp * ncs, coolControl, 0); //!< single-threaded, conventional coolsystem
  //!< 		modulei->setSUs(CinsM, checkCells, true);
  //!< 		modulei->setRcontact(Rc3);
  //!< 		MinR[is] = std::move(modulei);
  //!< 		Rc2[is] = Rc_s;

  //!< 	} //!< loop to make the modules for one rack

  //!< 	//!< assemble the modules in series for a rack
  //!< 	auto racki = make<Module_s>("s" + std::to_string(ip), T_ENV, true, false, ncp * ncs * nmodules, coolControl, 2); //!< single-threaded, pass through coolsystem
  //!< 	racki->setSUs(MinR, checkCells, true);
  //!< 	racki->setRcontact(Rc2);
  //!< 	RinB[ip] = std::move(racki);
  //!< 	Rc1[ip] = Rc_p;

  //!< } //!< loop to make the racks

  //!< //!< Assemble the racks in the battery compartment (bc)
  //!< auto bc = make<Module_p>("p", T_ENV, true, true, ncp * ncs * nmodules * nracks, coolControl, 2); //!< multithreaded parallel module, pass through coolsystem
  //!< bc->setSUs(RinB, checkCells, true);
  //!< bc->setRcontact(Rc1);

  //!< make the battery
  auto bat = make<Battery>(settings_str + "_ParTest"); //!< battery, gets HVAC coolsystem
  bat->setModule(std::move(mip));

  std::cout << "The voltage of the battery is " << bat->V() << " and it has a capacity of " << bat->Cap() << '\n';
  //!< 20 * 15 cells in series -> 20*15*3.68V = 1,104V
  //!< 9 strings in parallel -> 9*3Ah = 27Ah
  return bat;
}
} // namespace slide