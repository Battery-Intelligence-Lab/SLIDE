/*
 * Procedure.cpp
 *
 *  Created on: 3 Mar 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../cells/cells.hpp"
#include "../modules/modules.hpp"
#include "../system/Battery.hpp"
// #include "unit_tests.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"
#include "Cycler.hpp"
#include "Procedure.hpp"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ctime>
#include <typeinfo>
#include <vector>
#include <memory>
#include <mutex>

namespace slide {
//!< Global variables to be used
std::mutex MUTEX_procedure_checkUp; //!< the mutex to lock the file when writing checkup results

Procedure::Procedure(bool balance_, double Vbal_, int ndata_, bool unitTest_)
  : balance{ balance_ }, unitTest{ unitTest_ }, ndata{ ndata_ }, balance_voltage{ Vbal_ }
{
}

void Procedure::cycleAge(StorageUnit *su, int Ncycle, int Ncheck, int Nbal, bool testCV, double Ccha, double Cdis, double Vmax, double Vmin)
{
  /*
   * Simulate a cycle ageing experiment with the given parameters
   */

  //!< Initialise the Cycler
  std::string pref = "cycleAge"; //!< prefix appended at the start of the names of all files
  const bool diagnostic = true;  //!< stop (dis)charging when one of the cells reaches a voltage limit
  Cycler cyc;
  cyc.initialise(su, pref);
  cyc.setDiagnostic(diagnostic);

  //!< variables
  const double dt = 2;
  const double Icha = -Ccha * su->Cap();
  const double Idis = Cdis * su->Cap();
  double v{};
  Status succ{};
  double tlim, Ilim, Ah, Wh, dtime{};

  //!< ID for storing the steps in the producedure (i.e. which data refers to which action)
  constexpr int ID_CCcharge = 1;
  constexpr int ID_CVcharge = 2;
  constexpr int ID_CCdischarge = 3;
  constexpr int ID_CVdischarge = 4;

  //!< Make a clock to measure how long the simulation takes
  Clock clk;

  //!< loop for cycle ageing
  double Ahtot = 0;
  try {
    for (int i = 0; i < Ncycle; i++) {
      //!< check duration of the simulation
      v = su->V();

      if (!unitTest)
        std::cout << "SU " << su->getFullID() << " starting loop iteration " << i << " after "
                  << clk << " with V = " << v << ", T = " << K_to_Celsius(su->T())
                  << " and hot spot T = " << K_to_Celsius(su->getThotSpot()) << '\n';

      //!< Balance (in this thread) and do a check-up (on a separate thread) if needed
      balanceCheckup(su, balance && (i % Nbal == 0), (i % Ncheck == 0), Ahtot, i, pref);

      //!< CC charge
      tlim = std::numeric_limits<double>::max();
      Ilim = 0.1;


      try {
        dtime = Ah = Wh = 0; //!< reset in case error in CC and Ah gets not changed (without reset it would keep its old value)
        succ = cyc.CC(Icha, Vmax, tlim, dt, ndata, Ah, Wh, dtime);
      } catch (int e) {
        std::cout << "Error in CycleAge when charging in cycle "
                  << i << ", stop cycling.\n";
        break;
      }
      Ahtot += std::abs(Ah);
      storeThroughput(ID_CCcharge, Ah, Wh, su); //!< increase the throughput parameters

      if (!diagnostic)
        assert(succ == Status::ReachedVoltageLimit);

      if (Ah >= 0) {
        std::cerr << "Error in Procedure::cycleAge, we didn't manage to charge "
                  << "any energy. Stopping the cycling after " << i << " cycles.\n";
        break;
      }

      //!< CV charge
      if (testCV) {
        try {
          dtime = Ah = Wh = 0; //!< reset in case error in CC and Ah gets not changed (without reset it would keep its old value)
          succ = cyc.CV(Vmax, Ilim, tlim, dt, ndata, Ah, Wh, dtime);
        } catch (int e) {
          std::cout << "Error in CycleAge when CV charging in cycle "
                    << i << ", stopping cycling.\n";
          break;
        }

        Ahtot += std::abs(Ah);
        storeThroughput(ID_CVcharge, Ah, Wh, su);
        if (!diagnostic) //!< #TODO what is diagnostic?
          assert(succ == Status::ReachedCurrentLimit);
      }

      //!< CC discharge
      try {
        dtime = Ah = Wh = 0; //!< reset in case error in CC and Ah gets not changed (without reset it would keep its old value)
        succ = cyc.CC(Idis, Vmin, tlim, dt, ndata, Ah, Wh, dtime);
      } catch (int e) {
        std::cout << "Error in CycleAge when discharging in cycle " << i << ", stop cycling.\n";
        break;
      }
      Ahtot += std::abs(Ah);
      storeThroughput(ID_CCdischarge, Ah, Wh, su);
      if (!diagnostic)
        assert(succ == Status::ReachedVoltageLimit);

      //!< CV discharge
      if (testCV) {
        try {
          dtime = Ah = Wh = 0; //!< reset in case error in CC and Ah gets not changed (without reset it would keep its old value)
          succ = cyc.CV(Vmin, Ilim, tlim, dt, ndata, Ah, Wh, dtime);
        } catch (int e) //!< #TODO if we need to check any status codes here?
        {
          std::cout << "Error in CycleAge when CV discharging in cycle " << i << ", stop cycling.\n";
          break;
        }
        Ahtot += std::abs(Ah);
        storeThroughput(ID_CVdischarge, Ah, Wh, su);
        if (!diagnostic)
          assert(succ == Status::ReachedCurrentLimit);
      }

    }           //!< loop cycle ageing
  } catch (...) //!< #TODO why empty?
  {
    std::cout << "Why do we enter here?\n";
  }

  //!< final checkup
  //!< This checkup will write the usage statistics, so it must be called with the actual SU and not with a copy
  //!< because usage statistics of cells and modules (or their cooling systems) are not copied over in copy()
  try {
    if (!unitTest)
      std::cout << "Doing a final checkup after " << Ahtot << " Ah.\n";
    writeThroughput(su->getFullID(), Ahtot);

    auto task_indv = [&](int i) {
      if (i == 0)
        checkUp(su, Ahtot, Ncycle); //!< check-up of the cells
      else if (i == 1)
        checkMod(su); //!< check-up of the modules
    };

    run(task_indv, 2, 2); //!< #TODO if parallel in two cores is really needed.
  } catch (int e) {
    std::cout << "Error in CycleAge when getting a capacity checkup, skipping it\n";
  }

  //!< push a write such that if cells still have cycling data, this is written
  if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeTimeData)
    su->writeData(pref); //!< only do if cycling data. Usage statistics are written by the checkup
}

void Procedure::cycleAge(StorageUnit *su, bool testCV)
{
  /*
   * Test some simple cycle ageing for an arbitrarily large battery
   *
   */
  //!< Call the advanced function with standard settings
  const double Ccha{ 1 }, Cdis{ 1 };
  const double Vmax = su->Vmax();
  const double Vmin = su->Vmin();
  const unsigned Ncycle = 15000;
  const unsigned ncheck = 250; //!< do a checkup ever 250 cycles
  const unsigned nbal = 10;    //!< balance every 10 cycles

  cycleAge(su, Ncycle, ncheck, nbal, testCV, Ccha, Cdis, Vmax, Vmin);
}

void Procedure::useCaseAge(StorageUnit *su, int cool)
{
  /*
   * Simulate a realistic use case. The usage pattern is (starting from fully discharged at midnight)
   *
   * rest 4h fully discharged
   * 1C charge (4am to 5am)
   * rest 1h fully charged (5am to 6am)
   * 1C discharge (6 to 7am)
   * rest fully discharged from 7am to 11am
   * 0.5C charge during noon (11am to 1pm)
   * rest rully charged from 1pm to 5pm
   * 0.5C discharge from 5pm to 7pm
   * rest fully discharged until midnight
   */

  //!< settings
  const bool diagnostic = true; //!< stop (dis)charging when one of the cells reaches a voltage limit
  std::string pref = "useAge";  //!< prefix appended at the start of the names of all files
  const double dt = 1;
  constexpr unsigned Ncycle = 20; //!< 365 * 10; //!< 10 year
  constexpr unsigned Ncheck = 10; //!< do a checkup ever month
  constexpr unsigned Nbal = 7;    //!< balance every week

  //!< Variables
  Cycler cyc;
  double Ah, Wh, v;
  Status succ;
  cyc.initialise(su, pref);
  cyc.setDiagnostic(diagnostic);
  constexpr int ID_rest1 = 1;
  constexpr int ID_cha1 = 2;
  constexpr int ID_rest2 = 3;
  constexpr int ID_dis1 = 4;
  constexpr int ID_rest3 = 5;
  constexpr int ID_cha2 = 6;
  constexpr int ID_rest4 = 7;
  constexpr int ID_dis2 = 8;
  constexpr int ID_rest5 = 9;

  //!< fully discharge the SU at a C/2
  double tlim = std::numeric_limits<double>::max();
  try {
    double dtime;
    succ = cyc.CC(su->Cap() / 2.0, su->Vmin(), tlim, dt, ndata, Ah, Wh, dtime);
    //!< std::cout << getStatusMessage(succ) << '\n'; //!< #TODO if we put and test all conditions.
  } catch (int e) {
    std::cout << "Error in useAge when initially discharging the battery, continue as normal.\n";
  }

  //!< Make a clock to measure how long the simulation takes
  Clock clk;

  //!< loop for use case ageing
  double Ahtot = 0;
  for (unsigned i = 0; i < Ncycle; i++) {
    //!< check duration of the simulation
    v = su->V();

    if (!unitTest)
      std::cout << "SU " << su->getFullID() << " starting loop iteration " << i << " after "
                << clk << " with V = " << v << ", T = " << K_to_Celsius(su->T())
                << " and hot spot T = " << K_to_Celsius(su->getThotSpot()) << '\n';

    //!< Balance (in this thread) and do a check-up (on a separate thread) if needed
    balanceCheckup(su, balance && (i % Nbal == 0), (i % Ncheck == 0), Ahtot, i, pref);

    try {
      succ = cyc.rest(4 * 3600, dt, ndata, Ah, Wh); //!< 1) 4h rest
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting 1 in cycle " << i << '\n';
        break;
      }

      Ahtot += std::abs(Ah);
      storeThroughput(ID_rest1, Ah, Wh, su);

      //!< 2) 1C charge
      double dtime;
      cyc.CC(-su->Cap(), su->Vmax(), tlim, dt, ndata, Ah, Wh, dtime); //!< #TODO should have a condition.

      Ahtot += std::abs(Ah);
      storeThroughput(ID_cha1, Ah, Wh, su);

      succ = cyc.rest(1 * 3600, dt, ndata, Ah, Wh); //!< 3) 1h rest
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting2 in cycle " << i << '\n';
        break;
      }
      Ahtot += std::abs(Ah);
      storeThroughput(ID_rest2, Ah, Wh, su);

      cyc.CC(su->Cap(), su->Vmin(), tlim, dt, ndata, Ah, Wh, dtime); //!< 4) 1C discharge

      Ahtot += std::abs(Ah);
      storeThroughput(ID_dis1, Ah, Wh, su);

      succ = cyc.rest(4 * 3600, dt, ndata, Ah, Wh); //!< 5) 4h rest
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting3 in cycle " << i << '\n';
        break;
      }

      Ahtot += std::abs(Ah);
      storeThroughput(ID_rest3, Ah, Wh, su);

      cyc.CC(-su->Cap() / 2.0, su->Vmax(), tlim, dt, ndata, Ah, Wh, dtime); //!< 6) 0.5C charge
      Ahtot += std::abs(Ah);
      storeThroughput(ID_cha2, Ah, Wh, su);

      succ = cyc.rest(4 * 3600, dt, ndata, Ah, Wh); //!< 7) 4h rest
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting4 in cycle " << i << '\n';
        break;
      }

      Ahtot += std::abs(Ah);
      storeThroughput(ID_rest4, Ah, Wh, su);

      succ = cyc.CC(su->Cap() / 2.0, su->Vmin(), tlim, dt, ndata, Ah, Wh, dtime); //!< 8) 0.5C discharge

      Ahtot += std::abs(Ah);
      storeThroughput(ID_dis2, Ah, Wh, su);

      //!< 9) 5h rest
      succ = cyc.rest(5 * 3600, dt, ndata, Ah, Wh);
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting5 in cycle " << i << '\n';
        break;
      }

      Ahtot += std::abs(Ah);
      storeThroughput(ID_rest5, Ah, Wh, su);
    } catch (int e) {
      std::cout << "An unexpected problem! May be CCs. See what is happening.\n";
      break;
    }
  } //!< loop cycle ageing

  //!< final checkup
  //!< This checkup will write the usage statistics, so it must be called with the actual SU and not with a copy
  //!< because usage statistics of cells and modules (or their cooling systems) are not copied over in copy()
  try {
    writeThroughput(su->getFullID(), Ahtot);
    auto task_indv = [&](int i) {
      if (i == 0)
        checkUp(su, Ahtot, Ncycle); //!< check-up of the cells #TODO if we can use the previous definition of this task_indv
      else if (i == 1)
        checkMod(su); //!< check-up of the modules
    };

    run(task_indv, 2, 2); //!< #TODO if parallel in two cores is really needed.
  } catch (int e) {
    std::cout << "Error in CycleAge when getting a capacity checkup, skipping it.\n";
  }

//!< push a write such that if cells still have cycling data, this is written
#if DATASTORE_CELL == 2
  su->writeData(pref); //!< only do if cycling data. Usage statistics are written by the checkup
#endif
}

void Procedure::storeThroughput(int ID, double Ah, double Wh, StorageUnit *su)
{
  /*
   * Function to store data of this action in the procedure.
   * Data is stored in the Vectors of this Procedure.
   * The function writethroughput will write these vectors to a file
   *
   * IN
   * ID 	integer with an indentification of what this step was
   * Ah	charge throughput of of the cells in this step
   * Wh 	energy throughput of the cells in this step
   * su 	StorageUnit after this step
   */

  //!< Data storage is slightly different depending on whether this function is done with a Cell, Module or Battery

  double coolSystemLoad{ 0 }, convloss{ 0 }; //!< Zero initialization is important.

  if (auto m = dynamic_cast<Module *>(su)) //!< pointer to module to cast su to a Module
  {                                        //!< if this is a module, add the power needed to operate the cooling system
    coolSystemLoad = m->getCoolingLoad() / 3600.0;
  } else if (auto b = dynamic_cast<Battery *>(su)) //!< pointer to module to cast su to a Battery
  {
    coolSystemLoad = b->getCoolingLoad() / 3600.0;
    convloss = b->getAndResetConvLosses() / 3600.0;
  }

  throughput.push_back({ ID, Ah, Wh, coolSystemLoad, convloss });
}

void Procedure::balanceCheckup(StorageUnit *su, bool balance, bool checkup, double Ahtot, int nrCycle, std::string pref)
{
  /*
   * Wrapper to call balancing and check-up functions
   */
  if (!unitTest) //!< Write battery state
  {
    const std::string name = su->getFullID() + "_state.csv";
    std::ofstream file(PathVar::results + name);

    if (file.is_open()) {
      std::vector<double> s; //!< #TODO if vector is removable.
      s.clear();
      su->getStates(s);

      file << nrCycle << '\n';

      for (auto s_i : s)
        file << s_i << '\n';

      file.close();
    } else {
      std::cout << "Error in Procedure::balanceCheckup. File:\n"
                << PathVar::results + name
                << " could not be opened.\n";
    }
  }
  //!< balance
  if (balance) {
    try {
      if (!unitTest)
        std::cout << "Starting balancing.\n";
      su->setBlockDegAndTherm(true);
      rebalance(su);
      su->setBlockDegAndTherm(false); //!< is skipped if error in rebalance
      if (!unitTest)
        std::cout << "Finish balancing.\n";
    } catch (int e) {
      std::cout << "Error in CycleAge when balancing the cells, skipping it.\n";
    }

    su->setBlockDegAndTherm(false); //!< ensure degradation is always turned on

    //!< #TODO if degradation should be toggled to previous state not just make it false.
  }

  //!< do a capacity checkup in a separate thread
  if (checkup) {
    //!< push a write such that if cells still have cycling data, this is written
    if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeTimeData)
      if (Ahtot > 0) {
        std::cout << "flush cycling data stored in the cells.\n";
        try {
          su->writeData(pref); //!< only do if cycling data. Usage statistics are written by the checkup
        } catch (...) {
          std::cerr << "Error when flushing the data from the cells, skipping for now\n";
          //!< note: this might result in corrupt data collection until the next write since some cells finished
          //!< writing (index = 0) while other dont (index != 0)
          //!< but since every cell writes its own data, that will not be a major problem.
          //!< worst case, it results in a writeData because the index of the cells which had not yet written data fills up
          //!< The cell which caused the error might be a bigger problem since it might have part-written results but not reset the index to 0
          //!< so it can write some duplicate data, but that can be sorted out later by comparing the time in column 7 (and removing duplicates)
        }
      }

    //!< update the statistics of the cooling system
    //!< Note that this cannot be done with a copy, so we have to do it in this thread
    if (!unitTest)
      std::cout << "Write thermal stats of the modules.\n";

    checkMod(su);

    //!< write the throughput parameters
    writeThroughput(su->getFullID(), Ahtot);

    //!< do the main checkup for the cells
    std::unique_ptr<StorageUnit> su2{ su->copy() }; //!< copy of the battery which can be used for a check-up on a separate thread
    if (!unitTest)
      std::cout << "Do a check-up on a separate thread after Ah = " << Ahtot << '\n';
    try {
      checkUp(su, Ahtot, nrCycle);

      // checkUp(su2.get(), Ahtot, nrCycle);
      //!< if you want to do it in this thread with the SU itself instead of copy to a separate thread
      //!< normally it was su #TODO ???
      //!< #TODO copy the SU to a new one (with same state), so we can do the capcheck with this copy and the main SU can keep cycling
    } catch (int e) {
      std::cout << "Error in CycleAge when getting a capacity checkup, skipping it.\n";
    }
  }
}

void Procedure::rebalance(StorageUnit *su)
{
  /*
   * Bring all lowest-level cells to the same voltage and 0 current.
   * The voltage we bring them to is given by Vset
   *
   * Note, this is a very crude balancing algorithm.
   * Normally you would go to the mean of all cells in a module, and use hierarchical balancing between modules
   */

  const double Cset = 1.0 / 2.0;    //!< use a C/2 current in the CC phase of the rebalance
  const double Clim = 1.0 / 1000.0; //!< crate for CV limit current
  double dt = 2;
  const int ndata = 0;
  double Ahi, Whi;

  //!< check if this is a p or s module, the type of the object stored at su is ... (hence *su)
  auto m = dynamic_cast<Module *>(su);
  //!< if it is a module, recursively call this function on all children to rebalance the cells
  if (m) {
    for (auto &cs_k : m->getSUs()) {
      try {
        rebalance(cs_k.get());
      } catch (int e) {
        const auto k = std::distance(m->getSUs()[0].get(), cs_k.get());
        std::cout << "Error when rebalancing of SU " << k << " whose ID = "
                  << cs_k->getFullID() << ", error " << e << ".\n";
      }
    }

    //!< if it is a parallel module, redistributeCurrent to ensure we end up with a valid voltage
    auto mp = dynamic_cast<Module_p *>(su);
    if (mp)
      mp->redistributeCurrent(true, true); //!< check the voltages are valid
  }

  //!< if it is an SPM cell,
  if (typeid(*su) == typeid(Cell_SPM)) {
    Cycler cyc(su, "rebalance"); //!< #TODO if we can do it without a cycler since it initialises a string?
    double Iset = su->Cap() * Cset;
    double Ilim = su->Cap() * Clim;
    double dtime;
    cyc.CCCV(Iset, balance_voltage, Ilim, dt, ndata, Ahi, Whi, dtime);
    auto status = su->setCurrent(0, true, true); //!< set the current of the cell to 0

    if (isStatusBad(status))
      throw 5555; //!< #TODO -> we are doing this because previously setCurrent was throwing.
  }
}

void Procedure::checkUp(StorageUnit *su, double Ah, int nrCycle)
{
  /*
   * Main function for a checkup.
   * It will call all the subfunctions to do the things we need
   * It will write one document throughout the degradation procedure
   *
   *
   * The document has one column per cell
   * The first column indicates what is on this row
   *		1			ID
   *		1			variation
   *		9			separator marker
   *		2			total battery Ah
   *		3			cell throughput
   *		4			cell capacity
   *		10-10+nstate cell state [1 to x]
   *
   * A separate document has the most recent cell statistics
   *		5			histogram I
   *		6			histogram V
   *		7			histogram T
   *
   * IN
   * su 		storage unit which needs a check-up
   * Ah		cumulative charge throughput from the total SU until now
   * nrCycle 	number of cycles done so far
   */

  //!< Name of the file which will be written
  std::string name = su->getFullID() + "_checkUp.csv";
  std::string name_stats = su->getFullID() + "_cellStats.csv";
  int separator = 9999999; //!< marker to separate different things in the document
  std::ofstream file, file_stats;

  //!< bring to correct voltage
  try {
    checkUp_prep(su);
  } catch (...) {
  };

  //!< get a vector with pointers to the cells
  std::vector<Cell *> cells;
  checkUp_getCells(su, cells);

  //!< Get the mutex to start writing in the results file
  std::lock_guard<std::mutex> writeGuard(MUTEX_procedure_checkUp); //!< locks the mutex if it is free
  //!< #TODO remove mutex somehow.
  //!< write the usage stats of all cells in a separate document
  if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeHistogramData) {
    file_stats.open(name_stats); //!< open and clear whatever was there since we want to write the most recent stats only
    try {
      checkUp_writeStats(cells, file_stats);
    } catch (...) {
    }

    file_stats.close();
  }

  //!< If this is the first checkup, start by writing the cell IDs and cell-to-cell parameters
  //!< else write the cumulative charge throughput
  if (Ah == 0) {
    file.open(name); //!< open from scratch, clear whatever was in the file before
    checkUp_writeInitial(cells, file, separator);
  } else
    file.open(name, std::ios_base::app); //!< append to the existing file

  //!< Write the charge throughput and cycle number of the entire battery (identical for all cells)
  file << 2 << ','; //!< identification for this row
  for (size_t i = 0; i < cells.size(); i++)
    file << Ah << ',';
  file << '\n'
       << 2 << ','; //!< identification for this row
  for (size_t i = 0; i < cells.size(); i++)
    file << nrCycle << ',';
  file << '\n';

  //!< do the main checkup, which writes the throughput, capacity and state of every cell
  try {
    checkUp_writeMain(cells, file);
  } catch (...) {
    std::cout << "checkUp_writeMain has throwed! We continue but worth to check probably.\n";
  }

  //!< close the file
  file.close();

  //!< The mutex will be automatically released when writeGuard goes out of scope at the end of this function
  //!< see also https://en.cppreference.com/w/cpp/thread/lock_guard

  if (!unitTest)
    std::cout << "Finishing the check-up on a separate thread after Ah = "
              << Ah << '\n';
}

void Procedure::checkUp_prep(StorageUnit *su)
{
  /*
   * Prepare the top-level SU for a capacity check.
   * We do two things
   * 		ensure we are in the middle V range
   * 			problem in parallel modules if you are <3V then the state might be illegal cause dV is very large if there was no CV phase due to steep OCV curve
   * 			and if initial state is illegal, the capcheck will fail since you always restore sini (i.e. this illegal state)
   * 		this will also allow the battery to cool down to a normal temperature
   */

  const double vwindow = (su->Vmax() - su->Vmin());
  const double V = su->Vmax() - 0.25 * vwindow; //!< go to a voltage at about 75% SOC
  const double I = -su->Cap() / 25.0;           //!< use a C/25 rate to charge, so P modules have plenty of time to equalise the voltages
  const double tlim = 999999;
  const double dt = 2;
  int ndata = 0;
  double ahi, whi;

  //!< charge to middle voltage
  Cycler cyc;
  cyc.initialise(su, "pre-checkUp");
  try {
    double dtime;
    cyc.setDiagnostic(true); //!< stop when the voltage of one cell has reached the maximum or minimum
    cyc.CC(I, V, tlim, dt, ndata, ahi, whi, dtime);
  } catch (int e) {
    std::cout << "Error in checkUp_prep when trying to charge SU " << su->getFullID()
              << " to " << V << " at I = " << I << ". Skipping.\n";
  }
}

void Procedure::checkUp_getCells(StorageUnit *su, std::vector<Cell *> &cells)
{
  /*
   * This functions makes an array with pointers to the cells
   *
   * IN
   * su			pointer to a storage unit
   * cells 		std::vector with pointers to the cells found so far
   *
   * OUT
   * cells 		append pointers to the newly found cells at the end
   *
   */

  //!< If SU is a cell, simply return a pointer to it (after casting it to a cell)
  if (auto c = dynamic_cast<Cell *>(su))
    cells.push_back(c);
  else if (auto b = dynamic_cast<Battery *>(su)) //!< if SU is a battery, recursively call this function on the module with the cells
    checkUp_getCells(b->getCells(), cells);
  else if (auto m = dynamic_cast<Module *>(su)) //!< If su is a module, recursively call this function on its children
  {
    auto &cs = m->getSUs();
    for (auto &su_ptr : cs)
      checkUp_getCells(su_ptr.get(), cells);
  }
}

void Procedure::checkUp_writeInitial(std::vector<Cell *> &cells, std::ofstream &file, int separator)
{
  /*
   * Write the cell IDs and the parameters of the cell-to-cell variations in each cell of the vector.
   *
   * This function starts the document 'name' from scratch and erases whatever was in there
   */

  int IDnumber = 1;     //!< number in first column to indicate what is on this row [ID string and variations]
  int ID_separator = 9; //!< number in first column to indicate what is on this row [separator]

  if (file.is_open()) {
    //!< first row are cell IDs
    file << IDnumber << ',';
    for (auto &cell : cells)
      file << cell->getFullID() << ',';
    file << '\n';

    //!< rows 2-5 is are the variation-parameter (capacity, resistance, SEI degradation rate, LAM degradation rate spread)
    // #TODO Print variations here!

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (auto &_ : cells)
      file << separator << ',';
    file << '\n';
  } else
    std::cerr << "Error in checkUp_writeInitial, the file is no open. Skipping this part.\n";
}

void Procedure::checkUp_writeMain(std::vector<Cell *> &cells, std::ofstream &file)
{
  /*
   * Write the throughput, capacity and cell state of every cell
   *
   * this function appends at the end of a file, so can be called multiple times
   */
  const double dt = 2;

  //!< if you want to calculate the total heat generation from all cells:
  /*double Qgen = 0;
std::shared_ptr<Cell_SPM> c;
for(auto& cell : cells){
  //!< cast to SPM cell
  c = dynamic_pointer_cast<Cell_SPM>(cell); //!< cast to a smart pointer
  Qgen += c->thermal_getTotalHeat();
}
std::cout<<"Total heat generation of all cells is "<< Qgen<<endl;*/

  if (!file.is_open()) {
    std::cerr << "Error in checkUp_writeMain, the file is not open,"
                 " could not write to it. Skipping this checkup.\n";
    return;
  }

  // Capacity checking protocol?
  auto checkCap = [this, dt](auto &cell) {
    double cap{ 0 };
    double dAh, whtot, dtime; // Dummy variables.
    //!< get every cell to the middle voltage (checkup_prep did total SU, so in series string a small cell might not be at the middle)
    auto status = cell->setCurrent(0, false, true); //!< ensure we start a check-up with a 0 current

    if (isStatusBad(status)) {
      std::cout << "Error when measuring the capacity of Cell " << cell->getFullID()
                << ", error " << getStatusMessage(status) << ". Recording as 0.\n";
      return cap;
    }

    const double Vwindow = cell->Vmax() - cell->Vmin();
    const double V = cell->Vmax() - 0.25 * Vwindow; //!< go to a voltage at about 75% SOC
    Cycler cyc(cell, "checkUp");
    cell->setBlockDegAndTherm(true);
    cyc.CCCV(cell->Cap() / 5.0, V, cell->Cap() / 5.0, dt, ndata, dAh, whtot, dtime);
    status = cell->setCurrent(0, false, true); //!< ensure we start a check-up with a 0 current

    if (isStatusBad(status)) {
      std::cout << "Error when measuring the capacity of Cell " << cell->getFullID()
                << ", error " << getStatusMessage(status) << ". Recording as 0.\n";
      return cap;
    }

    cap = cyc.testCapacity(dAh, dtime);
    cell->setBlockDegAndTherm(false);

    return cap;
  };

  const size_t start = (typeid(*cells[0]) == typeid(Cell_SPM)) ? 2 * settings::nch : 0;
  const auto nstate = cells[0]->viewStates().size();
  //!< Write the throughput of each cell
  file << "time [s]" << ',' << "Current throughput [Ah]" << ',' << "Energy throughput [Wh]"
       << "Capacity [Ah]" << ',' << "States" << '\n';
  for (auto &cell : cells) {
    const auto th = cell->getThroughputs();
    const auto cap = checkCap(cell);

    file << th.time << ',' << th.Ah << ',' << th.Wh << ','
         << cap << ',';


    const auto st_view = cell->viewStates(); //!< write the state of a cell, if an SPM cell, skip the concentration states since they don't give info
    for (size_t i{ start }; i < nstate; i++) //!< loop for each state variable (rows)
    {
      if (i != start)
        file << ',';
      file << st_view[i];
    }

    file << '\n'
         << '\n';
  }
}

void Procedure::checkUp_writeStats(std::vector<Cell *> &cells, std::ofstream &file)
{
  /*
   * Write the usage statistics of every cell
   * Remember that in copy(), usage stats are not copied across to the new cell.
   * So this function will write all zeros if called with a copy of the SU
   */

  if constexpr (true) {
    std::cout << "Not implemented: checkUp_writeStats!\n";
    return;
  }

  const int ID_histI = 5;
  const int ID_histV = 6;
  const int ID_histT = 7;
  const int ID_separator = 9;
  const int seperator = 123;

  if (file.is_open()) {
    //!< write statistics of the current
    double *histI, *histV, *histT;
    for (int i = 0; i < settings::DATASTORE_NHIST; i++) { //!< loop for each bin in the stats of I (rows)
      file << ID_histI << ',';
      for (auto &_ : cells) { //!< loop for each cell (columns)
        //!< cell->getHist(&histI, &histV, &histT, nout);
        file << histI[i] << ',';
      }
      file << '\n';
    }

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    // for (auto &_ : cells)
    //   file << separator << ',';
    file << '\n';

    //!< statistics of the voltage
    for (int i = 0; i < settings::DATASTORE_NHIST; i++) { //!< loop for each bin in the stats of V (rows)
      file << ID_histV << ',';
      for (auto &_ : cells) { //!< loop for each cell (columns)
        //	cell->getHist(&histI, &histV, &histT, nout);
        file << histV[i] << ',';
      }
      file << '\n';
    }

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    // for (auto &_ : cells)
    //   file << separator << ',';
    file << '\n';

    //!< statistics of the temperature
    for (int i = 0; i < settings::DATASTORE_NHIST; i++) { //!< loop for each bin in the stats of T (rows)
      file << ID_histT << ',';
      for (auto &_ : cells) { //!< loop for each cell (columns)
        //!< cell->getHist(&histI, &histV, &histT, nout); #TODO -> This needs to be activated sometime.
        file << histT[i] << ',';
      }
      file << '\n';
    }

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    // for (auto &_ : cells)
    //   file << separator << ',';
    file << '\n';

  } //!< if file is open
  else
    std::cerr << "Error in checkUp_writeStats, the file is not open, could not write to it. Skipping this checkup.\n";
}

void Procedure::checkMod(StorageUnit *su)
{
  /*
   * Overall function to do a check-up of the modules
   */

  //!< Name of the file which will be written
  std::string name = su->getFullID() + "_checkModules.csv";
  int separator = 9999999; //!< marker to separate different things in the document

  //!< get a vector with pointers to the modules
  std::vector<Module *> mods;
  checkMod_getModules(su, mods);

  //!< get the Battery if there is one
  auto b = dynamic_cast<Battery *>(su);

  //!< create file and write module ID names
  std::ofstream file(PathVar::results + name); //!< open from scratch, clear whatever was in the file before
  checkMod_writeInitial(mods, b, file, separator);

  //!< write the usage stats of the cooling systems
  if constexpr (settings::DATASTORE_COOL == 1) {
    try {
      checkMod_writeCoolStats(mods, b, file, separator);
    } catch (...) {
    };
  }

  //!< close the file
  file.close();

  if (!unitTest)
    std::cout << "Finishing the check-up of the modules\n";
}

void Procedure::checkMod_getModules(StorageUnit *su, std::vector<Module *> &mods)
{
  /*
   * This functions makes an array with pointers to the cells
   *
   * IN
   * su			pointer to a storage unit
   * modules 		std::vector with pointers to the modules found so far
   *
   * OUT
   * modules 		append pointers to the newly found modules at the end
   *
   */

  //!< if SU is a battery, recursively call this function on the module with the cells
  if (auto b = dynamic_cast<Battery *>(su))
    checkMod_getModules(b->getCells(), mods);

  //!< check if this SU is a module
  if (auto m = dynamic_cast<Module *>(su)) //!< cast it to a module and add a pointer to it to the vector
  {
    mods.push_back(m);

    for (auto &cs_k : m->getSUs()) //!< recursively call on its children, to check if they are modules too
      checkMod_getModules(cs_k.get(), mods);
  } //!< else, do nothing
}

void Procedure::checkMod_writeInitial(std::vector<Module *> mods, Battery *batt, std::ofstream &file, int separator)
{
  /*
   * Write the module IDs and the number of cells connected to this module
   *
   * This function starts the document 'name' from scratch and erases whatever was in there
   */

  const int IDnumber = 1;      //!< number in first column to indicate what is on this row [ID string and variations]
  const int ID_separator = 99; //!< number in first column to indicate what is on this row [separator]
  const auto ncells = mods.size() + 1;
  const bool b{ batt != nullptr }; //!< boolean indicating if there is a Battery or not

  if (file.is_open()) {
    //!< first row are module IDs
    file << IDnumber << ',';
    if (b)
      file << batt->getFullID() << ','; //!< value of total battery

    for (auto &mod : mods)
      file << mod->getFullID() << ',';

    file << '\n';

    //!< second row are number of cells connected to this module
    file << IDnumber << ',';
    if (b)
      file << batt->getCells()->getNcells() << ','; //!< value of total battery

    for (auto &mod : mods)
      file << mod->getNcells() << ',';
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';
  } else
    std::cerr << "Error in checkUp_writeInitial, the file is no open. Skipping this part.\n";
}

struct GetHistograms //!< #TODO how do we remove these helper functions in constexpr?
{
  const Histogram<> &Q{ EmptyHistogram }, &flr{ EmptyHistogram }, &E{ EmptyHistogram };
  const Histogram<> &T{ EmptyHistogram }, &Qac{ EmptyHistogram }, &Eac{ EmptyHistogram };

  template <typename T>
  GetHistograms(T &) {}

  GetHistograms(CoolSystem_HVACHist_t &data) : Qac{ data.Qac }, Eac{ data.Eac }
  {
  }

  GetHistograms(CoolSystemHist_t &data) : Q{ data.Q }, flr{ data.flr }, E{ data.E }, T{ data.T }
  {
  }
};

void Procedure::checkMod_writeCoolStats(std::vector<Module *> mods, Battery *batt, std::ofstream &file, int separator)
{
  /*
   * Write the statistics of the cooling system
   * Remember that in copy(), usage stats are not copied across to the new module.
   * So this function will write all zeros if called with a copy of the SU
   */
  const auto ncells = mods.size() + 1; //!< +1 for the total battery
  const int ID_tot = 2;
  const int ID_histT = 5;
  const int ID_histQ = 6;
  const int ID_histflr = 7;
  const int ID_histE = 8;
  const int ID_histQac = 9;
  const int ID_histEac = 10;
  const int ID_separator = 99;
  const bool b{ batt != nullptr }; //!< boolean indicating if there is a Battery or not

  if (file.is_open()) {
    //!< Write overall performance of the coolsystem (total heat evacuated or absorbed
    file << ID_tot << ',';
    if (b)
      file << batt->getCoolSystem()->getTotalTime() << ','; //!< value of total battery

    for (auto &mod : mods)                                 //!< loop for each cell (columns)
      file << mod->getCoolSystem()->getTotalTime() << ','; //!< total time active [s]

    file << '\n';

    file << ID_tot << ',';
    if (b)
      file << batt->getCoolSystem()->getHeatEvac() << ','; //!< value of total battery

    for (auto &mod : mods)                                //!< loop for each cell (columns)
      file << mod->getCoolSystem()->getHeatEvac() << ','; //!< total heat evacuated from the children of this coolsystem [J]

    file << '\n';

    file << ID_tot << ',';
    if (b)
      file << batt->getCoolSystem()->getHeatabsorbed() << ','; //!< value of total battery

    for (auto &mod : mods)                                    //!< loop for each cell (columns)
      file << mod->getCoolSystem()->getHeatabsorbed() << ','; //!< total heat absorbed in the coolant of this coolsystem [J]

    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< write statistics of the temperature
    file << ID_histT << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->coolData.tData).T << '\n';

    for (auto &mod : mods) //!< loop for each cell (columns)
      file << GetHistograms(mod->getCoolSystem()->coolData.tData).T << '\n';
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< statistics of the heat extracted from children per cell
    file << ID_histQ << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->coolData.tData).Q << '\n';

    for (auto &mod : mods) //!< loop for each cell (columns)
      file << GetHistograms(mod->getCoolSystem()->coolData.tData).Q << '\n';
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< statistics of the flow rate per cell

    file << ID_histflr << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->coolData.tData).flr << '\n';

    for (auto &mod : mods) //!< loop for each cell (columns)
      file << GetHistograms(mod->getCoolSystem()->coolData.tData).flr << '\n';
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< statistics of the operating power for the fan per cell
    file << ID_histE << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->coolData.tData).E << '\n';

    for (auto &mod : mods) //!< loop for each cell (columns)
      file << GetHistograms(mod->getCoolSystem()->coolData.tData).E << '\n';
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< Check if this module has an HVAC coolsystem. If so, also write cooling power of the AC system
    file << ID_histQac << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->HVACdata.tData).Qac << '\n';
    for (auto &mod : mods) //!< loop for each cell (columns)
    {
      if (auto c = dynamic_cast<CoolSystem_HVAC *>(mod->getCoolSystem()))
        file << GetHistograms(c->HVACdata.tData).Qac << '\n';
      else
        file << 0 << '\n';
    }
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

    //!< Check if this module has an HVAC coolsystem. If so, also write operating power of the AC system
    file << ID_histEac << ',';
    if (b)
      file << GetHistograms(batt->getCoolSystem()->HVACdata.tData).Eac << '\n';
    for (auto &mod : mods) //!< loop for each cell (columns)
    {
      if (auto c = dynamic_cast<CoolSystem_HVAC *>(mod->getCoolSystem()))
        file << GetHistograms(c->HVACdata.tData).Eac << '\n';
      else
        file << 0 << '\n';
    }
    file << '\n';

    //!< write a row with markers to indicate we have finished this section
    file << ID_separator << ',';
    for (size_t i = 0; i < ncells; i++)
      file << separator << ',';
    file << '\n';

  } //!< if file is open
  else
    std::cerr << "Error in checkUp_writeStats, the file is not open, could not write to it. Skipping this checkup.\n";
}

void Procedure::writeThroughput(std::string SUID, double Ahtot)
{
  /*
   * Write a file with the throughput per action.
   * Then clear those entries from the vector to avoid it getting too long
   *
   * The first time, the file is written from scratch (i.e. clears any existing content in this file)
   * On later times, the new content is added.
   *
   * IN
   * SUID 	ID string of the storageUnit (su->getFullID())
   * Ahtot 	total charge throughput until now
   */

  //!< Open the file
  std::string name = SUID + "_throughput.csv";
  std::ofstream file;
  if (Ahtot == 0) {
    file.open(PathVar::results + name); //!< open from scratch, clear whatever was in the file before
    file << "ID number" << ',' << "cells charge throughput" << ','
         << "Cells energy throughput" << ','
         << "total energy to operate the thermal management system" << ','
         << "losses in the converter" << '\n';
  } else
    file.open(PathVar::results + name, std::ios_base::app); //!< append to the existing file

  //!< Write the data
  for (const auto &th : throughput)
    file << th.ID << ',' << th.charge << ',' << th.energy << ','
         << th.coolSystemLoad << ',' << th.convloss << '\n';

  //!< close the file
  file.close();

  //!< erase the vectors
  throughput.clear();
}
} // namespace slide