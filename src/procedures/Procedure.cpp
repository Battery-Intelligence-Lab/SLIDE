/*
 * Procedure.cpp
 *
 * Created on: 3 Mar 2020
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */
#include "Procedure.hpp"
#include "procedure_util.hpp"
#include "Cycler.hpp"
#include "../cells/cells.hpp"
#include "../modules/modules.hpp"
#include "../system/Battery.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ctime>
#include <typeinfo>
#include <vector>
#include <memory>

namespace slide {

struct GetHistograms //!< #TODO how do we remove these helper functions in constexpr?
{
  const Histogram<> &Q{ EmptyHistogram }, &flr{ EmptyHistogram }, &E{ EmptyHistogram };
  const Histogram<> &T{ EmptyHistogram }, &Qac{ EmptyHistogram }, &Eac{ EmptyHistogram };

  template <typename T>
  GetHistograms(T &) {}
  GetHistograms(CoolSystem_HVACHist_t &data) : Qac{ data.Qac }, Eac{ data.Eac } {}
  GetHistograms(CoolSystemHist_t &data) : Q{ data.Q }, flr{ data.flr }, E{ data.E }, T{ data.T } {}
};

Procedure::Procedure(bool balance_, double Vbal_, int ndata_, bool unitTest_)
  : balance{ balance_ }, unitTest{ unitTest_ }, ndata{ ndata_ }, balance_voltage{ Vbal_ } {}

void Procedure::cycleAge(StorageUnit *su, int Ncycle, int Ncheck, int Nbal, bool testCV, double Ccha, double Cdis, double Vmax, double Vmin)
{
  /*
   * Simulate a cycle ageing experiment with the given parameters
   */
  const std::string pref{ "cycleAge" }; //!< prefix appended at the start of the names of all files
  constexpr auto diagnostic = true;     //!< stop (dis)charging when one of the cells reaches a voltage limit
  auto cyc = Cycler(su, pref).setDiagnostic(diagnostic);

  //!< variables
  constexpr double dt = 2;
  const double Ilim{ 0.1 }, Icha{ -Ccha * su->Cap() }, Idis{ Cdis * su->Cap() };
  Status succ{};
  ThroughputData th{};

  //!< Make a clock to measure how long the simulation takes
  Clock clk{};

  //!< loop for cycle ageing
  for (int i = 0; i < Ncycle; i++) {
    if (!unitTest)
      std::cout << "SU " << su->getFullID() << " starting loop iteration " << i << " after "
                << clk << " with V = " << su->V() << ", T = " << K_to_Celsius(su->T())
                << " and hot spot T = " << K_to_Celsius(su->getThotSpot()) << '\n';

    //!< Balance (in this thread) and do a check-up (on a separate thread) if needed
    balanceCheckup(su, balance && (i % Nbal == 0), (i % Ncheck == 0), th.Ah(), i, pref);

    const auto th_before = th.Ah();
    succ = cyc.CC(Icha, Vmax, TIME_INF, dt, ndata, th); //!< CC charge

    if (!isVoltageLimitReached(succ)) {
      std::cout << "Error in CycleAge when charging in cycle "
                << i << ", stop cycling. " << getStatusMessage(succ);
      break;
    }

    storeThroughput(th, su); //!< increase the throughput parameters

    if ((th.Ah() - th_before) < 1e-5) {
      std::cerr << "Error in Procedure::cycleAge, we didn't manage to charge "
                << "any energy. Stopping the cycling after " << i << " cycles.\n";
      break;
    }

    //!< CV charge
    if (testCV) {
      succ = cyc.CV(Vmax, Ilim, TIME_INF, dt, ndata, th);
      if (!isCurrentLimitReached(succ)) {
        std::cout << "Error in CycleAge when CV charging in cycle "
                  << i << ", stopping cycling.\n"
                  << getStatusMessage(succ);
        break;
      }

      storeThroughput(th, su); // #TODO why do we do this?
    }


    succ = cyc.CC(Idis, Vmin, TIME_INF, dt, ndata, th); //!< CC discharge
    storeThroughput(th, su);
    if (!isVoltageLimitReached(succ)) {
      std::cout << "Error in CycleAge when discharging in cycle " << i << ", stop cycling.\n";
      break;
    } // #TODO depending on diagnostics a bad status may be skipped.

    //!< CV discharge
    if (testCV) {
      //!< reset in case error in CC and Ah gets not changed (without reset it would keep its old value)
      succ = cyc.CV(Vmin, Ilim, TIME_INF, dt, ndata, th);
      if (!isCurrentLimitReached(succ)) { //!< #TODO -> if diagnostic on the individual cell limits are respected. Otherwise system level.
        std::cout << "Error in CycleAge when CV discharging in cycle " << i << ", stop cycling.\n";
        break;
      }
      storeThroughput(th, su);
    }

  } //!< loop cycle ageing

  //!< final checkup
  //!< This checkup will write the usage statistics, so it must be called with the actual SU and not with a copy
  //!< because usage statistics of cells and modules (or their cooling systems) are not copied over in copy()
  if (!unitTest) std::cout << "Doing a final checkup after " << th.Ah() << " Ah.\n";
  writeThroughput(su->getFullID(), th.Ah());
  checkUp(su, th.Ah(), Ncycle); //!< check-up of the cells
  checkMod(su);                 //!< check-up of the modules

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
  const unsigned Ncycle = 15;  // #TODO it should be 15000
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
  constexpr bool diagnostic = true;  //!< stop (dis)charging when one of the cells reaches a voltage limit
  const std::string pref = "useAge"; //!< prefix appended at the start of the names of all files
  constexpr double dt = 1;
  constexpr unsigned Ncycle = 20; //!< 365 * 10; //!< 10 year
  constexpr unsigned Ncheck = 10; //!< do a checkup ever month
  constexpr unsigned Nbal = 7;    //!< balance every week

  //!< Variables
  auto cyc = Cycler(su, pref).setDiagnostic(diagnostic);
  ThroughputData th{};
  Status succ;

  //!< fully discharge the SU at a C/2
  succ = cyc.CC(su->Cap() / 2.0, su->Vmin(), TIME_INF, dt, ndata, th);

  if (!isLimitsReached(succ))
    std::cout << "Error in useAge when initially discharging the battery, continue as normal.\n"
              << getStatusMessage(succ) << '\n';

  //!< Make a clock to measure how long the simulation takes
  Clock clk;

  //!< loop for use case ageing
  for (unsigned i = 0; i < Ncycle; i++) {
    if (!unitTest)
      std::cout << "SU " << su->getFullID() << " starting loop iteration " << i << " after "
                << clk << " with V = " << su->V() << ", T = " << K_to_Celsius(su->T())
                << " and hot spot T = " << K_to_Celsius(su->getThotSpot()) << '\n';

    //!< Balance (in this thread) and do a check-up (on a separate thread) if needed
    balanceCheckup(su, balance && (i % Nbal == 0), (i % Ncheck == 0), th.Ah(), i, pref);

    try {
      succ = cyc.rest(4 * 3600, dt, ndata, th); //!< 1) 4h rest
      storeThroughput(th, su);                  // #TODO why do we store this!
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting 1 in cycle " << i << '\n';
        break;
      }

      //!< 2) 1C charge
      cyc.CC(-su->Cap(), su->Vmax(), TIME_INF, dt, ndata, th); //!< #TODO should have a condition.
      storeThroughput(th, su);

      succ = cyc.rest(1 * 3600, dt, ndata, th); //!< 3) 1h rest
      storeThroughput(th, su);
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting2 in cycle " << i << '\n';
        break;
      }

      cyc.CC(su->Cap(), su->Vmin(), TIME_INF, dt, ndata, th); //!< 4) 1C discharge
      storeThroughput(th, su);                                // #TODO sometimes we don't check limits.

      succ = cyc.rest(4 * 3600, dt, ndata, th); //!< 5) 4h rest
      storeThroughput(th, su);
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting3 in cycle " << i << '\n';
        break;
      }

      cyc.CC(-su->Cap() / 2.0, su->Vmax(), TIME_INF, dt, ndata, th); //!< 6) 0.5C charge
      storeThroughput(th, su);

      succ = cyc.rest(4 * 3600, dt, ndata, th); //!< 7) 4h rest
      storeThroughput(th, su);
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting4 in cycle " << i << '\n';
        break;
      }

      succ = cyc.CC(su->Cap() / 2.0, su->Vmin(), TIME_INF, dt, ndata, th); //!< 8) 0.5C discharge
      storeThroughput(th, su);
      if (succ != Status::ReachedVoltageLimit) {
        std::cout << "Error in useAge when 0.5C discharge in cycle " << i << '\n';
        break;
      }

      //!< 9) 5h rest
      succ = cyc.rest(5 * 3600, dt, ndata, th);
      storeThroughput(th, su);
      if (succ != Status::ReachedTimeLimit) {
        std::cout << "Error in useAge when resting5 in cycle " << i << '\n';
        break;
      }

    } catch (...) {
      std::cout << "An unexpected problem! May be CCs. See what is happening.\n";
      break;
    }
  } //!< loop cycle ageing

  //!< final checkup
  //!< This checkup will write the usage statistics, so it must be called with the actual SU and not with a copy
  //!< because usage statistics of cells and modules (or their cooling systems) are not copied over in copy()
  writeThroughput(su->getFullID(), th.Ah());
  checkUp(su, th.Ah(), Ncycle); //!< check-up of the cells #TODO if we can use the previous definition of this task_indv
  checkMod(su);                 //!< check-up of the modules

  //!< push a write such that if cells still have cycling data, this is written
  if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeTimeData)
    su->writeData(pref); //!< only do if cycling data. Usage statistics are written by the checkup
}

void Procedure::storeThroughput(ThroughputData th, StorageUnit *su)
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

  if (auto m = dynamic_cast<Module *>(su)) //!< if this is a module, add the power needed to operate the cooling system
    coolSystemLoad = m->getCoolingLoad() / 3600.0;
  else if (auto b = dynamic_cast<Battery *>(su)) //!< pointer to module to cast su to a Battery
  {
    coolSystemLoad = b->getCoolingLoad() / 3600.0;
    convloss = b->getAndResetConvLosses() / 3600.0;
  }

  throughput.push_back({ th.Ah(), th.Wh(), coolSystemLoad, convloss });
}

void Procedure::balanceCheckup(StorageUnit *su, bool balance, bool checkup, double Ahtot, int nrCycle, std::string pref)
{
  /*
   * Wrapper to call balancing and check-up functions
   */
  if (!unitTest) //!< Write battery state
  {
    const std::string name = su->getFullID() + "_state.csv";
    std::ofstream file(PathVar::results / name);

    if (file.is_open()) {
      std::vector<double> s; //!< #TODO if vector is removable.
      s.clear();
      su->getStates(s);

      file << nrCycle << '\n';

      for (auto s_i : s) file << s_i << '\n';

      file.close();
    } else {
      std::cout << "Error in Procedure::balanceCheckup. File:\n"
                << PathVar::results / name << " could not be opened.\n";
    }
  }
  //!< balance
  if (balance) {
    if (!unitTest) std::cout << "Starting balancing.\n";
    su->setBlockDegAndTherm(true);
    const auto status = rebalance(su);
    if (!unitTest) std::cout << "Finish balancing.\n";

    if (isStatusBad(status))
      std::cout << "Error in CycleAge when balancing the cells, skipping it.\n";

    su->setBlockDegAndTherm(false); //!< ensure degradation is always turned on

    //!< #TODO degradation should be toggled to previous state not just make it false.
  }

  //!< do a capacity checkup in a separate thread
  if (checkup) {
    //!< push a write such that if cells still have cycling data, this is written
    if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeTimeData)
      if (Ahtot > 0) {
        su->writeData(pref); //!< only do if cycling data. Usage statistics are written by the checkup
        try {
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
    if (!unitTest) std::cout << "Write thermal stats of the modules.\n";

    checkMod(su);
    //!< write the throughput parameters
    writeThroughput(su->getFullID(), Ahtot);
    checkUp(su, Ahtot, nrCycle);
  }
}

Status Procedure::rebalance(StorageUnit *su)
{
  /*
   * Bring all lowest-level cells to the same voltage and 0 current.
   * The voltage we bring them to is given by Vset
   *
   * Note, this is a very crude balancing algorithm.
   * Normally you would go to the mean of all cells in a module, and use hierarchical balancing between modules
   *
   */
  constexpr double Cset = 1.0 / 2.0;    //!< use a C/2 current in the CC phase of the rebalance
  constexpr double Clim = 1.0 / 1000.0; //!< crate for CV limit current
  constexpr double dt = 1;
  constexpr int ndata = 0;

  //!< if it is a module, recursively call this function on all children to rebalance the cells
  if (auto m = dynamic_cast<Module *>(su)) {
    for (auto &cs_k : m->getSUs()) {
      const auto status = rebalance(cs_k.get());
      if (!isStatusSuccessful(status)) {
        const auto k = std::distance(m->getSUs()[0].get(), cs_k.get());
        std::cout << "Error when rebalancing of SU " << k << " whose ID = "
                  << cs_k->getFullID() << ", error " << getStatusMessage(status) << ".\n";

        return status;
      }
    }

  } else if (typeid(*su) == typeid(Cell_SPM)) { // #TODO why only SPM?
    const double Iset = su->Cap() * Cset;
    const double Ilim = su->Cap() * Clim;
    ThroughputData th{};
    Cycler(su).CCCV(Iset, balance_voltage, Ilim, dt, ndata, th); //!< #TODO if we can do it without a cycler since it initialises a string?
    auto status = su->setCurrent(0, true, true);                 //!< set the current of the cell to 0

    if (isStatusBad(status)) return status; //!< #TODO -> we are doing this because previously setCurrent was throwing.
  }

  return Status::Success;
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
  std::ofstream file, file_stats;

  checkUp_prep(su); //!< bring to correct voltage

  //!< get a vector with pointers to the cells
  std::vector<Cell *> cells;

  auto copyCells = [&cells](auto *su_now) {
    if (auto c = dynamic_cast<Cell *>(su_now))
      cells.push_back(static_cast<Cell *>(c->copy()));
  };

  visit_SUs(su, copyCells);

  //!< write the usage stats of all cells in a separate document
  // if constexpr (settings::DATASTORE_CELL == settings::cellDataStorageLevel::storeHistogramData)
  //  #TODO normally we were writing cell stats here but not necessary I think.
  //  open and clear whatever was there since we want to write the most recent stats only

  //!< do the main checkup, which writes the throughput, capacity and state of every cell
  const size_t start = (typeid(*cells[0]) == typeid(Cell_SPM)) ? 2 * settings::nch : 0; // #TODO : Consider having different type of cells.
  //!< Write the throughput of each cell

  //!< If this is the first checkup, start by writing the cell IDs and cell-to-cell parameters
  //!< else write the cumulative charge throughput
  file.open(name, std::ios_base::app); //!< open from scratch, clear whatever was in the file before

  file << "ID,var_cap,var_R,var_degSEI,var_degLAM,Ah,CycleNumber";
  file << "time [s],Current throughput [Ah],Energy throughput [Wh],Capacity [Ah],States\n";

  for (auto &cell : cells) {
    file << cell->getFullID() << ','; //!< first column is cell IDs

    //!< columns 2-5 is are the variation-parameter (capacity, resistance, SEI degradation rate, LAM degradation rate spread)
    for (auto var : cell->getVariations())
      file << var << ',';


    //!< Write the charge throughput and cycle number of the entire battery (identical for all cells)
    file << Ah << ',' << nrCycle << ',';

    const auto th = cell->getThroughputs();       // #TODO if throughputs are saved then they should be included in states. Therefore this should not be needed.
    const auto cap = Cycler(cell).testCapacity(); // #TODO this should be done in parallel.

    file << th.time() << ',' << th.Ah() << ',' << th.Wh() << ',' << cap;

    const auto st_view = cell->viewStates();         //!< write the state of a cell, if an SPM cell, skip the concentration states since they don't give info
    for (size_t i{ start }; i < st_view.size(); i++) //!< loop for each state variable (rows)
      file << ',' << st_view[i];

    file << '\n'
         << '\n';
  }

  //!< close the file
  file.close();

  if (!unitTest) std::cout << "Finishing the check-up after Ah = " << Ah << '\n';
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
  const double V = 0.75 * su->Vmax() + 0.25 * su->Vmin(); //!< go to a voltage at about 75% SOC
  const double I = -su->Cap() / 25.0;                     //!< use a C/25 rate to charge, so P modules have plenty of time to equalise the voltages
  const double dt = 2;
  int ndata = 0;
  double ahi, whi;

  //!< charge to middle voltage, stop when the voltage of one cell has reached the maximum or minimum
  ThroughputData th{};
  Cycler(su, "pre-checkUp").setDiagnostic(true).CC(I, V, TIME_INF, dt, ndata, th);
}

void Procedure::checkMod(StorageUnit *su)
{
  /*
   * Overall function to do a check-up of the modules
   */

  //!< Name of the file which will be written
  std::string name_overall = su->getFullID() + "_checkModules_overall.csv";
  std::string name_hist = su->getFullID() + "_checkModules_histograms.csv";

  //!< create file and write module ID names
  std::ofstream file_overall(PathVar::results / name_overall); //!< open from scratch, clear whatever was in the file before
  std::ofstream file_histograms(PathVar::results / name_hist); //!< open from scratch, clear whatever was in the file before

  if (!file_overall.is_open())
    std::cerr << "Error in Procedure::checkMod, the file " << PathVar::results / name_overall
              << " is not open!. Skipping this checkup.\n";

  if (!file_histograms.is_open())
    std::cerr << "Error in Procedure::checkMod, the file " << PathVar::results / name_hist
              << " is not open!. Skipping this checkup.\n";

  file_overall << "Full ID,";              //!< Module IDs
  file_overall << "Number of cells,";      //!< Number of cells inside SU
  file_overall << "Total time,";           //!< total time active [s]
  file_overall << "Total evacuated heat,"; //!< total heat evacuated from the children of this coolsystem [J]
  file_overall << "Total absorbed heat\n"; //!< total heat absorbed in the coolant of this coolsystem [J]

  auto writeModData = [&file_overall, &file_histograms](auto *su_now) {
    const auto su_fullID = su_now->getFullID();

    file_overall << su_fullID << ',';
    file_overall << su_now->getNcells() << ',';
    file_overall << su_now->getCoolSystem()->getTotalTime() << ',';
    file_overall << su_now->getCoolSystem()->getHeatEvac() << ',';
    file_overall << su_now->getCoolSystem()->getHeatabsorbed() << '\n';

    if constexpr (settings::DATASTORE_COOL == 1) {
      // Temperature statistics:
      file_histograms << "Temperature statistics:," << su_fullID << '\n';
      file_histograms << "Time histogram:\n";
      file_histograms << GetHistograms(su_now->getCoolSystem()->coolData.tData).T << '\n';
      file_histograms << "Cooling power per cell histogram:\n";
      file_histograms << GetHistograms(su_now->getCoolSystem()->coolData.tData).Q << '\n';
      file_histograms << "Flow rate per cell histogram:\n";
      file_histograms << GetHistograms(su_now->getCoolSystem()->coolData.tData).flr << '\n';
      file_histograms << "Fan power per cell histogram:\n";
      file_histograms << GetHistograms(su_now->getCoolSystem()->coolData.tData).E << '\n';

      if (auto c = dynamic_cast<CoolSystem_HVAC *>(su_now->getCoolSystem())) {
        file_histograms << "Qac histogram:\n";
        file_histograms << GetHistograms(c->HVACdata.tData).Qac << '\n';

        file_histograms << "Eac histogram:\n";
        file_histograms << GetHistograms(c->HVACdata.tData).Qac << '\n';
      }
    }
  };

  // #TODO eliminate use of two lambdas.
  // #TODO move the above logic into data storage so polymorphism takes care of pointer cast
  auto writeModMain = [&](auto *su_now) {
    if (auto *b = dynamic_cast<Battery *>(su_now))
      writeModData(b);
    else if (auto *m = dynamic_cast<Module *>(su_now))
      writeModData(m);
  };

  visit_SUs(su, writeModMain);

  //!< close the files
  file_overall.close();
  file_histograms.close();

  if (!unitTest) std::cout << "Finished checking modules.\n";
}

void Procedure::writeThroughput(const std::string &SUID, double Ahtot)
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
    file.open(PathVar::results / name); //!< open from scratch, clear whatever was in the file before
    file << "ID number" << ',' << "cells charge throughput" << ','
         << "Cells energy throughput" << ','
         << "total energy to operate the thermal management system" << ','
         << "losses in the converter" << '\n';
  } else
    file.open(PathVar::results / name, std::ios_base::app); //!< append to the existing file

  //!< Write the data
  for (const auto &th : throughput)
    file << th.charge << ',' << th.energy << ','
         << th.coolSystemLoad << ',' << th.convloss << '\n';

  //!< close the file
  file.close();

  //!< erase the vectors
  throughput.clear();
}
} // namespace slide