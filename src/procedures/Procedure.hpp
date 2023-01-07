/*
 * Procedure.hpp
 *
 *  Created on: 3 Mar 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "Cycler.hpp"
#include "../cells/Cell.hpp"
#include "../modules/Module.hpp"
#include "../system/Battery.hpp"
#include "../types/data_storage/cell_data.hpp"

#include <vector>
#include <iostream>
#include <fstream>

namespace slide {

class Procedure
{
protected:
  bool balance{ true };
  bool unitTest{ false };
  int ndata{ 0 };
  double balance_voltage{ 3.65 };

  std::vector<ProcedureThroughputData> throughput;
  void storeThroughput(ThroughputData th, StorageUnit *su);

public:
  Procedure() = default;
  Procedure(bool balance, double Vbal, int ndata, bool unitTest = false);

  ~Procedure() = default;

  void cycleAge(StorageUnit *su, bool testCV);
  void cycleAge(StorageUnit *su, int Ncycle, int Ncheck, int Nbal, bool testCV, double Ccha, double Cdis, double Vmax, double Vmin);
  void useCaseAge(StorageUnit *su, int cool);

  //!< function calls
  void balanceCheckup(StorageUnit *su, bool balance, bool checkup, double Ah, int nrCycle, std::string pref);

  //!< balancing
  Status rebalance(StorageUnit *su);

  //!< check-up procedures
  void checkUp(StorageUnit *su, double Ah, int nrCycle); //!< main checkup function which will call the others
  void checkUp_prep(StorageUnit *su);                    //!< bring the SU to a good voltage
  void checkUp_writeInitial(std::vector<Cell *> &cells, std::ofstream &file);

  //!< cooling system procedures
  void checkMod(StorageUnit *su); //!< main function to do a checkup on the modules

  //!< write the charge and energy throughput
  void writeThroughput(const std::string &SUID, double Ahtot);
};
} // namespace slide