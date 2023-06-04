/*
 * running_Cell_Bucket.hpp
 *
 *  Benchmark file for Cell_Bucket
 *
 *  Created on: 07 Aug 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>

namespace slide::benchmarks {

inline void run_LP_case_SmallPack()
{
  constexpr double T_ENV = 15.0_degC;

  std::string ID = "PyBAMM_LP_SmallPack"; // + std::to_string(Crate) + '_'

  auto c = Cell_SPM();
  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;
  c.setBlockDegAndTherm(true);

  Deep_ptr<StorageUnit> cs[] = { Deep_ptr<StorageUnit>(c.copy()),
                                 Deep_ptr<StorageUnit>(c.copy()) };

  auto module = Module_p("SmallPack", T_ENV, true, false, std::size(cs), 1, 1);
  module.setSUs(cs, false);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  cyc.CC(2, 2.7, 5 * 60, 1, 10, th);
  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_LP_case_MediumPack()
{
  constexpr double T_ENV = 15.0_degC;

  std::string ID = "PyBAMM_LP_MediumPack"; // + std::to_string(Crate) + '_'

  auto c = Cell_SPM();
  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;
  c.setBlockDegAndTherm(true);


  //!< Make the battery
  constexpr size_t ns = 10; //!< number of cells in modules
  constexpr size_t np = 32; //!< number of modules in strings
  constexpr bool checkCells = true;

  const double Rc_s = 1e-4;
  const double Rc_p = 2e-4;

  Deep_ptr<StorageUnit> MinS[np]; //!< array with modules in one string
  Deep_ptr<StorageUnit> CinM[ns]; //!< array with cells in one module

  double Rc2[np], Rc3[ns];               //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t is = 0; is < np; is++) {   //!< loop for modules in string
    for (size_t ic = 0; ic < ns; ic++) { //!< loop for cells in modules
      CinM[ic] = Deep_ptr<StorageUnit>(c.copy());
      Rc3[ic] = Rc_s; //!< in series module, so every cell has a resistance of Rc
    }
    auto mi = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ns, 1, 1); //!< print warning messages, single-threaded

    mi->setSUs(CinM, checkCells, true);
    mi->setRcontact(Rc3);
    MinS[is] = std::move(mi);
    Rc2[is] = Rc_p; //!< in series module, so every cell has a resistance of Rc
  }

  auto module = Module_p("MediumPack", T_ENV, true, false, ns * np, 1, 1);
  module.setSUs(MinS, false);
  module.setRcontact(Rc2);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  cyc.CC(32, 2.7, 5 * 60, 1, 10, th);
  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_LP_case_LargePack()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "PyBAMM_LP_LargePack"; // + std::to_string(Crate) + '_'

  auto c = Cell_SPM();
  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;
  c.setBlockDegAndTherm(true);


  //!< Make the battery
  constexpr size_t ns = 10; //!< number of cells in modules
  constexpr size_t np = 64; //!< number of modules in strings
  constexpr bool checkCells = true;

  const double Rc_s = 1e-4;
  const double Rc_p = 2e-4;

  Deep_ptr<StorageUnit> MinS[np]; //!< array with modules in one string
  Deep_ptr<StorageUnit> CinM[ns]; //!< array with cells in one module

  double Rc2[np], Rc3[ns];               //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t is = 0; is < np; is++) {   //!< loop for modules in string
    for (size_t ic = 0; ic < ns; ic++) { //!< loop for cells in modules
      CinM[ic] = Deep_ptr<StorageUnit>(c.copy());
      Rc3[ic] = Rc_s; //!< in series module, so every cell has a resistance of Rc
    }
    auto mi = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ns, 1, 1); //!< print warning messages, single-threaded

    mi->setSUs(CinM, checkCells, true);
    mi->setRcontact(Rc3);
    MinS[is] = std::move(mi);
    Rc2[is] = Rc_p; //!< in series module, so every cell has a resistance of Rc
  }

  auto module = Module_p("LargePack", T_ENV, true, false, ns * np, 1, 1);
  module.setSUs(MinS, false);
  module.setRcontact(Rc2);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  cyc.CC(64, 2.7, 5 * 60, 1, 10, th);
  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_LP_case_LargePackLong()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "PyBAMM_LP_LargePackLong"; // + std::to_string(Crate) + '_'

  auto c = Cell_SPM();
  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;
  c.setBlockDegAndTherm(true);


  //!< Make the battery
  constexpr size_t ns = 10; //!< number of cells in modules
  constexpr size_t np = 64; //!< number of modules in strings
  constexpr bool checkCells = true;

  const double Rc_s = 1e-4;
  const double Rc_p = 2e-4;

  Deep_ptr<StorageUnit> MinS[np]; //!< array with modules in one string
  Deep_ptr<StorageUnit> CinM[ns]; //!< array with cells in one module

  double Rc2[np], Rc3[ns];               //!< arrays with the contact resistances for the different levels (3 = lowest level)
  for (size_t is = 0; is < np; is++) {   //!< loop for modules in string
    for (size_t ic = 0; ic < ns; ic++) { //!< loop for cells in modules
      CinM[ic] = Deep_ptr<StorageUnit>(c.copy());
      Rc3[ic] = Rc_s; //!< in series module, so every cell has a resistance of Rc
    }
    auto mi = make<Module_s>("s" + std::to_string(is), T_ENV, true, false, ns, 1, 1); //!< print warning messages, single-threaded

    mi->setSUs(CinM, checkCells, true);
    mi->setRcontact(Rc3);
    MinS[is] = std::move(mi);
    Rc2[is] = Rc_p; //!< in series module, so every cell has a resistance of Rc
  }

  auto module = Module_p("LargePackLong", T_ENV, true, false, ns * np, 1, 1);
  module.setSUs(MinS, false);
  module.setRcontact(Rc2);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  for (size_t i_repeat{}; i_repeat < 3; i_repeat++) {
    cyc.CC(-64, 2.7, 20 * 60, 1, 10, th);
    cyc.rest(15 * 60, 1, 10, th);
    cyc.CC(64, 2.7, 20 * 60, 1, 10, th);
    cyc.rest(30 * 60, 1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}


} // namespace slide::benchmarks