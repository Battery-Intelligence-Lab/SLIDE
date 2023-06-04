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

inline void run_Cell_ECM_single_default_pulse()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_single_default_pulse"; // + std::to_string(Crate) + '_'
  auto c = Cell_ECM();
  c.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&c, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CC(16, 2.7, 5 * 60, 0.1, 5, th);
    cyc.CC(-16, 4.2, 5 * 60, 0.1, 5, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_Cell_ECM_single_default_CCCV()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_single_default_CCCV"; // + std::to_string(Crate) + '_'
  auto c = Cell_ECM();
  c.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&c, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CCCV(16, 2.7, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
    cyc.CCCV(16, 4.2, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}


inline void run_Cell_ECM_2_RC_single_default_pulse()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_2_RC_single_default_pulse"; // + std::to_string(Crate) + '_'
  auto c = Cell_ECM<2>();
  c.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&c, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CC(16, 2.7, 5 * 60, 0.1, 5, th);
    cyc.CC(-16, 4.2, 5 * 60, 0.1, 5, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_Cell_ECM_2_RC_single_default_CCCV()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_2_RC_single_default_CCCV"; // + std::to_string(Crate) + '_'
  auto c = Cell_ECM<2>();
  c.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&c, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CCCV(16, 2.7, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
    cyc.CCCV(16, 4.2, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}


inline void run_Cell_ECM_parallel_3_default_pulse()
{
  constexpr double T_ENV = 15.0_degC;
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_parallel_3_default_pulse"; // + std::to_string(Crate) + '_'

  double capin{ 16 }, SOCin{ 0.5 }, Rdc_{ 2e-3 };
  constexpr double Cp0 = 38e3; // first parallel capacitance
  constexpr double Rp0 = 15.8e-3;
  std::array<double, 1> Rp_{ Rp0 }, inv_tau{ 1.0 / (Rp0 * Cp0) };

  Deep_ptr<StorageUnit> cs[3] = {
    make<Cell_ECM<1>>("1", capin, SOCin, 1e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("2", capin, SOCin, Rdc_, Rp_, inv_tau),
    make<Cell_ECM<1>>("3", capin, SOCin, Rdc_, Rp_, inv_tau)
  };

  auto module = Module_p("Par3", T_ENV, true, false, std::size(cs), 1, 1);
  module.setSUs(cs, false);

  module.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);


  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CC(16, 2.7, 5 * 60, 0.1, 5, th);
    cyc.CC(-16, 4.2, 5 * 60, 0.1, 5, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_Cell_ECM_parallel_3_default_CCCV()
{
  constexpr double T_ENV = 15.0_degC;

  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_parallel_3_default_CCCV"; // + std::to_string(Crate) + '_'
  double capin{ 16 }, SOCin{ 0.5 }, Rdc_{ 2e-3 };
  constexpr double Cp0 = 38e3; // first parallel capacitance
  constexpr double Rp0 = 15.8e-3;
  std::array<double, 1> Rp_{ Rp0 }, inv_tau{ 1.0 / (Rp0 * Cp0) };

  Deep_ptr<StorageUnit> cs[] = {
    make<Cell_ECM<1>>("1", capin, SOCin, 1e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("2", capin, SOCin, 3e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("3", capin, SOCin, Rdc_, Rp_, inv_tau)
  };

  auto module = Module_p("Par3", T_ENV, true, false, std::size(cs), 1, 1);
  module.setSUs(cs, false);

  module.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CCCV(16, 2.7, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
    cyc.CCCV(16, 4.2, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}

inline void run_Cell_ECM_parallel_3_withRcontact_CCCV()
{
  constexpr double T_ENV = 15.0_degC;
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_parallel_3_withRcontact_CCCV"; // + std::to_string(Crate) + '_'
  double capin{ 16 }, SOCin{ 0.5 }, Rdc_{ 2e-3 };
  constexpr double Cp0 = 38e3; // first parallel capacitance
  constexpr double Rp0 = 15.8e-3;
  std::array<double, 1> Rp_{ Rp0 }, inv_tau{ 1.0 / (Rp0 * Cp0) };

  std::vector<double> Rcontact{ 0.5e-3, 1e-3, 0.7e-3 };

  Deep_ptr<StorageUnit> cs[3] = {
    make<Cell_ECM<1>>("1", capin, SOCin, 1e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("2", capin, SOCin, 3e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("3", capin, SOCin, Rdc_, Rp_, inv_tau)
  };

  auto module = Module_p("Par3", T_ENV, true, false, std::size(cs), 1, 1);
  module.setSUs(cs, false);
  module.setRcontact(Rcontact);

  module.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);

  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CCCV(16, 2.7, 50e-3, 0.1, 1, th);
    cyc.rest(10 * 60, 0.1, 1, th);
    cyc.CCCV(16, 4.2, 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}


inline void run_Cell_ECM_series_3_withRcontact_CCCV()
{
  constexpr double T_ENV = 15.0_degC;
  // Benchmark with default parameters:
  std::string ID = "Cell_ECM_series_3_withRcontact_CCCV"; // + std::to_string(Crate) + '_'
  double capin{ 16 }, SOCin{ 0.5 }, Rdc_{ 2e-3 };
  constexpr double Cp0 = 38e3; // first parallel capacitance
  constexpr double Rp0 = 15.8e-3;
  std::array<double, 1> Rp_{ Rp0 }, inv_tau{ 1.0 / (Rp0 * Cp0) };

  std::vector<double> Rcontact{ 0.5e-3, 1e-3, 0.7e-3 };


  Deep_ptr<StorageUnit> cs[3] = {
    make<Cell_ECM<1>>("1", capin, SOCin, 1e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("2", capin, SOCin, 3e-3, Rp_, inv_tau),
    make<Cell_ECM<1>>("3", capin, SOCin, Rdc_, Rp_, inv_tau)
  };

  auto module = Module_s("Ser3", T_ENV, true, false, std::size(cs), 1, 1);
  module.setSUs(cs, false);
  module.setRcontact(Rcontact);

  module.setBlockDegAndTherm(true);

  ThroughputData th{};
  auto cyc = Cycler(&module, ID);


  Clock clk;
  constexpr size_t Nrepeat = 3;
  for (size_t i = 0; i < Nrepeat; i++) {
    cyc.CCCV(16, module.Vmin(), 50e-3, 0.1, 1, th);
    cyc.rest(10 * 60, 0.1, 1, th);
    cyc.CCCV(16, module.Vmax(), 50e-3, 0.1, 5, th);
    cyc.rest(10 * 60, 0.1, 10, th);
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  cyc.writeData();
}


inline void run_Cell_ECM_SmallPack()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "Cell_ECM_SmallPack"; // + std::to_string(Crate) + '_'

  auto c = Cell_ECM();

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

inline void run_Cell_ECM_MediumPack()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "Cell_ECM_MediumPack"; // + std::to_string(Crate) + '_'

  auto c = Cell_ECM();
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

inline void run_Cell_ECM_LargePack()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "Cell_ECM_LargePack"; // + std::to_string(Crate) + '_'

  auto c = Cell_ECM();
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

inline void run_Cell_ECM_LargePackLong()
{
  constexpr double T_ENV = 15.0_degC;
  std::string ID = "Cell_ECM_LargePackLong"; // + std::to_string(Crate) + '_'

  auto c = Cell_ECM();

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