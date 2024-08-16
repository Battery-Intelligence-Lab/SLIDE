/**
 * @file module_tests.hpp
 * @brief Example module test functions
 * @author Volkan Kumtepeli
 * @date 01 Aug 2024
 */

#pragma once

#include "../src/slide.hpp"

#include <vector>

namespace slide::examples {

void module_p_ECM()
{
  using namespace slide;
  // Make a module with N cells and a contact resistance
  constexpr int Ncells = 3;
  double Rc = 0.01;

  using cell_type = Cell_ECM<1>;

  std::vector<Module_p::SU_t> cs;
  std::vector<double> Rcs{};

  for (int i{}; i < Ncells; i++) {
    cs.push_back(make<cell_type>("cell" + std::to_string(i)));
    Rcs.push_back(Rc);
  }

  std::string n = "parECM";

  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = make<Module_p>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  mp->setRcontact(Rcs);

  //!< total resistance:
  //!< 		-Rp+-Rp-+-Rp-|
  //!< 		   |    |    |
  //!< 		  Rs    Rs   Rs
  //!< where Rp = contact resistance (value Rc) and Rs = cell resistance = 0.01;

  Cycler cyc(mp, "Cycler1");

  ThroughputData th{};
  double dt = 0.01;
  std::cout << "Voltage: " << mp->V() << " I: " << mp->I() << " A.\n";


  cyc.CC(-16, 4.2, 3600, dt, 1, th);

  std::cout << "Voltage: " << mp->V() << " I: " << mp->I() << " A.\n";

  cyc.writeData();
}

} // namespace slide::examples