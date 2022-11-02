/*
 * drive_cycles.hpp
 *
 *  Example drive cycle applications
 *
 *  Created on: 02 Nov 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>

namespace slide::examples {

inline void drive_cycle_artemis()
{
  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 5;

  auto c = Cell_SPM();
  auto cyc = Cycler(&c, ID);

  std::cout << "size SPM: " << sizeof(c) << '\n';

  c.setBlockDegAndTherm(true);


  double Ah, Wh, dtime;
  for (size_t i{ 0 }; i < N; i++) {
    cyc.CCCV(1, 4, 0.1, 1, 10, Ah, Wh, dtime);
    cyc.CCCV(1, 3, 0.1, 1, 10, Ah, Wh, dtime);
  }

  c.writeData("drive_cycle");
  std::cout << "V: " << c.V() << '\n';
  std::cout << "Wh: " << Wh << '\n';
  std::cout << "Finished drive_cycle example in " << clk << ".\n";
};

} // namespace slide::examples