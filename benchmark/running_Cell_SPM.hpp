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

inline void run_Cell_SPM_1(double Crate)
{
  std::string ID = "PyBAMM_1_CC_Crate"; // + std::to_string(Crate) + '_'
  auto c = Cell_SPM();

  const auto Idisch = Crate * c.Cap();

  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;

  std::cout << "Vbefore: " << c.V() << '\n';
  c.setBlockDegAndTherm(true);
  ThroughputData th{};

  auto cyc = Cycler(&c, ID);

  Clock clk;
  cyc.CC(Idisch, 2.7, 3600, 0.5, 2, th);
  std::cout << "Finished " << ID << " in " << clk << ".\n";

  std::cout << "Vafter: " << c.V() << '\n';

  cyc.writeData();
};

} // namespace slide::benchmarks