/**
 * @file benchmark_Cell_SPM.cpp
 * @brief Benchmark file for Cell_SPM
 * @author Volkan Kumtepeli
 * @date 07 Aug 2022
 */

#include "slide.hpp"

#include <string>
#include <fstream>
#include <iostream>

using namespace slide;

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
}


inline void run_Cell_SPM_2(double Crate)
{
  std::string ID = "PyBAMM_2_CC_Crate_Cell_SPM_cellData"; // + std::to_string(Crate) + '_'
  auto c = Cell_SPM();

  const auto Idisch = 1.0; // Crate * c.Cap();

  c.getStateObj().rDCn() = 0;
  c.getStateObj().rDCp() = 0;
  c.getStateObj().rDCcc() = 0;

  std::cout << "Vbefore: " << c.V() << '\n';
  c.setBlockDegAndTherm(true);
  // ThroughputData th{};

  c.setCurrent(Idisch);

  const double T_end{ 3600 }, dt{ 0.5 };
  double t_now{ 0 };
  Clock clk;

  std::ofstream file{ PathVar::results / (ID + ".csv"), std::ios::out };

  file << "\n\n\n";
  while (t_now <= T_end && c.V() > 2.7) {
    file << t_now << ',' << c.V() << ',' << c.SOC() << ',' << c.T() << ','
         << c.getThroughputs().time() << ',' << c.getThroughputs().Ah() << ','
         << c.getThroughputs().Wh() << ',' << c.getOCV() << "\n";

    c.timeStep_CC(dt, 2);

    t_now += dt * 2;
  }

  std::cout << "Finished " << ID << " in " << clk << ".\n";
  std::cout << "Vafter: " << c.V() << '\n';

  file.close();
}

int main()
{
  double Crate = 0.5;

  run_Cell_SPM_1(Crate);
  run_Cell_SPM_2(Crate);

  return EXIT_SUCCESS;
}