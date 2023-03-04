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

inline void run_Cell_Bucket_single_default_pulse()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_Bucket_single_default_pulse"; // + std::to_string(Crate) + '_'
  auto c = Cell_Bucket();
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

inline void run_Cell_Bucket_single_default_CCCV()
{
  // Benchmark with default parameters:
  std::string ID = "Cell_Bucket_single_default_CCCV"; // + std::to_string(Crate) + '_'
  auto c = Cell_Bucket();
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


inline void run_Cell_Bucket()
{
  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 2;

  auto c = Cell_Bucket();
  auto cyc = Cycler(&c, ID);

  for (size_t i{ 0 }; i < N; i++) {
    ThroughputData th;

    cyc.CCCV(1, 4, 0.1, 1, 1, th);
    cyc.CCCV(1, 3, 0.1, 1, 1, th);
  }

  c.writeData("test");
  std::cout << "V: " << c.V() << "\n";

  std::cout << "Finished run_Cell_Bucket in " << clk << ".\n";
}

inline void run_Cell_ECM()
{
  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 2;

  auto c = Cell_ECM();
  auto cyc = Cycler(&c, ID);

  for (size_t i{ 0 }; i < N; i++) {
    ThroughputData th;
    cyc.CCCV(1, 4, 0.1, 1, 1, th);
    cyc.CCCV(1, 3, 0.1, 1, 1, th);
  }

  c.writeData("test");
  std::cout << "V: " << c.V() << "\n";

  std::cout << "Finished run_Cell_ECM in " << clk << ".\n";
}

inline void run_Cell_SPM()
{

  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 15;

  auto c = Cell_SPM();
  std::cout << "V: " << c.V() << '\n';
  // auto cyc = Cycler(&c, ID);

  // double Ah, Wh, dtime;
  // for (size_t i{ 0 }; i < N; i++) {


  //   cyc.CCCV(1, 4, 0.1, 1, 10, Ah, Wh, dtime);
  //   cyc.CCCV(1, 3, 0.1, 1, 10, Ah, Wh, dtime);
  // }

  // c.writeData("test");
  // std::cout << "V: " << c.V() << '\n';
  // std::cout << "Wh: " << Wh << '\n';
  // std::cout << "Finished run_Cell_SPM in " << clk << ".\n";
}

} // namespace slide::benchmarks