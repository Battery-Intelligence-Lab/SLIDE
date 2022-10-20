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
inline void run_Cell_Bucket()
{
  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 2;

  auto c = Cell_Bucket();
  auto cyc = Cycler(&c, ID);

  for (size_t i{ 0 }; i < N; i++) {
    double Ah, Wh, dtime;

    cyc.CCCV(1, 4, 0.1, 1, 1, Ah, Wh, dtime);
    cyc.CCCV(1, 3, 0.1, 1, 1, Ah, Wh, dtime);
  }

  c.writeData("test");
  std::cout << "V: " << c.V() << "\n";

  std::cout << "Finished run_Cell_Bucket in " << clk << ".\n";
};

inline void run_Cell_ECM()
{
  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 2;

  auto c = Cell_ECM();
  auto cyc = Cycler(&c, ID);

  for (size_t i{ 0 }; i < N; i++) {
    double Ah, Wh, dtime;

    cyc.CCCV(1, 4, 0.1, 1, 1, Ah, Wh, dtime);
    cyc.CCCV(1, 3, 0.1, 1, 1, Ah, Wh, dtime);
  }

  c.writeData("test");
  std::cout << "V: " << c.V() << "\n";

  std::cout << "Finished run_Cell_ECM in " << clk << ".\n";
};

inline void run_Cell_SPM()
{

  std::string ID = "temp";
  Clock clk;

  constexpr size_t N = 5;

  auto c = Cell_SPM();
  auto cyc = Cycler(&c, ID);

  double Ah, Wh, dtime;
  for (size_t i{ 0 }; i < N; i++) {


    cyc.CCCV(1, 4, 0.1, 1, 10, Ah, Wh, dtime);
    cyc.CCCV(1, 3, 0.1, 1, 10, Ah, Wh, dtime);
  }

  c.writeData("test");
  std::cout << "V: " << c.V() << '\n';
  std::cout << "Wh: " << Wh << '\n';
  std::cout << "Finished run_Cell_SPM in " << clk << ".\n";
};

} // namespace slide::benchmarks