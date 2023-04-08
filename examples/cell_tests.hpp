/*
 * cell_tests.hpp
 *
 * Example cell test functions;
 *
 * Created on: 04 Apr 2022
 * Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>
#include <memory>
#include <fstream>
#include <span>

namespace slide::examples {

inline auto GITT_test()
{
  // Note: Entropic effect must be added!
  std::string ID = "temp";
  Clock clk;

  // double Tref = 21.0_degC; // Temperature at which the characterisation should be done [K]
  // Our data is between 23.5 and 25.9 with mean 24.2 C temperature. 26.44 for test data.
  if (settings::T_MODEL != 0) {
    std::cerr << "GITT_test example works with T_MODEL=0 but it is not!\n";
    throw 1234;
  }
  slide::DEG_ID deg{};

  auto c = Cell_SPM("cell_ancillary", deg, 1, 1, 1, 1);
  c.setBlockDegAndTherm(true);
  c.setT(21.0_degC);

  double Cmaxpos{ 51385 };
  double Cmaxneg{ 30555 };
  double cps{}, cns{};

  auto &st = c.getStateObj();
  auto cyc = Cycler(&c, "charge");

  auto d = c; // Copy cell for discharge.
  auto dcyc = Cycler(&d, "discharge");

  // Make cell empty!
  ThroughputData th{};
  cyc.CCCV(1, 2.7, 0.0001, 1, 0, th);
  cyc.rest(100, 1, 0, th);

  dcyc.CCCV(1, 4.2, 0.0001, 1, 0, th);
  dcyc.rest(100, 1, 0, th);


  // Start GITT test 20x0.05C pulse and 2 hr rest:
  const auto N_repeat{ 20 };     // Repeat 20 times.
  const auto t_pulse = 1 * 3600; // 1 hr pulse time.
  const auto t_rest = 2 * 3600;  // 2 hr rest time.
  const auto dt = 1;             // 1 seconds time step.
  const auto Crate = 0.05;
  auto current = Crate * c.Cap();

  std::ofstream out_GITT{ PathVar::results / "GITT_20x0.05C_1h_rest_2h.csv" };
  out_GITT << "Time [s],"
           << "Current [A],"
           << "Terminal voltage [V],"
           << "Current [A],"
           << "Terminal voltage [V]\n";

  double t_all{};
  out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
           << d.I() << ',' << d.V() << '\n';

  for (int i{}; i < N_repeat; i++) {
    c.setCurrent(-current);
    d.setCurrent(current);

    for (int j{}; j < 3600; j++) {
      out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
               << d.I() << ',' << d.V() << '\n';

      c.timeStep_CC(1, 1);
      d.timeStep_CC(1, 1);

      t_all += 1;
    }

    c.setCurrent(0);
    d.setCurrent(0);

    for (int j{}; j < 7200; j++) {
      out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
               << d.I() << ',' << d.V() << '\n';

      c.timeStep_CC(1, 1);
      d.timeStep_CC(1, 1);

      t_all += 1;
    }
  }
  out_GITT.close();
}

} // namespace slide::examples