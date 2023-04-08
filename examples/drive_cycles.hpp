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
#include <memory>
#include <fstream>
#include <span>

namespace slide::examples {
inline auto get_exampleCell()
{

  // Make one cell with standard parameters// Degradation settings
  slide::DEG_ID deg;
  // Kinetic SEI + porosity + Dai/Laresgoiti LAM

  deg.SEI_id.add_model(4); // add model 4.
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0); // add model 0.
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  auto example_cell = make<Cell_SPM>("cell_ancillary", deg, 1, 1, 1, 1);

  // Specify the OCV parameters (calculated by determineOCV::estimateOCVparameters)
  OCVparam ocvfit;
  ocvfit.elec_surf = (3.4 / 2.65) * 0.0982; // electrode surface Cap_ratio
  ocvfit.ep = 0.5;                          // volume fraction of active material in the cathode
  ocvfit.en = 0.5;                          // volume fraction of active material in the anode
  ocvfit.thickp = 70e-6;                    // thickness of the cathode
  ocvfit.thickn = 73.5e-6;                  // thickness of the anode
  ocvfit.lifracpini = 0.6862;               // lithium fraction in the cathode at 50% soC
  ocvfit.lifracnini = 0.4843;               // lithium fraction in the anode at 50% SOC
  ocvfit.cmaxp = 51385;                     // maximum lithium concentration in the cathode [mol m-3]
  ocvfit.cmaxn = 30555;                     // maximum lithium concentration in the anode [mol m-3]
  ocvfit.cap = 3.4;                         // the capacity of the cell [Ah]
  ocvfit.Vmax = 4.2;                        // maximum voltage of the cell [V]
  ocvfit.Vmin = 2.6;                        // minimum voltage of the cell [V]

  // c.setOCVcurve(ocvfit.namepos, ocvfit.nameneg);
  example_cell->setInitialConcentration(ocvfit.cmaxp, ocvfit.cmaxn, ocvfit.lifracpini, ocvfit.lifracnini);
  example_cell->setGeometricParameters(ocvfit.cap, ocvfit.elec_surf, ocvfit.ep, ocvfit.en, ocvfit.thickp, ocvfit.thickn);

  example_cell->setT(22.0_degC);    // set the temperature of the cell to the given value
  example_cell->setTenv(22.0_degC); // set the environmental temperature to the given value

  std::cout << "Voltage: " << example_cell->V() << " V.\n";
  std::cout << "Current: " << example_cell->I() << " A.\n";
  std::cout << "Capacity: " << example_cell->Cap() << " Ah.\n";

  return example_cell;
}


auto inline init_cell_manual(auto &su, double I_0, double V_0, double &Ah, double &ttot, double &dWh)
{
  Ah = 0;
  ttot = 0;
  dWh = 0;

  su->setBlockDegAndTherm(true);

  if (std::abs(I_0) < 0.01) {
    const double V_lim = V_0;
    const auto testV = su->V();
    if (V_lim > su->V()) {
      su->setCurrent(-0.02 * su->Cap(), true);
      while (V_lim > su->V()) {
        su->timeStep_CC(0.1, 1);
        Ah += 0.5 * 0.1;
        ttot += 0.1;
        dWh += 0.5 * 0.1 * su->V();
      }
    } else {
      su->setCurrent(0.02 * su->Cap(), true);
      while (V_lim < su->V()) {
        su->timeStep_CC(0.1, 1);
        Ah += 0.5 * 0.1;
        ttot += 0.1;
        dWh += 0.5 * 0.1 * su->V();
      }
    }

    su->setCurrent(I_0, true);
  } else {
    if (I_0 < 0) // This means we are charging. So empty the cell then chart.
    {
      const double V_lim = std::max(su->Vmin(), (V_0 - 0.2));

      if (V_lim < su->V()) // In case it is already empty.
        su->setCurrent(5, true);

      while (V_lim < su->V()) {
        su->timeStep_CC(2, 1);
        Ah += 5 * 2;
        ttot += 2;
        dWh += 5 * 2 * su->V();
      }

      su->setCurrent(I_0, true);

      while (V_0 > su->V()) {
        su->timeStep_CC(1, 1);
        Ah += std::abs(I_0);
        ttot += 1;
        dWh += std::abs(I_0) * su->V();
      }
    } else {
      const double V_lim = std::min(su->Vmax(), (V_0 + 0.2));

      if (V_lim > su->V()) // In case it is already full.
        su->setCurrent(-5, true);

      while (V_lim > su->V()) {
        su->timeStep_CC(2, 1);
        Ah += 5 * 2;
        ttot += 2;
        dWh += 5 * 2 * su->V();
      }

      su->setCurrent(I_0, true);

      while (V_0 < su->V()) {
        su->timeStep_CC(1, 1);
        Ah += std::abs(I_0);
        ttot += 1;
        dWh += std::abs(I_0) * su->V();
      }
    }
  }

  Ah /= 3600.0;
  dWh /= 3600.0;
  su->setBlockDegAndTherm(false);
  Ah = 0;
  ttot = 0;
  dWh = 0;
}


inline void drive_cycle_artemis()
{

  // Note: Entropic effect must be added!

  std::string ID = "temp";
  Clock clk;

  // double Tref = 21.0_degC; // Temperature at which the characterisation should be done [K]
  // Our data is between 23.5 and 25.9 with mean 24.2 C temperature. 26.44 for test data.
  if (settings::T_MODEL != 0) {
    std::cerr << "drive_cycle_artemis works with T_MODEL=0 but it is not!\n";
    throw 1234;
  }

  std::string profile_path{ "profiles/drive_cycles/ArtemisM_scaled.csv" };

  // Make one cell with standard parameters// Degradation settings
  slide::DEG_ID deg;
  // Kinetic SEI + porosity + Dai/Laresgoiti LAM

  deg.SEI_id.add_model(4); // add model 4.
  deg.SEI_porosity = 1;

  deg.CS_id.add_model(0); // add model 0.
  deg.CS_diffusion = 0;

  deg.LAM_id.add_model(1);
  deg.pl_id = 0;

  auto c = Cell_SPM("cell_ancillary", deg, 1, 1, 1, 1);
  c.setBlockDegAndTherm(true);
  c.setT(21.0_degC);

  double Cmaxpos{ 51385 };
  double Cmaxneg{ 30555 };
  double cps, cns;


  auto &st = c.getStateObj();
  auto cyc = Cycler(&c, "discharge");

  ThroughputData th{};
  cyc.CCCV(1, 2.7, 0.0001, 1, 0, th);
  cyc.rest(100, 1, 0, th);

  std::ofstream out_fraction{ PathVar::results / "Kokam_NMC_16Ah_OCV.csv" };

  out_fraction << "Current [A],"
               << "Voltage [V],"
               << "Li+ fraction [-],"
               << "Li- fraction [-],"
               << "zp(3),"
               << "zn(3)\n";

  c.getCSurf(cps, cns, false);
  double V_old{ c.V() }, fp_old{ cps / Cmaxpos }, fn_old{ cns / Cmaxneg };
  out_fraction << c.I() << ',' << c.V() << ',' << fp_old << ',' << fn_old << ',' << st.zp(3) << ',' << st.zn(3) << '\n';

  c.setCurrent(-0.02);
  double threshold{ 1e-4 };
  while (c.V() < 4.2) {
    c.getCSurf(cps, cns, false);
    auto fp{ cps / Cmaxpos }, fn{ cns / Cmaxneg };

    if (std::abs(V_old - c.V()) > threshold || std::abs(fp_old - fp) > threshold || std::abs(fn_old - fn) > threshold) {
      V_old = c.V();
      fp_old = fp;
      fn_old = fn;
      out_fraction << c.I() << ',' << c.V() << ',' << fp << ',' << fn << ',' << st.zp(3) << ',' << st.zn(3) << '\n';
    }
    c.timeStep_CC(1, 1);
  }

  out_fraction.close();

  c.getCSurf(cps, cns, false);
  std::cout << "V: " << c.V() << " cps, cns : " << cps / Cmaxpos << ", " << cns / Cmaxneg << ',' << st.zp(3) << ',' << st.zn(3) << "\n";
  for (auto z_i : st.z())
    std::cout << z_i << ' ';
  std::cout << '\n';


  double SOC_vs_Volt[] = { 2.7, 3.4872, 3.5606, 3.6206, 3.6492, 3.6845, 3.7635, 3.8388, 3.9289, 4.0439, 4.2 };
  double SOC_vs_zn3[] = { 0.017456, 0.072543, 0.127650, 0.182757, 0.237865, 0.292971, 0.348078, 0.403186, 0.458292, 0.513400, 0.568507 };
  double SOC_vs_zp3[] = { 0.670246, 0.630914, 0.591568, 0.552221, 0.512874, 0.473528, 0.434182, 0.394835, 0.355489, 0.316142, 0.276795 };

  double SOC_vs_fp[] = { 0.970530, 0.913563, 0.856588, 0.799614, 0.742639, 0.685665, 0.628691, 0.571716, 0.514742, 0.457767, 0.400792 };
  double SOC_vs_fn[] = { 0.028906, 0.120170, 0.211423, 0.302676, 0.393928, 0.485181, 0.576434, 0.667687, 0.758939, 0.850192, 0.941444 };

  for (int i = 0; i < 11; i++) {
    c.setC(SOC_vs_fp[i], SOC_vs_fn[i]);
    std::cout << "i = " << i << " V = " << c.V() << '\n';
  }

  DynamicMatrix<double> profile;
  loadCSV_Ncol(PathVar::data / profile_path, profile);

  // Set the characterisation parameters of the cell to the ones given as input to this function
  auto experiment = [&](std::string name) {
    std::vector<double> voltage;
    std::vector<State_SPM> states;
    std::vector<double> estimated_SOC;
    double cap = c.Cap();
    for (auto c_rate : profile.data) {
      c.setCurrent(c_rate * cap, true);
      voltage.push_back(c.V());
      states.push_back(c.getStateObj());

      estimated_SOC.push_back(100.0 * ((st.zn(3) - SOC_vs_zn3[0]) / (SOC_vs_zn3[10] - SOC_vs_zn3[0])));

      c.timeStep_CC(1, 1);

      if (st.zn(3) < SOC_vs_zn3[1]) // Until 10%
        break;
    }

    std::string Vname = std::string("drive_voltage_from_") + name + "_percent.csv";
    std::ofstream out(PathVar::results / Vname, std::ios_base::out);


    out << "Current [A],Voltage[V],Estimated SOC[%]\n";

    for (size_t i = 0; i < voltage.size(); i++)
      out << profile.data[i] * cap << ',' << voltage[i] << ',' << estimated_SOC[i] << '\n';

    out.close();


    std::string Sname = std::string("drive_states_from_") + name + "_percent.csv";
    std::ofstream Sout(PathVar::results / Sname, std::ios_base::out);

    Sout << "I[A],V[V],T[T],delta,LLI,thickp,thickn,ep,en,ap,an,CS,Dp,Dn,delta_pl,"
         << "zp_0,zp_1,zp_2,zp_3,zp_4,zn_0,zn_1,zn_2,zn_3,zn_4,rDCp,rDCn,rDCcc\n";

    for (auto st_i : states) {
      for (auto s : std::span<double>(st_i.begin(), st_i.begin() + 28))
        Sout << s << ',';

      Sout << '\n';
    }

    Sout.close();
  };

  // Experiment from 90%
  // c.setCurrent(0.01);
  // while (c.V() > SOC_vs_Volt[9]) // 9 -> 90%
  //   c.timeStep_CC(1, 1);

  for (int i = 9; i > 1; i--) {
    c.setC(SOC_vs_fp[i], SOC_vs_fn[i]);
    experiment(std::to_string(i) + "0");
  }

  // c.writeData("drive_cycle");
  std::cout << "V: " << c.V() << '\n';
  std::cout << "SOC: " << 100 * c.SOC() << '\n';
  std::cout << "Finished drive_cycle example in " << clk << ".\n";
}

} // namespace slide::examples