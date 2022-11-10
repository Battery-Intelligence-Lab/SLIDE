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

  auto example_cell = std::make_unique<Cell_SPM>("cell_ancillary", deg, 1, 1, 1, 1);

  const double Cap_actual = 100;
  const double Cap_bat = example_cell->Cap();
  const double Cap_ratio = Cap_bat / Cap_actual;

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

  example_cell->setT(C_to_Kelvin(22));    // set the temperature of the cell to the given value
  example_cell->setTenv(C_to_Kelvin(22)); // set the environmental temperature to the given value

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
  std::string ID = "temp";
  Clock clk;


  double Tref = PhyConst::Kelvin + 20; // Temperature at which the characterisation should be done [K]
                                       // Pivot power data is between 23.5 and 25.9 with mean 24.2 C temperature. 26.44 for test data.


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

  auto c = std::make_unique<Cell_SPM>("cell_ancillary", deg, 1, 1, 1, 1);
  c->setBlockDegAndTherm(true);

  DynamicMatrix<double> profile;
  loadCSV_Ncol(PathVar::data + profile_path, profile);

  // Set the characterisation parameters of the cell to the ones given as input to this function


  // auto cyc = Cycler(c.get(), ID);
  double dAh{ 0 }, dtime{ 0 }, dWh{ 0 };

  std::cout << "V: " << c->V() << " [V] SOC: " << 100 * c->SOC() << " [%] \n";
  init_cell_manual(c, 0.0, 4.08, dAh, dtime, dWh); // 4.08 -> 90% SOC and  3.42 V  -> 10% SOC
  std::cout << "V: " << c->V() << " [V] SOC: " << 100 * c->SOC() << " [%] \n";
  init_cell_manual(c, 0.0, 4.2, dAh, dtime, dWh); // 4.08 -> 90% SOC and  3.42 V  -> 10% SOC
  std::cout << "V: " << c->V() << " [V] SOC: " << 100 * c->SOC() << " [%] \n";

  std::cout << "SOC: " << 100 * c->SOC() << '\n';
  std::vector<double> voltage;
  std::vector<State_SPM> states;
  double cap = c->Cap();
  for (auto c_rate : profile.data) {
    c->setCurrent(c_rate * cap, true);
    voltage.push_back(c->V());
    states.push_back(c->getStateObj());
    c->timeStep_CC(1, 1);
  }

  std::string Vname = "drive_voltage.csv";
  std::ofstream out(PathVar::results + Vname, std::ios_base::out);

  for (auto v : voltage)
    out << v << '\n';

  out.close();


  std::string Sname = "drive_states.csv";
  std::ofstream Sout(PathVar::results + Sname, std::ios_base::out);

  for (auto st : states) {
    for (auto s : st)
      Sout << s << ',';

    Sout << '\n';
  }

  Sout.close();


  // c->writeData("drive_cycle");
  std::cout << "V: " << c->V() << '\n';
  std::cout << "Wh: " << dWh << '\n';
  std::cout << "SOC: " << 100 * c->SOC() << '\n';
  std::cout << "Finished drive_cycle example in " << clk << ".\n";
};

} // namespace slide::examples