/**
 * @file Cell_SPM_diffusion.cpp
 * @brief Diffusion functions for Cell_SPM class
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 2021
 */

#include "Cell_SPM.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>
#include <range/v3/view/iota.hpp>

namespace slide {
double Cell_SPM::calcSurfaceConcentration(double molarFlux, double Dt, Domain dom) //!< Should not throw normally, except divide by zero?
{
  //!< Calculate the surface concentration at the positive particle
  double c_surf{ 0 };
  //!< 	cp_surf = M->Cp[0][:] * zp[:] + M->Dp*jp/Dpt
  for (size_t j{}; j < settings::nch; j++)
    c_surf += M->C[dom](0, j) * st.z(j, dom);

  c_surf += M->D[dom](0) * molarFlux / Dt;

  return c_surf;
}

std::pair<double, DPair> Cell_SPM::calcMolarFlux() //!< Should not throw normally, except divide by zero?
{
  const double i_app = I() / geo.elec_surf; //!< current density on the electrode [A m-2]
  DPair molarFlux;
  for (auto dom : { pos, neg })
    molarFlux[dom] = electrode[dom].molarFlux(i_app, st); //!< molar flux on the particle [mol m-2 s-1]

  return { i_app, molarFlux };
}

} // namespace slide