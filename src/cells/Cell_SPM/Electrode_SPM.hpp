/**
 * @file Electrode_SPM.hpp
 * @brief Electrode class file, currently only for SPM.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 11 Aug 2024
 */

#pragma once

#include "constants.hpp"
#include "State_SPM.hpp"
#include "enum_definitions.hpp"
#include "Pair.hpp"

#include <cmath>

namespace slide {

struct Electrode_SPM
{
  Domain dom{ 10 };    //!< Just to give an error in default state.
  double k{ 5e-11 };   //!< rate constant of main reaction at reference temperature
  double k_T{ 58000 }; //!< activation energy for the Arrhenius relation of k
  double D_T{};
  double Cmax{};
  double x_0{}, x_100{};

  //!< Calculate the rate constants at the cell's temperature using an Arrhenius relation
  auto kt(double Arrhenius) const noexcept { return k * std::exp(k_T * Arrhenius); }
  auto i0(double cs, double Arrhenius) //!< exchange current density
  {
    using namespace PhyConst;
    return kt(Arrhenius) * n * F * std::sqrt(C_elec * cs * (Cmax - cs));
  }

  auto x(double cs, double Arrhenius, double i_app, State_SPM &st)
  {
    const auto i0_now = i0(cs, Arrhenius);
    return 0.5 * sign(dom) * i_app / (st.a(dom) * st.thick(dom)) / i0_now;
  }

  auto eta(double cs, double Arrhenius, double i_app, State_SPM &st)
  {
    //!< Calculate the overpotential using the Butler-Volmer equation
    //!< if alpha is 0.5, the Butler-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
    //!< and asinh(x) = ln(x + sqrt(1+x^2) -> to asinh(x) function.
    using namespace PhyConst;
    const auto x_ = x(cs, Arrhenius, i_app, st);
    return (2 * Rg * st.T()) / (n * F) * std::asinh(x_);
  }

  auto overpotential(double cs, double Arrhenius, double i_app, State_SPM &st) { return eta(cs, Arrhenius, i_app, st); }

  auto molarFlux(double i_app, State_SPM &st)
  { //!< molar flux on the particle [mol m-2 s-1]
    using namespace PhyConst;
    return sign(dom) * i_app / (st.a(dom) * n * F * st.thick(dom));
  }

  auto Dt(double Arrhenius, State_SPM &st)
  {
    return st.D(dom) * std::exp(D_T * Arrhenius); //!< diffusion constant of the positive particle [m s-1]
  }

  bool check_c_surf(double cs) { return cs > 0 && cs < Cmax; } //!< check if the surface concentration is within the allowed range
  bool check_c_surf(DPair cs) { return cs[dom] > 0 && cs[dom] < Cmax; }

  double z_surf(double cs) { return cs / Cmax; }
  double z_surf(DPair cs) { return cs[dom] / Cmax; }
};


} // namespace slide