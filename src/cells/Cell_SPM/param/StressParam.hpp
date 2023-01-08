/*
 * StressParam.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

namespace slide::param {
struct StressParam
{
  //!< Constants for the stress model
  double omegap; //!< partial molar volume of positive electrode [m3 mol-1]
  double omegan; //!< partial molar volume of negative electrode [m3 mol-1]
  double Ep;     //!< Young's modulus of positive electrode [GPa]
  double En;     //!< Young's modulus of negative electrode [GPa]
  double nup;    //!< Poisson's ratio of positive electrode [-]
  double nun;    //!< Poisson's ratio of negative electrode [-]

  //!< values of the stress are often needed. Because it takes very long to calculate them, we calculate them once and store them so we don't need to repeat the same calculation twice
  bool s_dai{ false };          //!< do we need to calculate the stress according to Dai's model?
  bool s_lares{ false };        //!< do we need to calculate the stress according to Laresgoiti's model?
  bool s_dai_update{ false };   //!< boolean to indicate if Dai's stress are up to date with the battery state at this time step
  bool s_lares_update{ false }; //!< boolean to indicate if Dai's stress are up to date with the battery state at this time step
  double s_dai_p{};             //!< maximum hydrostatic stress in the positive particle according to Dai's stress model
  double s_dai_n{};             //!< maximum hydrostatic stress in the negative particle according to Dai's stress model
  double s_lares_n{};           //!< stress in the negative particle according to Laresgoiti's stress model
  double s_dai_p_prev{};        //!< maximum hydrostatic stress in the previous time step in the positive particle according to Dai's stress model
  double s_dai_n_prev{};        //!< maximum hydrostatic stress in the previous time step in the negative particle according to Dai's stress model
  double s_lares_n_prev{};      //!< stress in the previous time step in the negative particle according to Laresgoiti's stress model
  double s_dt{ 1 };             //!< time period between the 'previous' and 'current' stress [s]
};
} // namespace slide::param