/*
 * determine_characterisation.hpp
 *
 * Header file for the functions used to find the parameters which will match measured CCCV cycles.
 * This are the diffusion constants, rate constants and DC resistance.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "../cells/cells.hpp"
#include "../utility/utility.hpp"
#include "CyclerOld.hpp"
#include "Cycler.hpp"

#include <string>

namespace slide {
bool CCCV_fit(Cell_SPM c1, double Crate, double Ccut, double Tref, double Dp, double Dn, double kp,
              double kn, double R, const struct OCVparam &ocvfit, const struct slide::Model_SPM &M,
              slide::XYdata_vv &Vsim, slide::XYdata_vv &Tsim);

void CCCV(double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, const struct OCVparam &ocvfit,
          const struct slide::Model_SPM &M, slide::XYdata_vv &Vsim, slide::XYdata_vv &Tsim);

void fitDiffusionAndRate(int hierarchy, int ir, double R, slide::FixedData<double> Dp_space, slide::FixedData<double> Dn_space,
                         slide::FixedData<double> kp_space, slide::FixedData<double> kn_space,
                         std::vector<slide::XYdata_vv> &Vdata_all, double weights[],
                         double Crates[], double Ccuts[], double Tref, const struct OCVparam &ocvfit,
                         double *err, std::array<double, 5> &par);

void hierarchicalCharacterisationFit(int hmax, slide::FixedData<double> r_space, slide::FixedData<double> Dp_space,
                                     slide::FixedData<double> Dn_space, slide::FixedData<double> kp_space,
                                     slide::FixedData<double> kn_space, std::vector<slide::XYdata_vv> &Vdata_all,
                                     double weights[], double Crates[], double Ccuts[], double Tref,
                                     const struct OCVparam &ocvfit, double *err, std::array<double, 5> &par);

void writeCharacterisationParam(int h, const std::array<double, 5> &par, double err);

void estimateCharacterisation();
} // namespace slide