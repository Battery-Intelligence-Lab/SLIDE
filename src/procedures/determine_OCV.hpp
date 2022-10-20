/*
 * determineOCV.hpp
 *
 * Header file for the functions used to find the parameters which will match a measured OCV curve of a cell
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 *
 *  Created on: 19 Dec 2019
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../utility/utility.hpp"

#include <string>
#include <array>

namespace slide {
bool validOCV(bool checkRange, slide::XYdata_vv &data);

void readOCVinput(const std::string &namepos, const std::string &nameneg, const std::string &namecell,
                  slide::XYdata_vv &OCVp, slide::XYdata_vv &OCVn, slide::XYdata_vv &OCVcell);

void discharge(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, double cap, const double AMp, const double AMn,
               const double cmaxp, const double cmaxn, double sp, double sn, double Vend, slide::XYdata_vv &OCV,
               slide::XYdata_vv &OCVanode, slide::XYdata_vv &OCVcathode, double fp[], double fn[]);

auto discharge_noexcept(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, double cap, const double AMp, const double AMn,
                        const double cmaxp, const double cmaxn, double sp, double sn, double Vend, slide::XYdata_vv &OCV,
                        slide::XYdata_vv &OCVanode, slide::XYdata_vv &OCVcathode, double fp[], double fn[]);

double calculateError(bool bound, slide::XYdata_vv &OCVcell, slide::XYdata_vv &OCVsim);

void estimateOCVparameters();

void writeOCVParam(int h, const std::array<double, 4> &par);

void fitAMnAndStartingPoints(int hierarchy, int ap, slide::FixedData<double> AMn_space, slide::FixedData<double> sp_space,
                             slide::FixedData<double> sn_space, double cmaxp, double cmaxn,
                             double *err, std::array<double, 4> &par, std::array<int, 4> &parindex,
                             slide::XYdata_vv &OCVp, slide::XYdata_vv &OCVn, slide::XYdata_vv &OCVcell);

auto hierarchicalOCVfit(int hmax, slide::FixedData<double> AMp_space, slide::FixedData<double> AMn_space, slide::FixedData<double> sp_space,
                        slide::FixedData<double> sn_space, std::string namepos, std::string nameneg, std::string namecell, double cmaxp, double cmaxn);

} // namespace slide