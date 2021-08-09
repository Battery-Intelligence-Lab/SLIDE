/*
 * determineOCV.h
 *
 * Header file for the functions used to find the parameters which will match a measured OCV curve of a cell
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <string>
#include <array>

#include "slide_aux.hpp"

bool validOCV(bool checkRange, slide::vec_XYdata &data);

void readOCVinput(const std::string &namepos, const std::string &nameneg, const std::string &namecell,
				  slide::vec_XYdata &OCVp, slide::vec_XYdata &OCVn, slide::vec_XYdata &OCVcell);

void discharge(const slide::vec_XYdata &OCVp, const slide::vec_XYdata &OCVn, double cap, const double AMp, const double AMn,
			   const double cmaxp, const double cmaxn, double sp, double sn, double Vend, slide::vec_XYdata &OCV,
			   slide::vec_XYdata &OCVanode, slide::vec_XYdata &OCVcathode, double fp[], double fn[]);

auto discharge_noexcept(const slide::vec_XYdata &OCVp, const slide::vec_XYdata &OCVn, double cap, const double AMp, const double AMn,
						const double cmaxp, const double cmaxn, double sp, double sn, double Vend, slide::vec_XYdata &OCV,
						slide::vec_XYdata &OCVanode, slide::vec_XYdata &OCVcathode, double fp[], double fn[]);

double calculateError(bool bound, slide::vec_XYdata &OCVcell, slide::vec_XYdata &OCVsim);

void estimateOCVparameters();

void writeOCVParam(int h, const std::array<double, 4> &par);

void fitAMnAndStartingPoints(int hierarchy, int ap, slide::fixed_data<double> AMn_space, slide::fixed_data<double> sp_space,
							 slide::fixed_data<double> sn_space, double cmaxp, double cmaxn,
							 double *err, std::array<double, 4> &par, std::array<int, 4> &parindex,
							 slide::vec_XYdata &OCVp, slide::vec_XYdata &OCVn, slide::vec_XYdata &OCVcell);

auto hierarchicalOCVfit(int hmax, slide::fixed_data<double> AMp_space, slide::fixed_data<double> AMn_space, slide::fixed_data<double> sp_space,
						slide::fixed_data<double> sn_space, std::string namepos, std::string nameneg, std::string namecell, double cmaxp, double cmaxn);