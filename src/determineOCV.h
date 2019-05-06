/*
 * determineOCV.h
 *
 * Header file for the functions used to find the parameters which will match a measured OCV curve of a cell
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_DETERMINEOCV_H_
#define SRC_DETERMINEOCV_H_

using namespace std;

#include <string>

void loadCSV2(string name, int nin, double x[][2]);
bool validOCV(bool checkRange, int nin, double x[][2]);
void readOCVinput(string namepos, string nameneg, string namecell, int np, int nn, int ncell, double OCVp[][2], double OCVn[][2], double OCVcell[][2]);
double linInt(double dat[][2], int nin, double x);
void discharge(int np, int nn, double OCVp[][2], double OCVn[][2], double cap, double AMp, double AMn, double cmaxp, double cmaxn, double sp, double sn,
		double Vend, int nin, int* nout, double OCV[][2], double OCVanode[][2], double OCVcathode[][2], double fp[], double fn[]);
double calculateError(bool bound, int ncell, int nsim, double OCVcell[][2], double OCVsim[][2]);
void estimateOCVparameters();

#endif /* SRC_DETERMINEOCV_H_ */
