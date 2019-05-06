/*
 * ReadCSVfiles.h
 *
 * groups functions for reading csv files into arrays and matrices
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_READCSVFILES_H_
#define SRC_READCSVFILES_H_

#include <string>
#include "State.hpp"

using namespace std;

void loadCSV_1col(string name, int nin, double x[]);						// read a csv file with 1 column to one array
void loadCSV_2col(string name, int nin, double x[], double y[]);			// read a csv file with 2 columns to two arrays
void loadCSV_2colMatrix(string name, int nin, double x[][2]);				// read a csv file with 2 columns to a matrix with two columns
void loadCSV_mat(string name, int length, int width, double x[][5*nch]);	// read a csv file with many columns (fewer than 5*nch) to a matrix


#endif /* SRC_READCSVFILES_H_ */
