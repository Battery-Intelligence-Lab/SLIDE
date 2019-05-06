/*
 * Interpolation.h
 *
 * Groups functions to do linear interpolation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_INTERPOLATION_H_
#define SRC_INTERPOLATION_H_


double linInt(bool verbose, bool bound, double xdat[], double ydat[], int nin, double x); 	// interpolate between data from two arrays
double linInt_matrix(bool verbose, bool bound, double dat[][2], int nin, double x);			// interpolate between data in a matrix with two columns


#endif /* SRC_INTERPOLATION_H_ */
