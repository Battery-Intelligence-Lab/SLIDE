/*
 * Interpolation.cpp
 *
 * Groups functions to do linear interpolation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */


#include "Interpolation.h"
#include <iostream>
#include <math.h>
#include <cassert>

using namespace std;

double linInt(bool verbose, bool bound, double xdat[], double ydat[], int nin, double x){
	/*
	 * function for linear interpolation with the data points provided as two arrays
	 *
	 * IN
	 * verbose 	if false, no error message is written (but the error is still thrown)
	 * 			if true, an error message is written
	 * bound	boolean deciding what to do if the value of x is out of range of xdat
	 * 			if true, an error is thrown
	 * 			if false, the value closest to x is returned
	 * xdat 	x data points in strictly increasing order
	 * ydat 	y data points
	 * nin 		number of data points
	 * x 		x point at which value is needed
	 *
	 * OUT
	 * y 		y value corresponding to x
	 *
	 * THROWS
	 * 1 		if bound = true && if x is out of bounds, i.e. smaller than the smallest value of xdat or larger than the largest value of xdat
	 * 				i.e. bound AND (x < xdat [0] OR x > xdat[end])
	 */

	// check that x is within the limits of the data points
	if ( (x < xdat[0]) || (x > xdat[nin-1]) ){
		if(bound){						// throw an error
			if(verbose)
				cerr<<"ERROR in Interpolation::linInt: x is out of bounds. x = "<<x<<" while xmin = "<<xdat[0]<<" and xmax is "<<xdat[nin-1]<<endl<<flush;
			throw 1;
		}
		else{							// return the value of the closest point
			// if x is below the minimum value, return the y value of the first data point
			if ( x <= xdat[0])
				return ydat[0];

			// if x is above the maximum value, return the y value of the last data point
			if (x >= xdat[nin-1])
				return ydat[nin-1];
		}
	}

	// scan the data points
	double yy = nan("double");					// start with a nan value for y
	for (int i=0;i<nin;i++){					// loop through the data point
		if (x == xdat[i]){						// if the x-value is exactly the x-value of a data point
			yy = ydat[i];						// return the corresponding y-value
			break;								// stop scanning the data points
		}
		else if (xdat[i] > x) {					// the data points xdat are increasing, so find the first data point larger than x
			double xr = xdat[i];				// then that point is the point 'to the right' of x
			double yr = ydat[i];
			double xl = xdat[i-1];				// while the previous point is the point 'to the left' of x
			double yl = ydat[i-1];
			yy = yl + (yr-yl)*(x-xl)/(xr-xl);	// interpolate linearly between the points to the left and right
			break;								// stop searching
		}										// else go to the next point
	}

	return yy;
}


double linInt_matrix(bool verbose, bool bound, double dat[][2], int nin, double x){
	/*
	 * function for linear interpolation with the data points provided as a matrix
	 *
	 * IN
	 * verbose 	if false, no error message is written (but the error is still thrown)
	 * 			if true, an error message is written
	 * bound	boolean deciding what to do if the value of x is out of range of xdat
	 * 			if true, an error is thrown
	 * 			if false, the value closest to x is returned
	 * dat 		matrix with the data points
	 * 			the first column dat[][0] must have the x data points in strictly increasing order
	 *  		the second column dat[][1] must have the y data points
	 * ydat 	y data points
	 * nin 		number of data points
	 * x 		x point at which value is needed
	 *
	 * OUT
	 * y 		y value corresponding to x
	 *
	 * THROWS
	 * 1 		if bound = true && if x is out of bounds, i.e. smaller than the smallest value of xdat or larger than the largest value of xdat
	 * 				i.e. bound AND (x < xdat [0] OR x > xdat[end])
	 */

	// Split the matrix in two arrays
	double xdat[nin], ydat[nin];
	for(int i=0;i<nin;i++){
		xdat[i] = dat[i][0];
		ydat[i] = dat[i][1];
	}

	// Call the function for linear interpolation
	double y;
	try{
		y = linInt(verbose, bound, xdat, ydat, nin, x);
	}
	catch(int e){
		if(verbose)
			cout<<"Error in Interpolation::linInt_matrix when calling linInt: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// make the output variable
	return y;
}


