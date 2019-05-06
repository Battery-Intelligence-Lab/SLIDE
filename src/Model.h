/*
 * Model.h
 *
 * Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_MODEL_H_
#define SRC_MODEL_H_

#include "State.hpp"

// Define a structure with the matrices of the spatial discretisation of the solid diffusion PDE
// See the matlab script modelSetup.m
struct Model{
	double Input[4];		// array with the input parameters of the Matlab files

	double xch[nch]; 		// location of the chebyshev nodes in the positive domain EXCLUDING centre and surface

	// state space model
		//	dzpos/dt = Ap*zpos + Bp*jp		time derivative of (transformed) concentration at the inner nodes
		//	dzneg/dt = An*zneg + Bn*jn
		//	cp = Cp*zpos + Dp*jp			actual concentration [mol m-3] of the nodes (surface, inner)
		// 	cn = Cn*zpos + Dn*jn
	double Ap[nch]; 		// only main diagonal is non-zero, so only store those values
	double Bp[nch];
	double An[nch]; 		// only main diagonal is non-zero, so only store those values
	double Bn[nch];
	double Cp[nch+1][nch];
	double Cn[nch+1][nch];
	double Cc[nch+1];		// matrix to get the concentration at the centre node
	double Dp[nch];
	double Dn[nch];
	double Vn[nch][nch];	// inverse of the eigenvectors for the negative electrode
	double Vp[nch][nch];	// inverse of the eigenvectors for the positive electrode

	double Q[2*nch+3][2*nch+3];	// Matrix for Chebyshev integration
};


void Model_initialise(Model& M); 	// initialise the values of the matrices



#endif /* SRC_MODEL_H_ */
