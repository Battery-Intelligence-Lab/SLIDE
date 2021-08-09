/*
 * Model.h
 *
 * Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <array>

#include "state.hpp"
#include "slide_aux.hpp"

namespace slide
{

	// Define a structure with the matrices of the spatial discretisation of the solid diffusion PDE
	// See the matlab script modelSetup.m
	struct Model
	{
		std::array<double, 4> Input; // array with the input parameters of the Matlab files

		std::array<double, settings::nch> xch; // location of the Chebyshev nodes in the positive domain EXCLUDING centre and surface

		// state space model
		//	dzpos/dt = Ap*zpos + Bp*jp		time derivative of (transformed) concentration at the inner nodes
		//	dzneg/dt = An*zneg + Bn*jn
		//	cp = Cp*zpos + Dp*jp			actual concentration [mol m-3] of the nodes (surface, inner)
		// 	cn = Cn*zpos + Dn*jn

		std::array<double, settings::nch> Ap, An; // only main diagonal is non-zero, so only store those values
		std::array<double, settings::nch> Bp, Bn;

		slide::Matrix<double, settings::nch + 1, settings::nch> Cp, Cn;

		std::array<double, settings::nch + 1> Cc, Dp, Dn; // matrix to get the concentration at the centre node

		slide::Matrix<double, settings::nch, settings::nch> Vp, Vn; // inverse of the eigenvectors for the positive/negative electrode

		slide::Matrix<double, 2 * settings::nch + 3, 2 * settings::nch + 3> Q; // Matrix for Chebyshev integration

		Model();
	};

}