/*
 * Model.cpp
 *
 * Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include <iostream>

#include "model.h"
#include "read_CSVfiles.h"
#include "constants.hpp"

namespace slide
{

	Model::Model()
	{
		/*
	 * Constructor to initialise the matrices of the spatial discretisation of the solid diffusion PDE.
	 * The matrices are calculated by the Matlab function modelSetup.m and written to csv files.
	 * This function reads those csv files.
	 */
		try
		{
			// Read the data which have one column (arrays)

			loadCSV_1col(PathVar::data + "Cheb_Nodes.csv", xch);
			loadCSV_1col(PathVar::data + "Cheb_An.csv", An);
			loadCSV_1col(PathVar::data + "Cheb_Ap.csv", Ap);
			loadCSV_1col(PathVar::data + "Cheb_Bn.csv", Bn);
			loadCSV_1col(PathVar::data + "Cheb_Bp.csv", Bp);
			loadCSV_1col(PathVar::data + "Cheb_Cc.csv", Cc);
			loadCSV_1col(PathVar::data + "Cheb_Dn.csv", Dn);
			loadCSV_1col(PathVar::data + "Cheb_Dp.csv", Dp);
			loadCSV_1col(PathVar::data + "Cheb_input.csv", Input);

			// Read the data which have multiple columns (matrices)
			loadCSV_mat(PathVar::data + "Cheb_Cn.csv", Cn);
			loadCSV_mat(PathVar::data + "Cheb_Cp.csv", Cp);
			loadCSV_mat(PathVar::data + "Cheb_Vn.csv", Vn);
			loadCSV_mat(PathVar::data + "Cheb_Vp.csv", Vp);
			loadCSV_mat(PathVar::data + "Cheb_Q.csv", Q);
		}
		catch (int e)
		{
			std::cout << "Error in Model::Model() when reading the files: " << e << ". Throwing it on.\n";
			throw e;
		}
	}
}
