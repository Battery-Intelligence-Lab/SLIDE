/*
 * Model.cpp
 *
 * Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Model.h"
#include "ReadCSVfiles.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <assert.h>
#include <memory>
#include <thread>

void Model_initialise(Model& M){
	/*
	 * Function to initialise the matrices of the spatial discretisation of the solid diffusion PDE.
	 * The matrices are calculated by the Matlab function modelSetup.m and written to csv files.
	 * This function reads those csv files.
	 *
	 * IN
	 * M 	Model structure (defined in Cell.hpp) with the matrices for the model
	 */

	try{

		// Read the data which have one column (arrays)
		loadCSV_1col("Cheb_Nodes.csv", nch, M.xch);
		loadCSV_1col("Cheb_An.csv", nch, M.An);
		loadCSV_1col("Cheb_Ap.csv", nch, M.Ap);
		loadCSV_1col("Cheb_Bn.csv", nch, M.Bn);
		loadCSV_1col("Cheb_Bp.csv", nch, M.Bp);
		loadCSV_1col("Cheb_Cc.csv", nch+1, M.Cc);
		loadCSV_1col("Cheb_Dn.csv", nch, M.Dn);
		loadCSV_1col("Cheb_Dp.csv", nch, M.Dp);
		loadCSV_1col("Cheb_input.csv", 4, M.Input);

		// Read the data which have multiple columns (matrices)
		// In C++ you can't pass variable-size matrices, so we make a matrix A which is larger than what we need.
		// We then copy the rows and columns we need to the actual matrices.
		double A[5*nch][5*nch];
		int n,m;

		// M.Cn
		n = nch+1;								// number of rows
		m = nch;								// number of columns
		loadCSV_mat("Cheb_Cn.csv", n, m, A);	// read data into A
		for(int i=0;i<n;i++){					// copy the first n rows
			for(int j=0;j<m;j++)				// copy the first m columns
				M.Cn[i][j] = A[i][j];
		}

		// M.Cp
		n = nch+1;								// number of rows
		m = nch;								// number of columns
		loadCSV_mat("Cheb_Cp.csv", n, m, A);	// read data into A
		for(int i=0;i<n;i++){					// copy the first n rows
			for(int j=0;j<m;j++)				// copy the first m columns
				M.Cp[i][j] = A[i][j];
		}

		// M.Vn
		n = nch;								// number of rows
		m = nch;								// number of columns
		loadCSV_mat("Cheb_Vn.csv", n, m, A);	// read data into A
		for(int i=0;i<n;i++){					// copy the first n rows
			for(int j=0;j<m;j++)				// copy the first m columns
				M.Vn[i][j] = A[i][j];
		}

		// M.Vp
		n = nch;								// number of rows
		m = nch;								// number of columns
		loadCSV_mat("Cheb_Vp.csv", n, m, A);	// read data into A
		for(int i=0;i<n;i++){					// copy the first n rows
			for(int j=0;j<m;j++)				// copy the first m columns
				M.Vp[i][j] = A[i][j];
		}

		// M.Q
		n = 2*nch+3;							// number of rows
		m = 2*nch+3;							// number of columns
		loadCSV_mat("Cheb_Q.csv", n, m, A);		// read data into A
		for(int i=0;i<n;i++){					// copy the first n rows
			for(int j=0;j<m;j++)				// copy the first m columns
				M.Q[i][j] = A[i][j];
		}

	}
	catch(int e){
		cout<<"Error in Model::Model_initialise when reading the files: "<<e<<". Throwing it on"<<endl;
		throw e;
	}
}


