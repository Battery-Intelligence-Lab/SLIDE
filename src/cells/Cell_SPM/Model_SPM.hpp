/*
 * Model_SPM.hpp
 *
 * Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "../../types/matrix.hpp"
#include "../../settings/settings.hpp"
#include "../../utility/utility.hpp"

#include <array>

namespace slide {
//!< Define a structure with the matrices of the spatial discretisation of the solid diffusion PDE
//!< See the matlab script modelSetup.m
//!< This class has 326 elements. So it should be allocated in heap once.
struct Model_SPM
{
  constexpr static auto nch = settings::nch;
  std::array<double, 4> Input; //!< array with the input parameters of the MATLAB files

  std::array<double, nch> xch; //!< location of the Chebyshev nodes in the positive domain EXCLUDING centre and surface

  //!< state space model
  //	dzpos/dt = Ap*zpos + Bp*jp		time derivative of (transformed) concentration at the inner nodes
  //	dzneg/dt = An*zneg + Bn*jn
  //	cp = Cp*zpos + Dp*jp			actual concentration [mol m-3] of the nodes (surface, inner)
  //!< 	cn = Cn*zpos + Dn*jn

  std::array<double, nch> Ap, An; //!< only main diagonal is non-zero, so only store those values
  std::array<double, nch> Bp, Bn;

  slide::Matrix<double, nch + 1, nch> Cp, Cn;

  std::array<double, nch + 1> Cc, Dp, Dn; //!< matrix to get the concentration at the centre node

  slide::Matrix<double, nch, nch> Vp, Vn; //!< inverse of the eigenvectors for the positive/negative electrode

  slide::Matrix<double, 2 * nch + 3, 2 * nch + 3> Q; //!< Matrix for Chebyshev integration

  Model_SPM() //!< #TODO do not do file reading at the constructor.
  {
    /*
     * Constructor to initialise the matrices of the spatial discretisation of the solid diffusion PDE.
     * The matrices are calculated by the MATLAB function modelSetup.m and written to csv files.
     * This function reads those csv files.
     */
    try {
      //!< Read the data which have one column (arrays)
      loadCSV_1col(PathVar::data / "Cheb_Nodes.csv", xch);
      loadCSV_1col(PathVar::data / "Cheb_An.csv", An);
      loadCSV_1col(PathVar::data / "Cheb_Ap.csv", Ap);
      loadCSV_1col(PathVar::data / "Cheb_Bn.csv", Bn);
      loadCSV_1col(PathVar::data / "Cheb_Bp.csv", Bp);
      loadCSV_1col(PathVar::data / "Cheb_Cc.csv", Cc);
      loadCSV_1col(PathVar::data / "Cheb_Dn.csv", Dn);
      loadCSV_1col(PathVar::data / "Cheb_Dp.csv", Dp);
      loadCSV_1col(PathVar::data / "Cheb_input.csv", Input);

      //!< Read the data which have multiple columns (matrices)
      loadCSV_mat(PathVar::data / "Cheb_Cn.csv", Cn);
      loadCSV_mat(PathVar::data / "Cheb_Cp.csv", Cp);
      loadCSV_mat(PathVar::data / "Cheb_Vn.csv", Vn);
      loadCSV_mat(PathVar::data / "Cheb_Vp.csv", Vp);
      loadCSV_mat(PathVar::data / "Cheb_Q.csv", Q);
    } catch (int e) {
      std::cout << "Error in Model_SPM::Model_SPM() when reading the files: "
                << e << ".\n";
      throw e;
    }
  }

  static Model_SPM *makeModel() //!< #TODO make other type of models possible.
  {
    static Model_SPM M;
    return &M;
  }
};

} // namespace slide