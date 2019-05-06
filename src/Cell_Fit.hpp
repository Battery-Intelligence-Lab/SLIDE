/*
 * CellFit.hpp
 *
 * One of the child classes of Cell.
 *
 * This cell is used to fit the cycling parameters of a cell to data.
 * It is done by the functions in determineCharacterisation.cpp
 * This cell should never be used by other functions (e.g. it should never be used for degradation simulations)
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_CELL_FIT_HPP_
#define SRC_CELL_FIT_HPP_


#include "Cell.hpp"

class Cell_Fit : public Cell{
public:
	Cell_Fit(const struct Model&, int verbosei);								// standard constructor
	virtual ~Cell_Fit();

	void setOCVcurve(string namepos, string nameneg, int np, int nn);			// sets the OCV curve of the cell to the given value
	void setInitialConcentration(double cmaxp, double cmaxn, double lifracp, double lifracn); // sets the initial concentration
	void setGeometricParameters(double capnom, double elec_surf, double ep, double en, double thickp, double thickn); // sets the geometric parameters related to the amount of active material
	void setRamping(double Istep, double tstep);								// sets the ramping parameters

	void setCharacterisationParam(double Dp, double Dn, double kp, double kn, double Rdc); // sets the parameters related to the characterisation of the cell
};

#endif /* SRC_CELL_FIT_HPP_ */
