/*
 * Cell_KokamNMC.hpp
 *
 * One of the child classes of Cell which implements a real cell.
 * The cycling parameters are for a high power 18650 NMC cell manufactured by Kokam.
 * The degradation parameters are set such that each mechanism clearly affects the battery life (which is not the case in reality).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef Cell_KokamNMC_HPP_
#define Cell_KokamNMC_HPP_

#endif /* CELL_HPP_ */


#include <string>
#include <cstring>

#include "Cell.hpp"


class Cell_KokamNMC : public Cell{
private:

public:
	// constructors
		Cell_KokamNMC(const struct Model&, int verbosei);
		Cell_KokamNMC(const Model& M, const DEG_ID&, int verbosei);
	virtual ~Cell_KokamNMC();
};

