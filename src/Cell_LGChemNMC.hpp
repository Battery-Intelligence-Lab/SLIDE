/*
 * Cell_LGChemNMC.hpp
 *
 * One of the child classes of Cell which implements a real cell.
 * The cycling parameters are for a high energy 18650 NMC cell manufactured by LG Chem.
 * The degradation parameters are set such that each mechanism clearly affects the battery life (which is not the case in reality).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef Cell_LGCHEMNMC_HPP_
#define Cell_LGCHEMNMC_HPP_

#endif /* CELL_HPP_ */


#include <string>
#include <cstring>

#include "Cell_KokamNMC.hpp" 				// includes Cell.hpp, State.hpp and Model.hpp


class Cell_LGChemNMC : public Cell{
private:

public:
	// constructor
		Cell_LGChemNMC(const struct Model&, int verbosei);
		Cell_LGChemNMC(const Model& M, const DEG_ID&, int verbosei);

	virtual ~Cell_LGChemNMC();
};

