/*
 * Cell_user.hpp
 *
 * A child classes of Cell which the user can use to set his/her own parameters
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_CELL_USER_HPP_
#define SRC_CELL_USER_HPP_


#include <string>
#include <cstring>

#include "Cell_LGChemNMC.hpp"

namespace std {

class Cell_user : public Cell {
public:
	// constructors
		Cell_user(const struct Model&, int verbosei);
		Cell_user(const Model& M, const DEG_ID&, int verbosei);
	virtual ~Cell_user();
};

} /* namespace std */

#endif /* SRC_CELL_USER_HPP_ */
