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

#pragma once

#include "cell.hpp"

class Cell_KokamNMC : public Cell
{
private:
public:
	// constructors
	Cell_KokamNMC(const struct slide::Model &, int verbosei);
	Cell_KokamNMC(const slide::Model &M, const DEG_ID &, int verbosei);
};
