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

#pragma once

#include "Cell.hpp"

namespace slide {
class Cell_LGChemNMC : public Cell_SPM
{
private:
public:
  Cell_LGChemNMC(Model_SPM *, int verbosei); //!< constructor
  Cell_LGChemNMC(Model_SPM *M, DEG_ID &, int verbosei);
};
} // namespace slide
