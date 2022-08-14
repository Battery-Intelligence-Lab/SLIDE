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

#include "Cell_SPM.hpp"

namespace slide {
class Cell_KokamNMC : public Cell_SPM
{
public:
  // constructors
  Cell_KokamNMC(Model_SPM *, int verbosei);
  Cell_KokamNMC(Model_SPM *, DEG_ID &, int verbosei);
};
} // namespace slide
