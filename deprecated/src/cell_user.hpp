/*
 * Cell_user.hpp
 *
 * A child classes of Cell which the user can use to set his/her own parameters
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "Cell_LGChemNMC.hpp"
#include "settings/settings.hpp"

namespace slide {
class Cell_user : public Cell
{
public:
  //!< constructors
  Cell_user(const struct slide::Model_SPM &, int verbosei);
  Cell_user(const slide::Model_SPM &M, const DEG_ID &, int verbosei);
};

} // namespace slide
