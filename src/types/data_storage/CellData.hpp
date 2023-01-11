/*
 * CellData.hpp
 *
 *  Created on: 13 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "CellDataStorage.hpp"
#include "CellDataWriter.hpp"
#include "../../settings/enum_definitions.hpp"

#include <string>

namespace slide {
template <settings::cellDataStorageLevel N>
struct CellData : public CellDataStorage<N>
{
  auto writeData(auto &cell, const std::string &prefix)
  {
    CellDataWriter<N>::writeData(cell, prefix, *this);
  }
};
} // namespace slide
