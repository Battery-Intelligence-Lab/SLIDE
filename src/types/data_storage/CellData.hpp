/*
 * CellData.hpp
 *
 *  Created on: 13 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string>

#include "CellDataStorage.hpp"
#include "CellDataWriter.hpp"
#include "../../settings/enum_definitions.hpp"

namespace slide {
template <settings::cellDataStorageLevel N>
struct CellData : public CellDataStorage<N>
{
  CellDataWriter<N> writer;

  auto writeData(auto &cell, const std::string &prefix)
  {
    writer.writeData(cell, prefix, *this);
  }
};
} // namespace slide
