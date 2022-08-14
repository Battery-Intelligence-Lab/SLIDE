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

namespace slide {
template <int N>
struct CellData : public CellDataStorage<N>
{
  CellDataWriter<N> writer;

  auto writeData(auto *ths, const std::string &prefix)
  {
    writer.writeData(ths, prefix, *this);
  }
};
} // namespace slide
