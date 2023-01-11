/*
 * procedure_util.hpp
 *
 *  Created on: 08 Mar 2023
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../system/Battery.hpp"

namespace slide {
inline void visit_SUs(StorageUnit *su, auto &&fn)
{
  //!< visits storage units and applies function fn
  //!< function should be void as of now.
  if (auto c = dynamic_cast<Cell *>(su))
    fn(c); // It is cell, no submodules so apply function.
  else if (auto b = dynamic_cast<Battery *>(su)) {
    fn(b);                        // It is battery apply then call sub ones.
    visit_SUs(b->getCells(), fn); //!< if SU is a battery, recursively call this function on the module with the cells
  } else if (auto m = dynamic_cast<Module *>(su)) {
    //!< If su is a module, recursively call this function on its children
    fn(m);
    for (auto &su_ptr : m->getSUs())
      visit_SUs(su_ptr.get(), fn);
  }
}
} // namespace slide