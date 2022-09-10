/*
 * CoolSystem_open.cpp
 *
 *  Created on: 5 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "CoolSystem_open.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace slide {
CoolSystem_open::CoolSystem_open() : CoolSystem()
{
  flowrate = 0;
}
CoolSystem_open::CoolSystem_open(int Ncells, int control) : CoolSystem(Ncells, control), h(90)
{
  flowrate = 0;
}

double CoolSystem_open::getH() { return h; }

void CoolSystem_open::control(double Thot_local, double Thot_global)
{
  flowrate = 0;
}
} // namespace slide
