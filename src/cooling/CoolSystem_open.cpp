/**
 * @file CoolSystem_open.cpp
 * @brief Implementation of the CoolSystem_open class.
 *
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 5 Jun 2020
 */

#include "CoolSystem_open.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace slide {

/**
 * @brief Default constructor for the CoolSystem_open class.
 */
CoolSystem_open::CoolSystem_open() : CoolSystem()
{
  flowrate = 0;
}

/**
 * @brief Constructor for the CoolSystem_open class.
 * @param Ncells The number of cells.
 * @param control The control parameter.
 */
CoolSystem_open::CoolSystem_open(int Ncells, int control) : CoolSystem(Ncells, control), h(90)
{
  flowrate = 0;
}

/**
 * @brief Get the value of the h parameter.
 * @return The value of h.
 */
double CoolSystem_open::getH() { return h; }

/**
 * @brief Control the cooling system based on local and global hot temperatures.
 * @param Thot_local The local hot temperature.
 * @param Thot_global The global hot temperature.
 */
void CoolSystem_open::control(double Thot_local, double Thot_global)
{
  flowrate = 0;
}

} // namespace slide