/*
 * CoolSystem_open.hpp
 *
 * This class is a coolsystem with an 'open' heat exchanger.
 * I.e. there is no barrier between its child SUs and parent module.
 * It is a 'pass through' cool system.
 *
 * This is implemented by having a very high convective cooling constant, such that its temperature is almost equal to the temperature of its children.
 * Additionally the flowrate is 0, so there is no energy necessary to spin a fan.
 *
 * Note: in reality we would probably prefer the temperature to be close to the temperature of its parent instead of the children (i.e. the children access the air from the parent).
 * This is however difficult to implement, since it requires the cooling constant of its parent to be very high.
 * That violates the rules of object oriented programming, since the coolsystem of the parent should be independent of what we are implementing here.
 * therefore, chose the opposite and say this coolsystem is behaving like an aggregation of its children, and the parent cools something which is at the temperature of its children.
 *
 *  Created on: 5 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "CoolSystem.hpp"

namespace slide {
class CoolSystem_open : public CoolSystem
{
protected:
  double h{ 90 }; //!< cooling constant for perfect heat exchange
public:
  CoolSystem_open();
  CoolSystem_open(int Ncells, int control);

  double getH() override;
  void control(double Thot_local, double Thot_global) override;

  CoolSystem_open *copy() override { return new CoolSystem_open(*this); }
};
} // namespace slide