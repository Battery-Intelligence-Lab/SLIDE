/**
 * @file CoolSystem_open.hpp
 * @brief Declaration of the CoolSystem_open class.
 *
 * This class represents a coolsystem with an 'open' heat exchanger.
 * It is a 'pass through' cool system where there is no barrier between its child SUs and parent module.
 *
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 5 Jun 2020
 */

#pragma once

#include "CoolSystem.hpp"

namespace slide {

/**
 * @brief A coolsystem with an 'open' heat exchanger.
 *
 * This class implements a coolsystem with a very high convective cooling constant,
 * such that its temperature is almost equal to the temperature of its children.
 * Additionally, the flowrate is 0, so there is no energy necessary to spin a fan.
 *
 * Note: In reality, we would probably prefer the temperature to be close to the temperature
 * of its parent instead of the children (i.e., the children access the air from the parent).
 * However, this is difficult to implement since it requires the cooling constant of its parent
 * to be very high, which violates the rules of object-oriented programming.
 * Therefore, we chose the opposite approach and consider this coolsystem as an aggregation of
 * its children, and the parent cools something which is at the temperature of its children.
 */
class CoolSystem_open : public CoolSystem
{
protected:
  double h{ 90 }; ///< Cooling constant for perfect heat exchange.

public:
  /**
   * @brief Default constructor for the CoolSystem_open class.
   */
  CoolSystem_open();

  /**
   * @brief Constructor for the CoolSystem_open class.
   * @param Ncells The number of cells.
   * @param control The control parameter.
   */
  CoolSystem_open(int Ncells, int control);

  /**
   * @brief Get the value of the cooling constant (h).
   * @return The value of the cooling constant (h).
   */
  double getH() override;

  /**
   * @brief Control the cooling system based on local and global hot temperatures.
   * @param Thot_local The local hot temperature.
   * @param Thot_global The global hot temperature.
   */
  void control(double Thot_local, double Thot_global) override;

  /**
   * @brief Create a copy of the CoolSystem_open object.
   * @return A pointer to the newly created CoolSystem_open object.
   */
  CoolSystem_open *copy() override { return new CoolSystem_open(*this); }
};

} // namespace slide