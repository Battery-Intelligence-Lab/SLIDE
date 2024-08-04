/**
 * @file Converter.hpp
 * @brief Power electronic converter class
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 10 Jun 2020
 * @details Based on paper from Patsios: two-stage converter (variable DC / fixed DC and fixed DC to fixed AC)
 */

#pragma once

namespace slide {
class Converter
{
private:
  double Vdc;  //!< voltage of the DC bus
  double Vac;  //!< voltage of the AC link
  double Pnom; //!< nominal power, [Wh]
public:
  Converter();

  void setPower(double Pnom);
  double getLosses(double Vin, double Iin);
};
} // namespace slide
