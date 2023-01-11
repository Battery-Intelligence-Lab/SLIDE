/*
 * Converter.cpp
 *
 *  Created on: 10 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Converter.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ctime>

namespace slide {
Converter::Converter()
{
  Vdc = 1500;   //!< should be larger than max battery voltage (15*20*4.2 = 1260)
  Vac = 1500;   //!< random guess, 5 * 380V (line voltage)
  Pnom = 200e3; //!< 200 kWh from Schimpe etc.
}
void Converter::setPower(double Pn)
{
  Pnom = Pn;
}

double Converter::getLosses(double Vin, double Iin)
{
  using PhyConst::pi;
  /*
   * function to calculate the losses of the converter [W]
   * Note that even as the input power goes to 0, the losses won't go below 1.5 kW due to the switching losses and losses on the DC bar and in the filter
   * all of which are independent on the input current or voltage
   *
   * IN
   * Vin		input voltage [V], i.e. the voltage of the battery
   * Iin 		input current [A], i.e. the current of the battery
   *
   * OUT
   * loss 	power lost in the conversion process [W]
   */

  /*
   * Parameter values are from Patsios et al, Schimpe et al. or Schimpe's reference [69]
   * Their batteries are smaller, e.g. [69] is a 200 kWh inverter.
   * Assume conduction losses are the same formula (since they get multiplied by current or voltage), but scale switching losses up by the number of inverters we would need
   */

  //!< Unless stated otherwise, values are from Patsios et al
  //!< double Iinold = Iin; //!< for printing

  //!< We have muliple inverters, split the current between them
  constexpr double Pschimpe = 200 * 1e3;
  const int nInverters = static_cast<int>(std::ceil(Pnom / Pschimpe)); //!< the number of inverters needed to convert the power
  Iin /= nInverters;

  //!< voltage drop over IGBTs from [69] from Schimpe
  //!< https://www.infineon.com/dgdl/Infineon-FS150R12KE3-DS-v03_01-en_de.pdf?fileId=db3a304412b407950112b4311d745388
  //!< 	fig on p4, output characteristic IGBT at Vge = 15 and T = 125 (top right figure)
  //		linear approximation: V = 1 -> I = 30 while V = 3 -> I = 300
  //!< 							  V = 1 + (3-1)/(300-30)*(I - 30)
  //!< 		and below 30A: V = I/30
  //!< 	V = a + b (std::abs(I) - c)

  double a{ 0 }, b{ 1.0 / 30.0 }, c{ 0 };
  if (std::abs(Iin) > 30) //!< current through rectifier = current of battery
  {
    a = 1;
    b = 2.0 / 270.0;
    c = 30;
  }

  //!< Switching energy, based on same doc, right bottom
  //!< 	I = 50 -> E = 6 and I = 200 -> E = 21
  //!< 	so E = 6 + 15.0/150.0*(I-50) [miliJoule] for Eon.
  //!< 	Eoff = Eon + 1
  constexpr double as{ 6 }, bs{ 0.1 }, cs{ 50 };

  //!< Note, that ref is a one-stage converter.
  //!< 	The values for diodes are the by-pass diodes so might not be representative
  //!< 	so just use values for the IGBTs for DC/DC converter too

  //!< DC stage
  constexpr double f = 8000;                           //!< switching frequency, value from Schimpe
  double Vce = a + b * (std::abs(Iin) - c);            //!< voltage across the diode (i.e. resistive voltage drop across it, random guess about 0.1%)
  double D = 1 + 2 * Vin / Vdc;                        //!< modulation ratio
  double Eon = 0;                                      //!< switch-on loss. 0 since this is a diode
  double Eoff = (as + bs * std::abs(Iin - cs)) * 1e-3; //!< energy loss for switch off
  double DC_cond = Iin * Vce * D;                      //!< conduction losses
  double DC_switch = f * (Eoff + Eon);                 //!< switching losses
  const double Idc = Iin * Vin / Vdc;                  //!< current on the DC bus (ignoring losses)

  //!< AC stage
  if (std::abs(Idc) > 30) { //!< current through inverter = current of the DC link
    a = 1;
    b = 2.0 / 270.0;
    c = 30;
  } else {
    a = 0;
    b = 1.0 / 30.0;
    c = 0;
  }
  Vce = a + b * (std::abs(Idc) - c);
  D = 1 + 2 * Vdc / Vac; //!< modulation ratio
  Eon = (as + bs * (std::abs(Idc - cs))) * 1e-3;
  Eoff = Eon + 1e-3;                       //!< https://www.infineon.com/dgdl/Infineon-FS150R12KE3-DS-v03_01-en_de.pdf?fileId=db3a304412b407950112b4311d745388
  double AC_cond = Idc * Vce * D;          //!< conduction losses (already has full current so already account for all 6 IGBTs)
  double AC_switch = 6 * f * (Eoff + Eon); //!< switching losses: 6 IGBTs for 3 phase full bridge converter
  //!< double Iac = Idc * Vdc / Vac;		 //!< current on the AC side (ignoring losses)

  //!< filter
  constexpr double R = 13.6e-3;                              //!< Value from Schimpe (for resistance of grid interface)
  constexpr double Cdc = 3e-3;                               //!< DC link capacitance
  constexpr double dV = 1;                                   //!< DC ripple voltage, wild guess since I can't find a value #TODO if is correct?
  constexpr double Rg = 13.6e-3;                             //!< value from Schimpe (resistance to grid)
  constexpr double Ri = 13.6e-3;                             //!< value from Schimpe (resistance to inverter)
  constexpr double FIL_cdc = R * sqr(2 * pi * f * Cdc * dV); //!< losses in capacitor on DC bus
  const double FIL_lg = 3 * Rg * sqr(Idc);                   //!< losses in L at the output of the inverter
  const double FIL_lb = Ri * sqr(Iin);                       //!< losses in L at the input of the rectifier

  //!< total losses
  DC_cond = std::abs(DC_cond); //!< ensure all have same sign (some I^2, others I so if I < 0 they would have opposite sign)
  DC_switch = std::abs(DC_switch);
  AC_cond = std::abs(AC_cond);
  AC_switch = std::abs(AC_switch);
  const double loss_per_inverter = DC_cond + DC_switch + AC_cond + AC_switch + FIL_cdc + FIL_lg + FIL_lb;
  const auto loss = loss_per_inverter * nInverters;

  //!< cout<<"Converter Vin = "<<Vin<<", Iin = "<<Iinold<<", ninverters are "<<nInverters<<" losses are "<<loss/std::abs(Vin*Iinold)*100<<" %"<<endl;

  /*
  std::cout<<"\t Converter Pnom = "<<Pnom<<" giving "<<nInverters<<" inverters"<<endl;
  double fact = std::abs(100.0/(Vin*Iin)); //!< print relative losses in [%]
  std::cout<<"\t contributions DC: "<<DC_cond*fact<<" "<<DC_switch*fact<<" AC: "<< AC_cond*fact<<" "<<AC_switch*fact<<" filter: "<<FIL_cdc*fact<<" "<<FIL_lg*fact<<" "<< FIL_lb*fact<<" total "<<loss*fact<<endl;
  //fact = 1.0/1000.0; //!< print total losses in [kW]
  //cout<<"\t contributions DC: "<<DC_cond*fact<<" "<<DC_switch*fact<<" AC: "<< AC_cond*fact<<" "<<AC_switch*fact<<" filter: "<<FIL_cdc*fact<<" "<<FIL_lg*fact<<" "<< FIL_lb*fact<<" total "<<loss*fact<<endl;
*/
  return loss;
}
} // namespace slide
