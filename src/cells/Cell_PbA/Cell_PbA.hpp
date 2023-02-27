/*
 * Cell_PbA.hpp
 *
 *  Created on: 08 Jun 2022
 *   Author(s): Volkan Kumtepeli, Becky Perriment
 *
 * Model is taken from:
 * [1] Schiffer J, Sauer DU, Bindner H, Cronin T, Lundsager P, Kaiser R. Model prediction for ranking
 * lead-acid batteries according to expected lifetime in renewable energy systems and autonomous
 * power-supply systems. Journal of Power sources. 2007 May 25;168(1):66-78.
 */

#pragma once

#include "State_PbA.hpp"
#include "../Cell.hpp"
#include "../../settings/settings.hpp"
#include "../../utility/utility.hpp"

#include <string>
#include <vector>
#include <array>
#include <span>

namespace slide {
struct rho_param //!< from [1]
{
  using type = double;
  type a{ 0.4956 };   //!< g cm−3
  type b{ 0.2456 };   //!< g2 cm−6
  type c{ 77.53 };    //!< cm3
  type d{ -0.01278 }; //!< g cm−6
  type e{ 0.01289 };  //!< cm−3
  type f{ 0.0373 };   //!< g2 Ah−1 cm−6
};

class Cell_PbA : public Cell
{
protected:
  State_PbA st{ 0, 0.5, settings::T_ENV }; //!< I, T, SOC
  XYdata_vv k_OCVp;                        //!< OCVp vs corrosion speed mapping.
  double Rdc{ 2e-3 };                      //!< DC resistance [Ohm]

  //!< PbA parameters from [1]
  constexpr static double eps = 1e-11; //!< just a small number not to divide by zero.
  //!< SoC < 0  is possible.
  double offset = 0.86; //!< [V] offset parameter reported in the literature is typically 0.84–0.86 V.

  double SOC_infl = 10.0 / 13.0; //!< [-] SOC influence

  double Ucorr_0 = 1.75; //!< [V] Corrosion voltage of fully-charged battery without current flow, which is a function of the acid concentration
  double Igas_0 = 0.017; //!< [A] Normalized gassing current at Ugas,0 and Tgas,0
  //!< Typical value for the normalized gassing current for a new battery with antimony grid alloys. For VRLA battery typically 10 mA/100 Ah. For aged batteries the
  //!< normalized gassing can increase by a factor 5 or more.

  double c_u = 0.183;                  //!< [1/V] Voltage coefficient of gassing current
  double c_T = 0.06;                   //!< [1/K] Temperature coefficient of gassing current
  double Ugas_0 = 13.38;               //!< [V] Nominal voltage for gassing
  double Tgas_0{ 298 };                //!< [K] Nominal temperature for gassing
  double Tcorr_0{ 298 };               //!< [K] Nominal temperature for  corrosion
  double ks_T{ std::log(2.0) / 15.0 }; //!< [1 / K] Temperature coefficient of corrosion speed
  double cSOC_0{ 6.614e-5 };           //!< [1/h]
  double cSOC_min{ 3.307e-3 };         //!< [1/h]
  double SOC_limit{ 0.90 };            //!< [-] Minimum state-of-charge for bad charges
  double SOC_ref{ 0.95 };              //!< [-] Reference state-of-charge for bad charges
  double c_plus{ 1.0 / 30.0 };         //!< [-] Factor for increase of acid stratification
  double c_minus{ 0.1 };               //!< [-] Factor for decrease of acid stratification with gassing
  double Uref{ 2.5 };                  //!< [V] Reference voltage for decreasing acid stratification
  double Uacid_dec{ 2.3 };             //!< [V] Voltage at which gassing starts to remove acid stratification
  double D_H2SO4{ 20e-9 };             //!< [m^2/s] Diffusion constant for sulfuric acid
  double c_Z{ 5 };                     //!< [-] Exponent for calculation of capacity loss due to degradation
  double z_0{ 2.961e11 };              //!< [cm−3] Coefficient of number of sulfate crystals
  rho_param rhp{};

  //!< Battery dependent parameters:
  double Z_IEC{ 600 };                  //!< [cycles] Number of cycles under standard conditions (data sheet)
  double Lfloat{ 10 };                  //!< [years] Float lifetime (data sheet)
  double U0 = 3;                        //!< [V] open-circuit equilibrium cell voltage at the fully-charged state
  double g_OCV{ 0.076 };                //!< [V] electrolyte proportionality constant, Gradient of change in OCV with state-of-charge
  double rho_c{ 0.42 }, rho_d{ 0.699 }; //!< [Ohm*Ah] Effective internal resistance
  double M_c{ 0.888 }, M_d{ 0.0464 };   //!< [??] Resistance representing charge-transfer overvoltage coefficient which depends on state-of-charge
  double C_c{ 1.001 }, C_d{ 1.75 };     //!< [??] Normalised capacity of the battery.
  double Iref{ 55 };                    //!< [A] Normalized reference current for current factor //!< -55 in paper not needed.
  double h_batt{ 20 };                  //!< [cm] Height of battery (z in paper)

  double C_EOL{ 0.8 }; //!< End of life capacity. Cdeg,limit in paper. 80%.

  //!< Unknown parameters?

  double m_acid{};       //!< [??]  the weight (m) of the acid in the cell,
  double rho_nom{};      //!< nominal acid concentration in the battery
  double DeltaW_limit{}; //!< corrosion layer thickness when the battery has reached the end of its float lifetime (given in the battery datasheet).
  double C_corr_limit{}; //!< the limit of the loss of capacity by corrosion

  // From Becky:

  double V_w{ 17.5 }; //!< [cm^3/mol] molar volume of H2O
  double V_e{ 45 };   //!< [cm^3/mol] molar volume of H2SO4

  double M_w{ 18 }; //!< [g/mol] molar mass of H2O

public:
  Cell_PbA();
  Cell_PbA(std::string IDi, double capin, double SOCin);

  inline double SOC() override { return st.SOC(); }
  inline double I() override { return -st.I(); } //!< Inner representation have charging+ so we multiply with minus
  double V() override;

  void getStates(getStates_t s) override { s.insert(s.end(), st.begin(), st.end()); }         //!< returns the states of the cell collectively.
  std::span<double> viewStates() override { return std::span<double>(st.begin(), st.end()); } //!< returns the individual states.
  auto &getStateObj() { return st; }

  double getOCV(bool print = true) override; //!< crit is an optional argument
  //!< virtual int getNstates() { return S.size() + 1; } //!< +1 for current

  double getRtot() override { return Rdc; } //!< Return the total resistance, V = OCV - I*Rtot

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override;
  Status setSOC(double SOCnew, bool checkV = true, bool print = true) override;
  Status setStates(setStates_t s, bool checkV, bool print) override;

  //	void backupStates() override;  //!< Back-up states.
  //	void restoreStates() override; //!< restore backed-up states.

  //!< thermal model
  inline double T() override { return st.T(); }
  inline double getThotSpot() override { return T(); }
  double getThermalSurface() override { return 0; }; //!< Not implemented?
  inline void setT(double Tnew) override { st.T() = Tnew; }

  bool validStates(bool print = true) override;
  void timeStep_CC(double dt, int steps = 1) override;

  Cell_PbA *copy() override { return new Cell_PbA(*this); }

  //!< PbA specific functions:
  inline double DOD() { return (1 - SOC()); } //!< Returns the depth of discharge.
  double Ucorr();
  double rho_corr();
  double C_corr();
  double k_s();
  double k();

  double OCVp();

  double f_SOC();
  double f_acid();
  double f_stratification();

  double f_minus_gassing();
  double f_minus_diffusion();
  double f_minus() { return f_minus_gassing() + f_minus_diffusion(); }
  double f_plus();

  double f_I();

  double C_deg();

  double rho_empty(); //!< Eq. 4
  double I_gas();
  double Delta_n();

  bool isFullyCharged()
  {
    constexpr double fullyChargedSOC = 0.9999;
    return st.SOC() > fullyChargedSOC;
  }

  //!< virtual std::shared_ptr<StorageUnit> copy();

  //!< dataStorage
  //!< int increaseBin(double binedge[], double bin[], double data, int previndex);
  //!< virtual void storeData();
  //!< virtual void writeData(std::string prefix){}; //!< #TODO implement.
};

} // namespace slide