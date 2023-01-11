/*
 * Cell_SPM_degradation.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell_SPM.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>

namespace slide {
void Cell_SPM::SEI(double OCVnt, double etan, double *isei, double *den)
{
  /*
   * Function to calculate the degradation effects of growth of the SEI layer
   *
   * IN
   * OCVnt 	the OCV of the negative electrode at the battery temperature [V]
   * etan 	the overpotential at the negative electrode [V]
   *
   * OUT
   * isei 	current density for the SEI side-reaction [A m-2]
   * den 		decrease in the active volume fraction as a result of SEI growth [sec-1]
   *
   * THROWS
   * 106 		illegal value in id or por
   */

  using namespace PhyConst;
  using std::exp;

  //!< variables
  double is{ 0 }; //!< SEI side reaction current density of all models combined

  const auto ArrheniusCoeff = calcArrheniusCoeff();
  //!< Loop for each model to use
  for (const auto &sei_id : deg_id.SEI_id) {
    //!< Use a switch to calculate the magnitude of the SEI growth according to this degradation model
    switch (sei_id) {
    case 0: //!< no SEI growth
      is += 0;
      break;
    case 1: //!< kinetics and diffusion according to Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
    {
      const auto kseit = sei_p.sei1k * exp(sei_p.sei1k_T * ArrheniusCoeff); //!< Arrhenius relation for the rate parameter at the cell temperature
      is += nsei * F * kseit * exp(-nsei * F / (Rg * st.T()) * alphasei * (OCVnt + etan - OCVsei + rsei * st.delta() * I()));
      //!< Add the effect of this model
      //!< eta_sei = OCVneg + etaneg - OCVsei + rsei*I
      //!< isei = nFk exp(-nF/RT alpha eta_sei)
      //!< on charge, I < 0 and etan < 0.
      //!< so higher charging current -> more negative term in exponential -> larger isei
    } break;
    case 2: //!< Kinetic model according to Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004 #TODO -> In slidepack case1 paper and case2 paper are swapped.
    {
      const auto kseit = sei_p.sei2k * exp(sei_p.sei2k_T * ArrheniusCoeff); //!< Arrhenius relation for the rate parameter at the cell temperature
      const auto Dseit = sei_p.sei2D * exp(sei_p.sei2D_T * ArrheniusCoeff); //!< Arrhenius relation for the diffusion constant at the cell temperature

      //!< derivation of the formula:
      //!< start equations (use the symbols from Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
      //!< but with opposite sign for j
      //!< j = nFk c exp(..)
      //!< j/nF = - D/delta (c - c0)
      //!< j/nF = D/delta (- j/(nFk exp(..)) + c0)
      //!< j * ( 1/(nFk exp(..)) + delta/DnF ) = c0
      //!< j = c0 / ( 1/(nFk exp(..)) + delta/DnF )
      const auto isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * st.T()) * alphasei * (OCVnt + etan - OCVsei + rsei * st.delta() * I()));
      const auto isei3 = st.delta() / (nsei * F * Dseit);
      is += c_elec0 / (1.0 / isei2 + isei3); //!< Add the effects of this model
    } break;
    case 3: //!< model from Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
    {
      const auto kseit = sei_p.sei3k * exp(sei_p.sei3k_T * ArrheniusCoeff); //!< Arrhenius relation for the rate parameter at the cell temperature
      const auto Dseit = sei_p.sei3D * exp(sei_p.sei3D_T * ArrheniusCoeff); //!< Arrhenius relation for the diffusion constant at the cell temperature
      //!< Use equation [22] from the paper
      constexpr double a_L_K = 0.134461; //!< the parameter a_L_K is set to 0.134461 but this constant can be lumped into the rate- and diffusion constants
      const auto isei1 = a_L_K * exp(-nsei * F * (etan + rsei * st.delta() * I()) / (Rg * st.T()));
      const auto isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * st.T()) * alphasei * (OCVnt - OCVsei));
      const auto isei3 = st.delta() / (nsei * F * Dseit);
      is += isei1 / (1.0 / isei2 + isei3); //!< Add the effects of this model
    } break;
    case 4: //!< model from the optimisation in the paper
    {
      const auto kseit = sei_p.sei4k * exp(sei_p.sei4k_T * ArrheniusCoeff); //!< Arrhenius relation for the rate parameter at the cell temperature
      const auto Dseit = sei_p.sei4D * exp(sei_p.sei4D_T * ArrheniusCoeff); //!< Arrhenius relation for the diffusion constant at the cell temperature
      //!< Use equation [22] from the paper
      constexpr double a_L_K = 0.134461;                                //!< the parameter a_L_K is set to 0.134461 but this constant can be lumped into the rate- and diffusion constants
      const auto isei1 = a_L_K * exp(-nsei * F * etan / (Rg * st.T())); //!< note: model used in optimisation had kpt/knt or vice versa, here a fixed value
      const auto isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * st.T()) * alphasei * (OCVnt - OCVsei));
      const auto isei3 = st.delta() / (nsei * F * Dseit);
      is += isei1 / (1.0 / isei2 + isei3); //!< Add the effects of this model
    } break;
    default: //!< unknown degradation model
      std::cerr << "ERROR in Cell_SPM::SEI, unknown SEI degradation model with identifier "
                << static_cast<int>(sei_id) << ". Only values 0 to 3 are allowed. Throw an error.\n";
      throw 106;
      break;
    }
  } //!< end loop for all the models you want to use

  //!< Make the output for the SEI side reaction current density
  *isei = is;

  //!< Calculate how much we decrease the volume fraction due to SEI growth
  if (deg_id.SEI_porosity == 0) //!< don't decrease volume fraction
    *den = 0;
  else if (deg_id.SEI_porosity == 1) {                                 //!< decrease volume fraction according to Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
    double jn = I() / geo.elec_surf / (st.an() * n * F * st.thickn()); //!< molar flux on the negative particle
    *den = -sei_p.sei_porosity * (jn * Vmain + *isei * Vsei);
    //!< note: they use J = volumetric current [A/m3] -> they multiply with 'an' but we already have density [A/m2]
    //!< - because they use the porosity while we use the volume fraction
  } else { //!< unknown degradation model
    std::cerr << "ERROR in Cell_SPM::SEI, unknown value for decreasing the volume fraction "
              << deg_id.SEI_porosity << ". Only values 0 or 1 are allowed. Throw an error.\n";
    throw 106;
  }
}

void Cell_SPM::CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn)
{
  /*
   * function to calculate the degradation effect of surface cracking due to fatigue
   *
   * IN
   * OCVnt 			the OCV of the negative electrode at the battery temperature [V]
   * etan 			the overpotential at the negative electrode [V]
   *
   * OUT
   * isei_multiplyer 	Extra SEI side reaction due to crack growth as a fraction of the original SEI side reaction current [-].
   * 						i.e. total SEI growth = (1+isei_multiplyer)*isei
   * dCS 				increase in surface area due to fatigue in this time step [m2 sec-1]
   * dDn 				decrease in the diffusion constant of the graphite due to surface cracks [m s-1 s-1]
   *
   * THROWS
   * 106 				illegal value in id or d
   * 107				too many degradation models
   * 108 				the stress values are not up to date
   */

  using namespace PhyConst;

  //!< output parameters
  double ism = 0;                     //!< isei multiplier from all models to be considered
  double dcs = 0;                     //!< increase in surface area from all models to be considered
  const auto ASn = getAnodeSurface(); //!< active surface area of the anode [m2]
                                      //!< this active area is used to translate absolute values (such as currents) to relative values (such as current densities)

  //!< Loop for each model we want to use
  for (const auto cs_id : deg_id.CS_id) {
    //!< a switch to calculate the effect according to model i
    switch (cs_id) {
    case 0: //!< no surface cracks
      break;
    case 1: //!< Laresgoiti's stress and crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
            //!< this model calculates crack growth due to temporal variations in the li-concentration
      //!< check the calculated stress values are up to date
      if (!sparam.s_lares_update) {
        std::cerr << "ERROR in Cell_SPM::CS. The stress values for Laresgoiti's stress model are not updated. Throwing an error.\n";
        //!< if you see this error, you have to call Cell_SPM::updateLaresgoitiStress(), which calculates the stress and stores the values, before you call this function
        throw 108;
      }

      //!< Implement the equation from the paper
      //!< capacity loss is m-power of the so-called stress amplitude (sigma_max - sigma_min)/2
      //!< sigma_max and sigma_min are the max and min stresses 'of the cyclic signal' i.e. within one charge/discharge
      //!< assume m = 1, then (max - min) = (max - t1) + (t1-t2) + (t2-t3) + ... + (tn - min)
      //!< so stress amplitude can be substituted by (the stress in the previous time step) - (the stress in this time step)
      dcs += csparam.CS1alpha * std::sqrt(std::abs(sparam.s_lares_n - sparam.s_lares_n_prev) / sparam.s_dt);
      //!< equations (22)+ (27) from the paper
      //!< current density on particle = I /(elec_surf * thick * a)
      //!< 		isei also acts on this scale since it is an extra boundary condition ( itot = (jn + isei) =  (surface gradient)/nF )
      //!< 		crack growth -> increase isei_on_particle -> (jn + isei*(initial+crack_surface)/initial)*nF
      //!< 			such that if the crack surface area is the same as the initial electrode surface area, we double isei
      ism += st.CS() / ASn; //!< increase SEI growth proportionally the crack surface
      break;
    case 2: //!< Laresgoiti's crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
            //!< but with stress model from Dai, Cai, White, Journal of Power sources 247, 2014
            //!< instead of Laresgoiti's stress from figure 5
            //!< this model calculates crack growth due to spatial (Dai) and temporal (Laresgoiti) variations in the li-concentration
      //!< ensure the stress values are up to date
      if (!sparam.s_dai_update) {
        std::cerr << "ERROR in Cell_SPM::CS. The stress values for Dai's stress model are not updated. Throwing an error.\n";
        //!< if you see this error, you have to call Cell_SPM::updateDaiStress(), which calculates the stress and stores the values before you call this function
        throw 108;
      }

      //!< Add the effects of this model
      dcs += csparam.CS2alpha * std::sqrt(std::abs(sparam.s_dai_n - sparam.s_dai_n_prev) / sparam.s_dt);
      //!< equations (22)+ (27) from the paper but with Dai's stress
      ism += st.CS() / ASn; //!< increase SEI growth proportionally the crack surface
      break;
    case 3: //!< model by Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
            //!< this model is adapted to calculate crack growth due to spatial variations in the li-concentration
      //!< get concentrations

      double cp[settings::nch + 2], cn[settings::nch + 2];
      getC(cp, cn);
      //!< Add the effects of this model
      dcs += csparam.CS3alpha * sqr((cn[0] - cn[settings::nch + 1]) / Cmaxneg);
      //!< equations (8) + (21)
      //!< Note that eqn (8) refers to the change with respect to time while here we use the spatial variation
      //!< This makes the model capture more of the spatial variation in stress
      //!< Laresgoiti's model already accounted for temporal variation, so simply use that if you are interested in temporal rather than spatial variation
      ism += st.CS() / ASn; //!< increase SEI growth proportionally the crack surface
      break;
    case 4: {
      //!< model from Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
      //!< equation (1a): CS = Amax(1 - exp(-m * Ah)) with m a fitting parameter and Ah the charge throughput up to now
      //!< 	dCS/dAh = m Amax exp(-m Ah) = m Amax - m Amax + m Amax exp(-m Ah) = m Amax - m CS = m(Amax - CS)
      //!< 	dCS/dt = dCS/dAh * dAh/dt = dCS/dAh * abs(I) = m(Amax - CS)*abs(I)
      //!< 	where we use the absolute value of I because 'Ah' is the total charge throughput, i.e. int ( abs(I) dt )
      const double Amax = std::max(csparam.CS4Amax, st.CS());
      //!< 'maximum crack surface area', a fitting parameters
      //!< avoid negative crack growth if the crack surface becomes slightly larger than Amax
      //!< this is possible due to discrete time steps: CS(t) is just smaller, but CS (t+1) = CS(t) + dCS*dt is just larger

      //!< Add the effects of this model
      dcs += csparam.CS4alpha * (Amax - st.CS()) * std::abs(I()); //!< see above, with m = csparam.CS4
      ism += st.CS() / ASn;                                       //!< increase SEI growth proportionally the crack surface
    } break;
    case 5: {
      //!< model from Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
      //!< overpotential for the crack-side-reaction = overpotential for the SEI reaction
      const double etasei = (OCVnt + etan - OCVsei + rsei * st.delta() * I()); //!< overpotential [V], equation (6)

      //!< get surface concentration
      double cps, cns;
      getCSurf(cps, cns, true); //!< get the surface lithium concentration //!< #TODO only cns is used.

      double kcr; //!< rate constant for the side reaction
      //!< Calculate the rate constant, equation (11) with an Arrhenius relation for the temperature (which wasn't considered by Ekstrom)
      if (isDischarging())
        kcr = 0;
      else if (cns / Cmaxneg < 0.3)
        kcr = 2 * csparam.CS5k * exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / st.T()));
      else if (cns / Cmaxneg < 0.7)
        kcr = 0;
      else
        kcr = csparam.CS5k * exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / st.T()));

      //!< Add the effects of this model
      dcs += nsei * F * kcr * exp(-alphasei * nsei * F / (Rg * st.T()) * etasei); //!< equation (9)
      ism += st.CS() / ASn;                                                       //!< increase SEI growth proportionally the crack surface

    } break;
    default: //!< unknown degradation model
      std::cerr << "ERROR in Cell_SPM::CS, unknown crack growth model with identifier "
                << cs_id << ". Only values 0 to 5 are allowed. Throw an error.\n";
      throw 106;
    }
  }
  //!< Make the output variables
  *isei_multiplyer = ism;
  *dCS = dcs;

  //!< Decrease the negative diffusion constant if needed
  if (deg_id.CS_diffusion == 0) //!< don't decrease the negative diffusion constant
    *dDn = 0;
  else if (deg_id.CS_diffusion == 1) { //!< decrease it according to Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
    //!< equation (2) D(t) = D0 (1 - CS)^gamma
    //!< 	but this can become negative is CS is larger than 1, we assume there should be a term /Amax in eqn (2): D(t) = D0 (1 - (CS/Amax))^gamma, which becomes 0 if CS = Amax
    //!< so dD/dt = - gamma D0 (1-CS/Amax)^(gamma-1) 1/Amax dCS/dt

    //!< 'maximum crack surface area', a fitting parameters
    //!< avoid increasing diffusion coefficient if the crack surface becomes larger than Amax
    //!< this is possible if the user chooses a different CS growth model, which does give larger crack surfaces (i.e. not CS4 which is Barai's crack growth model)
    const double Amax = std::max(csparam.CS4Amax, st.CS());
    //!< cap the decrease rate at a maximum value of 2e-7 (2e-5% = kill the battery in  about 1000h)
    const double Dnmax = std::min(2e-7, csparam.CS_diffusion * std::pow(1 - st.CS() / Amax, csparam.CS_diffusion - 1) / Amax * dcs);
    *dDn = -Dnmax * st.Dn();
  } else { //!< unknown degradation model
    std::cerr << "ERROR in Cell_SPM::CS, unknown value for decreasing the diffusion constant "
              << deg_id.CS_diffusion << ". Only values 0 or 1 are allowed. Throw an error.\n";
    throw 106;
  }
}

void Cell_SPM::LAM(bool print, double zp_surf, double etap, double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den)
{
  /*
   * Function to calculate the effect of loss of active material (LAM).
   *
   * LAM is simulated by decreasing the amount of active material.
   * This will increase the current density on the particle for the same overall cell current (because there is less material to 'spread' it over).
   * This will mean that for the same overall cell current,
   * 		there will be a larger change in lithium concentration,
   * 			so there is a larger change in open circuit voltage,
   * 			so a smaller capacity before a voltage limit is reached
   * 			so the capacity decreases,
   * 		and there will be a larger resistive voltage drop, so the 'effective' resistance increases
   *
   * There are a couple of variables describing the amount of active material:
   * 	elec_surf 	the (geometric) surface of the electrode, i.e. the product of the height and length of the electrode [m2]
   * 	thick		the thickness of the electrode material [m]
   * 	R 			radius of the particle of the single particle model [m]
   * 	e			the volume fraction of active material	[-]
   * 	a 			the effective surface, i.e. the surface per unit of electrode volume [m2 m-3]
   * 				a = 3*e/R
   *
   * 	The current density is calculated as:
   * 	i = I / (a * thick * elec_surf) = I / (3* e / R * thick * elec_surf)
   *
   * 	A decrease in any of these geometric parameters will have exactly the same effect on the battery:
   * 	you can double i by halving a, or by halving thick, or by halving elec_surf.
   * 	and you can halve a by halving e or doubling R.
   * 	So it doesn't really matter which of the geometric parameters you decrease to account for LAM.
   * 	However, this model doesn't allow to change elec_surf and R because these values are needed on multiple locations in the code.
   * 	E.g. the value of R is needed to calculate values of the matrices used for the diffusion state space model (because the spatial discretisation depends on R).
   * 	So if you were to change R, you have to re-calculate the matrices, which would needlessly complicate the model.
   * 	Therefore, degradation can only decrease the values of 'thick', 'a' and 'e'
   *
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   * zp_surf 		li-fraction at the surface of the positive particle [-]
   * etap 		overpotential at the positive electrode [V]
   *
   * OUT
   * dthickp 		change in electrode thickness of the positive electrode [m s-1]
   * dthickn		change in electrode thickness of the negative electrode [m s-1]
   * dap			change in effective electrode surface of the positive electrode [m2 m-3 s-1]
   * dan			change in effective electrode surface of the negative electrode [m2 m-3 s-1]
   * dep			change in volume fraction of active material in the positive electrode [s-1]
   * den			change in volume fraction of active material in the negative electrode [s-1]
   *
   * THROWS
   * 106 			illegal value in id
   * 108 			the stress values are not updated
   */

  using namespace PhyConst;
  using std::exp;

  //!< output parameters
  double dthickpp{}, dthicknn{}, dapp{}, dann{}, depp{}, denn{};

  //!< loop for each model to use
  for (const auto lam_id : deg_id.LAM_id) {
    const auto ArrheniusCoeff = calcArrheniusCoeff();
    //!< calculate the effect of this model
    switch (lam_id) {
    case 0: //!< no LAM
      break;
    case 1: //!< Stress model from Dai, Cai, White, Journal of Power sources 247, 2014
            //!< LAM equation similar to CS equation from Laresgoiti
            //!< (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
      //!< ensure the stress values are up to date
      if (!sparam.s_dai_update) {
        std::cerr << "ERROR in Cell_SPM::LAM. The stress values for Dai's stress model are not updated. Throwing an error.\n";
        //!< if you see this error, you have to call Cell_SPM::updateDaiStress(), which calculates the stress and stores the values before calling this function
        throw 108;
      }

      //!< Laresgoiti's equation to link stress to LAM
      dthickpp += -lam_p.lam1p * std::abs(sparam.s_dai_p - sparam.s_dai_p_prev) / sparam.s_dt; //!< with Dai's stress model (values stored in s_dai_p)
      dthicknn += -lam_p.lam1n * std::abs(sparam.s_dai_n - sparam.s_dai_n_prev) / sparam.s_dt; //!< #TODO sparam.s_dt was 2.0 in slide, why?
      //!< ageing fit
      //!< you need to divide by the time step to counter the effect of the time step
      //!< 	larger dt has double effect: increase the stress difference, and time integrate by larger time period
      //!< 	assuming stress changes are constant at ds per second, it would be ds (time_now - time_prev), or ds * dt
      //!< 	and to cover a period of T (so we need T/dt steps of dt each), the total effect of stress is (ds*dt) * dt * T/dt = ds*dt*T
      //!< 	so the effect of stress increases linearly with the size of the time setp
      //!< 	so divide by total time (dt*nstep) to cancel this out, i.e. ds is per second (and not per total time period)
      //!< assume the other effects are 0
      break;
    case 2: //!< Model by Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
    {       //!< Get the molar flux on each particle

      const auto [i_app, jp, jn] = calcMolarFlux(); //!< current density, molar flux on the pos/neg particle
      //!< Use Arrhenius relations to update the fitting parameters for the cell temperature
      const double ap = lam_p.lam2ap * exp(lam_p.lam2t * ArrheniusCoeff);
      const double an = lam_p.lam2an * exp(lam_p.lam2t * ArrheniusCoeff);
      const double bp = lam_p.lam2bp * exp(lam_p.lam2t * ArrheniusCoeff);
      const double bn = lam_p.lam2bn * exp(lam_p.lam2t * ArrheniusCoeff);

      //!< Add the effects of this model
      const auto abs_jp{ std::abs(jp) }, abs_jn{ std::abs(jn) };
      depp += ap * abs_jp + bp * std::sqrt(abs_jp); //!< equation (5) from the paper
      denn += an * abs_jn + bn * std::sqrt(abs_jn);
      //!< assume the other effects are 0
    } break;
    case 3: //!< Model by Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
    {
      double OCVpt; //!< cathode potential
      try {
        OCVpt = OCV_curves.OCV_pos.interp(zp_surf, print);
        //!< get OCV of positive electrode, throw error if out of bounds
        //!< this should be updated for the cell's temperature using the entropic coefficient of the cathode
        //!< but I couldn't find any data on this, so I have ignored the effect
      } catch (int e) {
        //!< std::cout << "Throw test: " << 40 << '\n';
        std::cout << "Error in Cell_SPM::LAM when calculating the cathode potential for LAM: "
                  << e << ".\n";
        throw e;
      }

      //!< overpotential for the NMC dissolution reaction
      const double etap_LAM = OCVpt + etap - OCVnmc; //!< equation (9) from the paper

      //!< temperature dependent rate constant
      const double kt = lam_p.lam3k * exp(lam_p.lam3k_T * ArrheniusCoeff); //!< Arrhenius law

      //!< current density of the NMC dissolution reaction
      const double idiss = std::max(-5e-6, -kt * exp(n * F / Rg / st.T() * etap_LAM) / (n * F)); //!< equation (8) from the paper
      //!< cap the effect at 5e-6 to avoid a very fast drop of capacity (which could  cause an error)
      //!< a value of 5e-6 gives dap = -3.5. The initial value is about 17000, so the cell is dead in 5,000 seconds
      //!< so this cap is quite high

      //!< Add the effects of this model
      depp += idiss;
      //!< assume the other effects are 0
    } break;
    case 4: //!< Model by Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
      //!< Add the effects of this model
      dapp += -lam_p.lam4p * st.ap(); //!< equation (7) from the paper
      dann += -lam_p.lam4n * st.an();
      //!< assume the other effects are 0
      break;
    default: //!< unknown degradation model
      std::cerr << "ERROR in Cell_SPM::LAM, unknown LAM degradation model with identifier "
                << lam_id << ". Only values 0 to 4 are allowed. Throw an error.\n";
      throw 106;
      break;
    }
  }

  //!< Make the output variables
  *dthickp = dthickpp;
  *dthickn = dthicknn;
  *dap = dapp;
  *dan = dann;
  *dep = depp;
  *den = denn;
}

double Cell_SPM::LiPlating(double OCVnt, double etan)
{
  /*
   * Function to simulate the effect of lithium plating
   *
   * IN
   * OCVnt 	the OCV of the negative electrode at the current battery temperature [V]
   * etan 	the overpotential at the negative electrode [V]
   *
   * OUT
   * ipl 		current density for the plating side-reaction [A m-2]
   *
   * THROWS
   * 106 		illegal value in id
   */
  using namespace PhyConst;
  using std::exp;

  //!< Arrhenius relation for temperature-dependent plating parameters
  const double kplt = pl_p.pl1k * std::exp(pl_p.pl1k_T * calcArrheniusCoeff()); //!< Rate constant

  switch (deg_id.pl_id) {
  case 0:
    return 0;
  case 1: //!< Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
  {
    const auto temporary_var = (OCVnt + etan - OCVpl + rsei * st.delta() * I());
    return npl * F * kplt * exp(-n * F / (Rg * T()) * alphapl * temporary_var);
  }
  default:
    std::cerr << "ERROR in Cell_SPM::LiPlating, illegal degradation model identifier "
              << deg_id.pl_id << ", only values 0 and 1 are allowed. Throwing an error.\n";
    throw 106;
  }
}

void Cell_SPM::getDaiStress(double *sigma_p, double *sigma_n, sigma_type &sigma_r_p, sigma_type &sigma_r_n, sigma_type &sigma_t_p, sigma_type &sigma_t_n, sigma_type &sigma_h_p, sigma_type &sigma_h_n) noexcept
{
  /*
   * Calculates the radial and tangential stress for each positive Chebyshev node according to the formula by
   * Dai, Cai, White, Journal of Power sources 247, 2014
   *
   * It takes quite long to calculate the stress, so only call this function when needed.
   *
   * OUT
   * sigma_p		maximum hydrostatic stress in the positive particle, can be both positive and negative [Pa]
   * sigma_n		maximum hydrostatic stress in the negative particle, can be both positive and negative [Pa]
   * sigma_r_p	array with the radial stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
   * sigma_r_n	array with the radial stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
   * sigma_t_p	array with the tangential stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
   * sigma_t_n	array with the tangential stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
   * sigma_h_p	array with the hydrostatic stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
   * sigma_h_n	array with the hydrostatic stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
   * 				[0]			stress at the surface of the sphere
   * 				[1 to nch]	stress at the inner nodes
   * 				[nch + 1]	stress at the centre of the sphere
   */

  using settings::nch, std::pow;

  //!< Get the locations of the Chebyshev nodes
  double xp[nch + 2]; //!< location (x-value) of the positive Chebyshev nodes

  xp[0] = 1;
  std::copy(M->xch.begin(), M->xch.end(), &xp[1]);
  xp[nch + 1] = 0;

  double xtot[2 * nch + 3];    //!< location (x-value) of the positive and negative Chebyshev nodes [-surface .. centre .. +surface]
  xtot[nch + 1] = xp[nch + 1]; //!< centre node

  for (size_t i = 0; i < nch + 1; i++) {
    xtot[i] = -xp[i];                //!< negative nodes
    xtot[nch + 2 + i] = xp[nch - i]; //!< positive nodes
  }

  //!< get concentrations at each Chebyshev node
  //!< Due to symmetry, the concentration at the negative point is the same as the concentration of the positive point: c(-x) = c(x)
  double cp[nch + 2], cn[nch + 2]; //!< positive and negative nodes, [+surface .. inner .. centre]
  getC(cp, cn);

  double CP[2 * nch + 3], CN[2 * nch + 3]; //!< concentrations at all nodes, [-surface .. inner .. centre .. inner .. +surface]

  CP[nch + 1] = cp[nch + 1]; //!< cathode centre node
  CN[nch + 1] = cn[nch + 1]; //!< anode centre node

  for (size_t i = 0; i < nch + 1; i++) {
    CP[i] = cp[i];                 //!< cathode negative points
    CN[i] = cn[i];                 //!< anode negative points
    CP[nch + 2 + i] = cp[nch - i]; //!< cathode positive points
    CN[nch + 2 + i] = cn[nch - i]; //!< anode positive points
  }

  //!< The formula's to calculate the stress have integrals.
  //!< Integrals of Chebyshev points can be calculated using the Q-matrix in the state space struct (M)
  //!< The integral from the negative surface to node i is given by row i of the product Q*f
  //!< 		with Q the Chebyshev integration matrix
  //!< 			 f the value of the function you want to integrate, evaluated at every node
  //!< All integrals have to start from the negative surface (because you must cover the entire Chebyshev domain)
  //!< 	so if you need the integral of a function f from the centre (x = 0) to a positive point in the sphere (x = i)
  //!< 		int(f, x = 0 .. i) = int(f, x=-1 .. i) - int(f, x=-1 .. 0)
  //!< E.g. to get the integral of the positive li-concentration from the centre until the 4th positive Chebyshev node:
  //!< 		F = Q * CP 				[-surface .. centre .. +surface]
  //!< 		int(c(x), x = 0 .. i(4)) = int(c(x), x=-1 .. i(4)) - int(c(x), x=-1 .. 0)
  //!< 							  	 = F[nch + 1 + 4] 		   - F[nch+1]
  //!< 		(remember that the centre node is at [nch+1])

  //!< Calculate the matrix-vector product of the required functions as given in the paper by Dai et al. (concentration * radius^2)
  //!< we need to remember the transformation of variables from x (-1<x<1) to r (-R<r<R)
  //!< 		int( c * r^2 dr) = int(c * (x*R)^2 * d(R*x)) = int(c x^2 R^3 dx)
  //!< so the matrix-vector product we need is F = Q * (concentration * x^2 * R^3)

  //!< Note: since Fp (Fn) is multiplied by Rp^3 (Rn^3) then divided by them they are eliminated to remove numerical problems.

  std::array<double, 2 * nch + 3> Fp{}, Fn{}; //!< arrays with the product for the positive and negative electrode
  for (size_t i = 0; i < 2 * nch + 3; i++)    //!< loop for each row (one row = one node)
  {
    //!< calculate the matrix-vector product for row i as you would do it by hand:
    //!< F(i) = sum(Q(i,j)*C(j)*x(j)^2*R^3, j=0..2*nch+2)
    for (size_t j = 0; j < 2 * nch + 3; j++) {    //!< loop through the columns to calculate the sum
      Fp[i] += M->Q[i][j] * CP[j] * sqr(xtot[j]); //!< 		Q(i,j)*C(j)*x(j)^2
      Fn[i] += M->Q[i][j] * CN[j] * sqr(xtot[j]);
    }
  }

  //!< Calculate the integral from the centre to the positive surface, which is a constant present in all equations
  const double ap = Fp[2 * nch + 2] - Fp[nch + 1]; //!< int( cp*r^2, r=0..Rp ) // Note r^3 is eliminated!
  const double an = Fn[2 * nch + 2] - Fn[nch + 1]; //!< int( cn*r^2, r=0..Rn )

  *sigma_p = *sigma_n = 0;                      //!< Reset the stress values.
  for (size_t i = 0; i < sigma_h_p.size(); i++) //!< loop for the positive nodes
  {
    //!< r(i) = R * x(i) radius of positive node i in the positive particle
    //!< radius of positive node i in the negative particle
    const auto x_cube = cube(xtot[nch + 1 + i]);                  // Temporary variable
    const double bp_x = (Fp[nch + 1 + i] - Fp[nch + 1]) / x_cube; //!< integral from the centre to positive node i int(cp*zp^2, zp=0..rp(i)) = F[nch+1+i] - F[nch+1]
    const double bn_x = (Fn[nch + 1 + i] - Fn[nch + 1]) / x_cube; //!< integral from the centre to positive node i int(cn*zn^2, zn=0..rn(i))

    //!< Implement the equations from Dai et al.
    //!< Flip all arrays to get the opposite order (now it is [centre .. +surface] and we want [+surface .. centre]
    //!< and store in the output arrays
    const auto i_rev = sigma_h_p.size() - 1 - i; // reverse index.

    //!< centre node -> special formula (31 & 33) in Dai, Cai, White
    if (i == 0) {
      sigma_r_p[i_rev] = 2 * sparam.omegap * sparam.Ep / (9 * (1 - sparam.nup)) * (3 * ap - CP[nch + 1]);
      sigma_r_n[i_rev] = 2 * sparam.omegan * sparam.En / (9 * (1 - sparam.nun)) * (3 * an - CN[nch + 1]);

      sigma_t_p[i_rev] = 2 * sparam.omegap * sparam.Ep / (9 * (1 - sparam.nup)) * (3 * ap - CP[nch + 1]);
      sigma_t_n[i_rev] = 2 * sparam.omegan * sparam.En / (9 * (1 - sparam.nun)) * (3 * an - CN[nch + 1]);
    } else {                                                                                   //!< other nodes -> equation 13 in Dai, Cai, White
      sigma_r_p[i_rev] = 2 * sparam.omegap * sparam.Ep / (3 * (1 - sparam.nup)) * (ap - bp_x); //!< ap = int (c x^2, x=0..R), bp = int (c x^2 , x=0..r)
      sigma_r_n[i_rev] = 2 * sparam.omegan * sparam.En / (3 * (1 - sparam.nun)) * (an - bn_x);

      sigma_t_p[i_rev] = sparam.omegap * sparam.Ep / (3 * (1 - sparam.nup)) * (2 * ap + bp_x - cp[i]);
      sigma_t_n[i_rev] = sparam.omegan * sparam.En / (3 * (1 - sparam.nun)) * (2 * an + bn_x - cn[i]);
    }

    //!< Make the hydrostatic stress sh = (sr + 2sp)/3
    sigma_h_p[i_rev] = (sigma_r_p[i_rev] + 2 * sigma_t_p[i_rev]) / 3; //!< calculate hydrostatic stress
    sigma_h_n[i_rev] = (sigma_r_n[i_rev] + 2 * sigma_t_n[i_rev]) / 3;

    //!< find the maximum (in absolute value) of the stresses
    if (std::abs(sigma_h_p[i_rev]) > std::abs(*sigma_p))
      *sigma_p = sigma_h_p[i_rev]; //!< node with the maximum hydrostatic stress in the positive/negative particle
    if (std::abs(sigma_h_n[i_rev]) > std::abs(*sigma_n))
      *sigma_n = sigma_h_n[i_rev];
  }
}

void Cell_SPM::updateDaiStress() noexcept
{
  /*
   * Function which will update the values stored in the stress variables relating with Dai's stress model
   */

  //!< Make variables to store the stress
  sigma_type sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n;

  //!< Get the stress // #TODO these variables are not used.
  getDaiStress(&sparam.s_dai_p, &sparam.s_dai_n, sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n);
  //!< indicate that the values in the class variables are updated
  sparam.s_dai_update = true;
}

void Cell_SPM::getLaresgoitiStress(bool print, double *sigma_n)
{
  /*
   * Calculate the stress according to Laresgoiti's stress model
   *
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   *
   * OUT
   * sigma_n	stress in the negative particle [MPa]
   *
   * THROWS
   * 101 		the surface concentration is out of bounds
   */

  //!< Arrays with the stress from Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015 (figure 5)
  constexpr std::array<double, 11> xx{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };        //!< li-fraction in the graphite
  constexpr std::array<double, 11> yy{ 0.0, 5.5, 8.5, 9.5, 10.0, 10.0, 10.5, 13.0, 16.5, 21.0, 23.5 }; //!< stress [MPa]

  constexpr bool is_xx_fixed = true; //!< Change this if xx is not fixed time step.

  //!< Get the surface concentration
  double cps, cns;
  getCSurf(cps, cns, print);              //!< get the surface lithium concentration [mol m-3]
  const double zn_surf = (cns / Cmaxneg); //!< lithium fraction on negative surface [0, 1]

  //!< check if the surface concentration is within the allowed range
  //!< 	0 < cp < Cmaxpos
  //!< 	0 < cn < Cmaxneg
  if (cps < 0 || cns < 0 || cps > Cmaxpos || cns > Cmaxneg) {
    if (print) {
      std::cerr << "ERROR in Cell_SPM::getLaresgoitiStress: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos
                << " and the negative lithium fraction is " << cns / Cmaxneg << "they should both be between 0 and 1.\n";
    }
    throw 101;
  }

  //!< Interpolate linearly to get the stress, Make the output variable
  *sigma_n = linInt(print, true, xx, yy, 11, zn_surf, is_xx_fixed);
}

void Cell_SPM::updateLaresgoitiStress(bool print) // #TODO get and update seems unnecesary.
{
  /*
   * Function which will update the values stored in the stress variables relating with Laresgoiti's stress model
   *
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   */

  double s;
  getLaresgoitiStress(print, &s);

  //!< Update the stored value
  sparam.s_lares_n = s;
  sparam.s_lares_update = true; //!< indicate that the values in the class variables are updated
}

} // namespace slide