/**
 * @file Module_p.cpp
 * @brief Implementation of parallel module class.
 * @author Jorn Reniers, Volkan Kumtepeli
 * @date 18 Dec 2019
 */

#include "Module_p.hpp"

#include "../settings/settings.hpp"
#include "../utility/utility.hpp"
#include "../cells/cells.hpp"

#include <Eigen/Dense>

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <array>
#include <algorithm>
#include <ctime>

namespace slide {

/**
 * @brief Calculate the total resistance of the module.
 *
 * Computes the total resistance using the formula for resistances in parallel:
 * 1/Rtot = sum(1/R_i)
 *
 * @return The total resistance of the module.
 */
double Module_p::getRtot() // #TODO -> This function seems to be very expensive.
{
  if (SUs.empty()) return 0; // If there are no cells connected, return 0

  // Start from the cell furthest away
  double rtot = Rcontact.back() + SUs.back()->getRtot();

  // Then iteratively come closer, updating resistance  Rcontact[i] + (Rcell[i] \\ Rtot)
  for (int i = static_cast<int>(SUs.size()) - 2; i >= 0; i--) // #TODO Also check if Rcontact.empty() or in future it will be array.
    rtot = Rcontact[i] + (SUs[i]->getRtot() * rtot) / (SUs[i]->getRtot() + rtot);

  return rtot;
}

/**
 * @brief Calculate the voltage of each SU accounting for contact resistance.
 *
 * Computes the terminal voltage Vt for each SU in the module while considering
 * contact resistances and stores the values in the Vall span.
 *
 * @param Vall A span to store the voltage values.
 * @param print Unused parameter (can be removed if not required).
 */
void Module_p::getVall(std::span<double> Vall, bool print) // #TODO span may not be the best container here.
{
  /*
   * Return the voltage of SU[i] as seen from the terminal while accounting for the contact resistance
   *
   * Contact resistances:
   * 		we assume the terminals of the parallel module are on either side (i.e. next to SUs[0] and SUs[N-1])
   * 		the contact resistances are in the 'horizontal' paths, and the current of all subsequent cells goes through it
   * 			i.e. the current through Rc[0] = i[0] + i[1] + ... + i[N-1}
   * 		The series-resistance of every cell is already included in the cells themselves (as Rdc), and in the cell voltage v[i]
   * 		The terminal voltage Vt must be the same for the paths to all cells, i.e.
   * 			Vt = v[0] - R[0] (I[0] + I[1] + ... + I[N-1]) = v[0] - R[0]*I[0:N-1]
   * 			   = V[1] - R[1]*I[1:N-1] - R[0]*I[0:N-1]
   * 			   = v[i] - sum{ R[j]*I[j..N-1], j=0..i }
   * 			   = v[i] - sum{ R[j]*sum(I[k], k=j..N-1), j=0..i }
   */

  if (SUs.empty()) return;

  if (Vall.size() < SUs.size()) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module::getVall, container is too small!\n";
    throw 10;
  }

  double I_cumulative{ 0 };
  for (size_t i{}; i < SUs.size(); i++) {
    const auto j = SUs.size() - 1 - i; // Inverse indexing.
    Vall[j] = SUs[j]->V();             // Store the voltage of the current SU

    I_cumulative += SUs[j]->I(); // Update cumulative current

    for (auto k{ j }; k < SUs.size(); k++) // Update the voltage values in Vall considering contact resistances
      Vall[k] -= I_cumulative * Rcontact[j];
  }
}

/**
 * @brief Redistribute the current among the SUs to balance the voltage.
 *
 * Iteratively adjusts the current of each SU to balance the voltages across them.
 *
 * @param checkV Unused parameter (can be removed if not required).
 * @param print Unused parameter (can be removed if not required).
 * @return The status of the operation.
 */
Status Module_p::redistributeCurrent(bool checkV, bool print)
{
  // New redistributeCurrent without PI control:
  //!< voltage and initial current of each cell //!< #TODO it is a constant value SU.
  if (getNSUs() <= 1) return Status::Success;

  return setCurrent(I(), true, true);
}

Status Module_p::setVoltage(double Vnew, bool checkI, bool print)
{
  constexpr int maxIteration = 50;
  const auto nSU = getNSUs();
  int iter{}; // Current iteration


  auto StatusNow = Status::Success;

  using A_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, settings::MODULE_NSUs_MAX, settings::MODULE_NSUs_MAX>;
  using b_type = Eigen::Array<double, Eigen::Dynamic, 1, 0, settings::MODULE_NSUs_MAX>;

  b_type b(nSU), Iolds(nSU), Ib(nSU), Va(nSU), Vb(nSU), r_est(nSU);

  const double tolerance = 1e-9;

  double Icumulative{}, error_voltage{}, error_current{};
  for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
    Icumulative += Ib[i] = SUs[i]->I();
    Vb[i] = SUs[i]->V();

    if (i != 0)
      error_voltage += std::abs(Vb[i] - Vb[i - 1] - Icumulative * Rcontact[i]);
    else
      error_voltage += std::abs(Vb[i] - Vnew - Icumulative * Rcontact[i]);
  }

  if (error_voltage < tolerance) return StatusNow;

  Eigen::PartialPivLU<A_type> LU = [&]() { // #TODO this will be static in future defined in parallel block
    // Perturb a bit:
    const double perturbation = 0.5; //
    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      Ib[i] += perturbation; // Perturbation
      SUs[i]->setCurrent(Ib[i]);

      const double Vnew = SUs[i]->V();
      r_est[i] = (Vb[i] - Vnew) / perturbation; // Estimated resistance.
      Vb[i] = Vnew;
    }

    A_type A(nSU, nSU);
    // Set up the A matrix in Ax = b
    for (size_t j = 0; j < nSU; j++)
      for (size_t i = 0; i < nSU; i++) {
        if (i == j)
          A(i, j) = -Rcontact[i] - r_est[i];
        else if (i == j + 1)
          A(i, j) = r_est[j];
        else if (j > i)
          A(i, j) = -Rcontact[i];
        else
          A(i, j) = 0;
      }

    Eigen::PartialPivLU<A_type> LU(A); // LU decomposition of A.
    return LU;
  }();

  while (iter < maxIteration) {
    double Icumulative{}, error{};
    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      Icumulative += Ib[i];

      if (i != 0)
        b(i) = Vb[i] - Vb[i - 1] - Icumulative * Rcontact[i];
      else
        b(i) = Vb[i] - Vnew - Icumulative * Rcontact[i];

      error += std::abs(b(i));
    }

    if (error < tolerance) break;
    iter++;

    b_type deltaI = LU.solve(b.matrix()).array();

    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      StatusNow = std::max(StatusNow, SUs[i]->setCurrent(Ib[i] - deltaI[i]));
      Ib[i] = SUs[i]->I();
      Vb[i] = SUs[i]->V();
    }
  }

  if (iter == maxIteration)
    StatusNow = Status::RedistributeCurrent_failed;

  if constexpr (settings::printNumIterations)
    if (iter > 3) std::cout << "setVoltage iterations: " << iter << '\n';

  return StatusNow; // #TODO add some voltage/current etc. Also return max iter condition!!!!
}

Status Module_p::setCurrent(double Inew, bool checkV, bool print)
{
  /*
   * Set the current of a parallel module
   * This function takes small steps adapting the current of each connected cell until the total current is reached
   * 	 step 1 	change I by a bit  -> measure V_cell -> derive Rtot
   * 	 step 2 	do 50% of Inew-I() by increasing current proportionally to this resistance
   * 	 step 3   	iteratively change I of the cell with the smallest V (charge) or biggest V (discharge)
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 15 	after setting the current, the voltage of the cells are too far apart
   */

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  constexpr int maxIteration = 50;
  int iter{}; // Current iteration
  const int nSU = getNSUs();

  auto StatusNow = Status::Success;

  using A_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, settings::MODULE_NSUs_MAX, settings::MODULE_NSUs_MAX>;
  using b_type = Eigen::Array<double, Eigen::Dynamic, 1, 0, settings::MODULE_NSUs_MAX>;


  b_type b(nSU), Iolds(nSU), Ib(nSU), Va(nSU), Vb(nSU), r_est(nSU);


  StatusNow = Status::Success; //!< reset at each iteration.

  const double tolerance = 1e-6;
  for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
    Ib[i] = SUs[i]->I();
    Vb[i] = SUs[i]->V();
  }


  double Icumulative{}, error{};
  for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
    Icumulative += Ib[i];

    if (i != 0)
      b(i) = Vb[i] - Vb[i - 1] - Icumulative * Rcontact[i];
    else
      b(i) = Inew - Icumulative; // #TODO????

    error += std::max(std::abs(b(i)), error);

    r_est[i] = 200e-3; // 100e-3; // Init the resistances high at the beginning. #TODO probably not robust for all cases.
  }

  if (error < tolerance) return StatusNow;


  // // Perturb a bit:
  // const double perturbation = 0.5; //
  // for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
  //   Ib[i] += perturbation; // Perturbation
  //   SUs[i]->setCurrent(Ib[i]);

  //   const double Vnew = SUs[i]->V();
  //   r_est[i] = (Vb[i] - Vnew) / perturbation; // Estimated resistance.
  //   Vb[i] = Vnew;
  // }

  auto getLU = [&]() { // #TODO this will be static in future defined in parallel block
    A_type A(nSU, nSU);
    // Set up the A matrix in Ax = b
    for (size_t j = 0; j < nSU; j++)
      for (size_t i = 0; i < nSU; i++) {
        if (i == 0)
          A(i, j) = -1;
        else if (i == j)
          A(i, j) = -Rcontact[i] - r_est[i];
        else if (i == j + 1)
          A(i, j) = r_est[j];
        else if (j > i)
          A(i, j) = -Rcontact[i];
        else
          A(i, j) = 0;
      }
    Eigen::PartialPivLU<A_type> LU(A); // LU decomposition of A.
    return LU;
  };

  b_type deltaI;
  Eigen::PartialPivLU<A_type> LU;
  while (iter < maxIteration) {
    double Icumulative{}, error{};
    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      Icumulative += Ib[i];

      if (i != 0)
        b(i) = Vb[i] - Vb[i - 1] - Icumulative * Rcontact[i];
      else
        b(i) = Inew - Icumulative; // #TODO????

      error = std::max(std::abs(b(i)), error);
    }

    if (error < tolerance) break;

    LU = getLU();

    deltaI = LU.solve(b.matrix()).array();


    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      double Vb_old = Vb[i];
      StatusNow = std::max(StatusNow, SUs[i]->setCurrent(Ib[i] - deltaI[i]));
      Ib[i] = SUs[i]->I();
      Vb[i] = SUs[i]->V();


      if (std::abs(deltaI[i]) > 1e-6 || iter == 0) {
        double r_est_i = (Vb[i] - Vb_old) / deltaI[i];

        if (r_est_i > 0)
          r_est[i] = r_est_i;
        else
          r_est[i] = 1e-3;
      }
    }

    iter++;
  }

  if (iter == maxIteration)
    StatusNow = Status::RedistributeCurrent_failed;

  if constexpr (settings::printNumIterations)
    if (iter > 5) std::cout << "setCurrent iterations: " << iter << std::endl;


  // #TODO set old currents back here!
  return StatusNow;
}


/**
 * @brief Perform a time step at a constant current.
 *
 * Takes a time step at a constant current by either explicitly solving the system of equations
 * or letting every cell take a CC time step and checking if the voltage equation is satisfied.
 *
 * @param dt The time step.
 * @param nstep Number of steps.
 */
void Module_p::timeStep_CC(double dt, int nstep)
{
  assert(dt > 0); // Check if the time step is valid (only in debug mode)

  //!< we simply take one CC time step on every cell
  auto task_indv = [&](int i) { SUs[i]->timeStep_CC(dt, nstep); };

  try {
    run(task_indv, getNSUs(), (par ? -1 : 1));
  } catch (int e) {
    std::cout << "Error in Module_p::timeStep_CC with module ID " << getFullID()
              << ". error " << e << ", throwing it on.\n";
    throw e;
  }

  //!< **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************

  //!< update the time since the last update of the thermal model
  if (!blockDegAndTherm) {
    therm.time += nstep * dt;

    //!< Increase the heat from the contact resistances
    double Ii = 0; //!< current through resistor I
    for (size_t i = 0; i < SUs.size(); i++) {
      for (size_t j = i; j < SUs.size(); j++) // #TODO very important! Not calculating current for Rcontact properly!!!!!!!!!!
        Ii += SUs[j]->I();                    //!< resistor i sees the currents through the cells 'behind' them

      therm.Qcontact += Rcontact[i] * sqr(Ii) * nstep * dt;
    }

    //!< If this module has a parent module, this parent will call the thermal model with the correct parameters
    //!< which will include heat exchange with the module's neighbours and cooling from the cooling system of the parent module.

    //!< if there is no parent, this module is the top-level.
    //!< It then directly exchanges heat with the environment at a fixed temperature
    //!< If there is no parent, assume we update T every nstep*dt. So update the temperature now
    if (!parent) {
      //!< double check that we have an HVAC coolsystem (see constructor)
      if (typeid(*getCoolSystem()) != typeid(CoolSystem_HVAC)) { // #TODO expression with side effects will be evaluated despite being used as an operand to 'typeid' [-Wpotentially-evaluated-expression]
        std::cerr << "ERROR in module_p::timeStep_CC in module " << getFullID() << ". this is a top-level"
                  << "module but does not have an HVAC coolsystem for active cooling with the environment.\n";
        throw 14;
      }

      //!< Call the thermal model without heat exchanges with neighbours or parents (since this module doesn't have any)
      //!< 	Heat exchange between this module and the environment is done by the AC-part of the HVAC coolsystem
      //!< 	heat exchange between this module and its children is done the conventional way by the HVAC coolsystem
      //!< set the new temperature since we have calculated all the temperatures so there is no risk for inconsistency
      double Tneigh[1], Kneigh[1], Aneigh[1];
      //!< Make the arrays even though they will not be used (length should be 0 but I don't think you can make an array of length 0)
      setT(thermalModel(0, Tneigh, Kneigh, Aneigh, therm.time)); //!< the 0 signals there are no neighbours or parents

      /*
      double Tneigh[1] = {settings::T_ENV};							//!< T of environment
      double Kneigh[1] = {cool->getH()}; 					//!< h to environment is same as to children, since the speed of the coolant is the same
      double Aneigh[1] = {therm.A};						//!< A to the environment is the A of this module
      //!< note that this will have more heat exchange than to children, since A = min(A_this, A_parent)
      //!< 	children will have a smaller A themselves, so the resulting A = A_child < A of this module
      setT(thermalModel(1, Tneigh, Kneigh, Aneigh, therm.time));*/
    }

    //!< control the cooling system
    double Tlocal = 0;
    for (auto &SU : SUs)
      Tlocal = std::max(Tlocal, SU->T());

    cool->control(Tlocal, getThotSpot());
  }
}

} // namespace slide