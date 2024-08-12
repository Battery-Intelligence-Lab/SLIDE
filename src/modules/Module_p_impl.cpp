/**
 * @file Module_p_impl.cpp
 * @brief Parallel branch solution implementations
 * @author Nilsu Atlan
 * @author Volkan Kumtepeli
 * @date 12 Agus 2024
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

typedef std::vector<double> state_type;

double ocv_eval(const std::vector<double> &ocv_coefs, double z)
{
  double ocv = 0.0;
  for (size_t i = 0; i < ocv_coefs.size(); ++i) {
    ocv += ocv_coefs[i] * std::pow(z, ocv_coefs.size() - 1 - i);
  }
  return ocv;
} // already existing function --> otomatik olarak değer fonksiyona gelecek


void calculate_inverse_A22(int n_par, Eigen::VectorXd r, Eigen::MatrixXd &m)
{
  // Initialize m as a zero matrix
  m = Eigen::MatrixXd::Zero(n_par, n_par);

  // Calculate the sum of 1/r_l for l from 1 to n
  double Rsum_inv = 0.0;
  for (int l = 0; l < n_par; ++l) {
    Rsum_inv += 1.0 / r(l);
  }

  // Calculate m_{k,k} for k = 1, ..., n
  for (int k = 0; k < n_par; ++k) {
    m(k, k) = 1.0 / r(k) * Rsum_inv;
  }

  // Calculate m_{j,k} for j = 1, ..., n-1 and j ≠ k
  for (int j = 0; j < n_par - 1; ++j) {
    for (int k = 0; k < n_par; ++k) {
      if (j != k) {
        m(j, k) = (r(k) / r(j)) * m(j, j);
      }
    }
  }

  // Calculate m_{k,n} for k = 1, ..., n
  for (int k = 0; k < n_par; ++k) {
    m(k, n_par - 1) = 1.0 / r(k) * (1.0 / Rsum_inv);
  }
}

void new_compute_A11_A12_A21_A22(int n_par, Eigen::VectorXd R, Eigen::VectorXd C, Eigen::VectorXd Q, Eigen::VectorXd tau, Eigen::VectorXd r, Eigen::MatrixXd &A11, Eigen::MatrixXd &A12, Eigen::MatrixXd &A21, Eigen::MatrixXd &A22, Eigen::MatrixXd &m)
{
  // Initialize A11 and A12
  A11 = Eigen::MatrixXd::Zero(2 * n_par, 2 * n_par);
  A12 = Eigen::MatrixXd::Zero(2 * n_par, n_par);

  A11.block(0, 0, 2, 2) << 0, 0, 0, -tau(0);
  // A12.block(0, 0, 2, 1) << 1 / Q(0), 0, 0, 1 / C(0);
  A12(0, 0) = 1 / Q(0);
  A12(1, 0) = 1 / C(0);

  for (int i = 1; i < n_par; ++i) {
    A11.block(2 * i, 2 * i, 2, 2) << 0, 0, 0, -tau(i);
    // A12.block(2 * i, i, 2, 1) << 1 / Q(i), 0, 0, 1 / C(i);
    A12(2 * i, i) = 1 / Q(i);
    A12(2 * i + 1, i) = 1 / C(i);
  }

  // Initialize A22
  A22 = Eigen::MatrixXd::Zero(n_par, n_par);
  Eigen::VectorXd main_diag = -(r + R);
  main_diag(0) = r(0);

  Eigen::VectorXd lower_diag = r.head(n_par - 1);
  lower_diag(0) = r(0);

  A22.diagonal() = main_diag;
  A22.diagonal(-1) = lower_diag;

  for (int i = 1; i < n_par; ++i) {
    for (int j = i + 1; j < n_par; ++j) {
      A22(i, j) = -R(i);
    }
  }

  A22.row(0).setOnes();

  // Initialize A21
  A21 = Eigen::MatrixXd::Zero((n_par - 1), n_par * 2);
  for (int i = 0; i < n_par - 1; ++i) {
    for (int j = 0; j < n_par; ++j) {
      A21(i, 2 * j + 1) = A22(i, j);
    }
  }

  // Calculate m (inverse of A22)
  /* m = Eigen::MatrixXd::Zero(n_par, n_par);
  double Rsum_inv = 1 / R.sum();

  for (int ktil = 0; ktil < n_par; ++ktil) {
    for (int i = 0; i < n_par - 1; ++i) {
      if (ktil - i - 1 == 0) {
        m(i + 1, i) = std::pow(1 / R(i + 1), 2) * Rsum_inv - 1 / R(i + 1);
      } else {
        m(ktil, i) = 1 / (R(ktil) * R(i + 1)) * Rsum_inv;
      }
    }
  }

  m(n_par - 1, n_par - 1) = 1 / (R(n_par - 1) * Rsum_inv);
  for (int i = 0; i < n_par - 1; ++i) {
    m(i, n_par - 1) = 1 / (R(i) * Rsum_inv);
  } */

  m = Eigen::MatrixXd::Zero(n_par, n_par);
  calculate_inverse_A22(n_par, r, m);

  Eigen::MatrixXd inv_A22 = A22.inverse();
  double error = (inv_A22 - m).norm();

  // std::cout << "Error: " << error << std::endl;
}

// Function to evaluate the vector field of the parallel pack model
Eigen::VectorXd parallel_model(double t, const Eigen::VectorXd &x, double I, const Eigen::MatrixXd &A11, const Eigen::MatrixXd &A12, const Eigen::MatrixXd &A21, const Eigen::MatrixXd &A22, const std::vector<double> &ocv_coefs)
{
  int n = x.size() / 2;
  Eigen::VectorXd SOC(n), OCV_branch(n), w(n), v_mod(n);

  for (int i = 0; i < n; ++i) {
    SOC(i) = x(2 * i);
    OCV_branch(i) = ocv_eval(ocv_coefs, SOC(i));
    w(i) = x(2 * i + 1);
    v_mod(i) = OCV_branch(i) + w(i);
  }

  Eigen::VectorXd delta_v(n - 1);
  for (int i = 1; i < n; ++i) {
    delta_v(i - 1) = v_mod(i) - v_mod(0);
  }

  Eigen::VectorXd phi(n);
  phi << -delta_v, -I;

  Eigen::VectorXd i_branch = -A22 * phi;
  Eigen::VectorXd xdot = A11 * x + A12 * (-A22 * phi);
  return xdot;
}

// Function to evaluate the vector field of the parallel pack model using precomputed matrix inverse m
Eigen::VectorXd parallel_model_no_mat_inverse(double t, const Eigen::VectorXd &x, double I, const Eigen::MatrixXd &A11, const Eigen::MatrixXd &A12, const Eigen::MatrixXd &A21, const Eigen::MatrixXd &A22, const std::vector<double> &ocv_coefs, const Eigen::MatrixXd &m)
{
  int n = x.size() / 2;
  Eigen::VectorXd SOC(n), OCV_branch(n), w(n), v_mod(n);

  for (int i = 0; i < n; ++i) {
    SOC(i) = x(2 * i);
    OCV_branch(i) = ocv_eval(ocv_coefs, SOC(i));
    w(i) = x(2 * i + 1);
    v_mod(i) = OCV_branch(i) + w(i);
  }

  Eigen::VectorXd delta_v(n - 1);
  for (int i = 1; i < n; ++i) {
    delta_v(i - 1) = v_mod(i) - v_mod(0);
  }

  Eigen::VectorXd phi(n);
  phi << -delta_v, -I;

  Eigen::VectorXd i_branch = -A22 * phi;
  Eigen::VectorXd xdot = A11 * x + A12 * (-m * phi);
  return xdot;
}

void parallel_model_dae5(double t, const state_type &x_std, double I, const Eigen::MatrixXd &A11, const Eigen::MatrixXd &A12, const Eigen::MatrixXd &A21, const Eigen::MatrixXd &A22, const std::vector<double> &ocv_coefs, const Eigen::VectorXd &R, const Eigen::VectorXd &r, state_type &xdot_std, Eigen::VectorXd &i_branch)
{

  constexpr bool printDebug = false;

  Eigen::Map<const Eigen::VectorXd> x(x_std.data(), x_std.size());
  Eigen::Map<Eigen::VectorXd> xdot(xdot_std.data(), xdot_std.size());

  if (printDebug)
    std::cout << "x : " << x.transpose() << std::endl;


  int n = x.size() / 2;
  static Eigen::VectorXd OCV_branch(n), w(n), v_mod(n);

  for (int i = 0; i < n; ++i) {
    const auto SOC = x(2 * i);
    OCV_branch(i) = ocv_eval(ocv_coefs, SOC);
    w(i) = x(2 * i + 1);
    v_mod(i) = OCV_branch(i) + w(i);
  }

  if (printDebug) {
    std::cout << "v_mod : " << v_mod.transpose() << std::endl;
    std::cout << "w : " << w.transpose() << std::endl;
    std::cout << "OCV_branch : " << OCV_branch.transpose() << std::endl;
  }


  // static Eigen::VectorXd i_branch(n)
  static Eigen::VectorXd theta(n), rho(n);
  theta(0) = 0;
  rho(0) = 0;

  for (int i = 1; i < n; ++i) {
    theta(i) = (R(i) + r(i)) / r(i - 1);
    rho(i) = 1 / r(i - 1);
  }

  if (printDebug) {
    std::cout << "theta : " << theta.transpose() << std::endl;
    std::cout << "rho : " << rho.transpose() << std::endl;
  }


  double cn = 1.0;
  for (int j = 1; j < n; ++j) {
    double prod_theta = 1.0;
    for (int k = j; k < n; ++k) {
      prod_theta *= theta(k);
    }
    cn += prod_theta;
  }

  if (printDebug)
    std::cout << "cn : " << cn << std::endl;


  static Eigen::VectorXd f(n);
  f.fill(0);

  for (int j = 0; j < n - 1; ++j) {
    f(j) = rho(j + 1) * (v_mod(j + 1) - v_mod(j));
    for (int i = j + 1; i < n - 1; ++i) {
      double prod_theta = 1.0;
      for (int k = i; k < n - 1; ++k) {
        prod_theta *= theta(k + 1);
      }
      f(j) += prod_theta * rho(i + 1) * (v_mod(i + 1) - v_mod(j));
    }
  }
  f(n - 1) = rho(n - 1) * (v_mod(n - 1) - v_mod(n - 2));

  if (printDebug)
    std::cout << "f : " << f.transpose() << std::endl;

  i_branch(n - 1) = (I - f.sum()) / cn;

  for (int j = n - 2; j >= 0; --j)
    i_branch(j) = theta(j + 1) * i_branch(j + 1) + rho(j + 1) * (v_mod(j + 1) - v_mod(j));

  xdot = A11 * x + A12 * i_branch;

  if (printDebug) {
    std::cout << "i_branch : " << i_branch.transpose() << std::endl;
    std::cout << "xdot : " << xdot.transpose() << std::endl;
  }
}

// NİLSU
Status Module_p::setCurrent_analytical_impl(double Inew, bool checkV, bool print)
{
  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold
  constexpr bool printDebug = false;

  constexpr int maxIteration = 50;
  int iter{}; // Current iteration
  const int nSU = getNSUs();

  Eigen::VectorXd Q(nSU);   // Battery Capacity in As
  Eigen::VectorXd F(nSU);   // RC devresinin bir elemani
  Eigen::VectorXd C(nSU);   // RC devresinin bir elemani
  Eigen::VectorXd r(nSU);   // Rdc.
  Eigen::VectorXd tau(nSU); // 1/(RC)

  Eigen::VectorXd R = Eigen::VectorXd::Constant(nSU, 0.0); // Rcontant

  Eigen::VectorXd i_branch(nSU);

  for (int i{}; i < nSU; i++) {
    auto cp = dynamic_cast<Cell_ECM<1> *>(SUs[i].get());

    Q(i) = cp->Cap() * 3600;
    r(i) = cp->getRtot();
    F(i) = cp->getRp(0); // 0 -> because we have one RC pair.
    C(i) = cp->getC(0);
    tau(i) = 1 / (F(i) * C(i));
  }

  //   std::cout << "Q:\n"
  //             << Q << "\n\n";
  //   std::cout << "r:\n"
  //             << r << "\n\n";
  //   std::cout << "R:\n"
  //             << R << "\n\n";
  //   std::cout << "C:\n"
  //             << C << "\n\n";
  //   std::cout << std::endl;

  /// ---------------------
  Eigen::MatrixXd A11, A12, A21, A22, m;
  new_compute_A11_A12_A21_A22(nSU, R, C, Q, tau, r, A11, A12, A21, A22, m);

  if (printDebug) {
    std::cout << "A22:\n"
              << A22 << std::endl;
    std::cout << "m:\n"
              << m << std::endl;
  }


  // ----------- parallel_model_dae5 -----
  static Eigen::VectorXd OCV_branch(nSU), w(nSU), v_mod(nSU);

  for (int i{}; i < nSU; i++) {
    auto cp = dynamic_cast<Cell_ECM<1> *>(SUs[i].get());
    const auto SOC = cp->SOC(); // [0-1]
    OCV_branch(i) = cp->getOCV();
    w(i) = cp->getVr(); // Relaxation voltage.
    v_mod(i) = OCV_branch(i) + w(i);
  }

  if (printDebug) {
    std::cout << "v_mod : " << v_mod.transpose() << std::endl;
    std::cout << "w : " << w.transpose() << std::endl;
    std::cout << "OCV_branch : " << OCV_branch.transpose() << std::endl;
  }

  // static Eigen::VectorXd i_branch(nSU)
  static Eigen::VectorXd theta(nSU), rho(nSU);
  theta(0) = 0;
  rho(0) = 0;

  for (int i = 1; i < nSU; ++i) {
    theta(i) = (R(i) + r(i)) / r(i - 1);
    rho(i) = 1 / r(i - 1);
  }

  if (printDebug) {
    std::cout << "theta : " << theta.transpose() << std::endl;
    std::cout << "rho : " << rho.transpose() << std::endl;
  }

  double cn = 1.0;
  for (int j = 1; j < nSU; ++j) {
    double prod_theta = 1.0;
    for (int k = j; k < nSU; ++k) {
      prod_theta *= theta(k);
    }
    cn += prod_theta;
  }

  if (printDebug)
    std::cout << "cn : " << cn << std::endl;


  static Eigen::VectorXd f(nSU);
  f.fill(0);

  for (int j = 0; j < nSU - 1; ++j) {
    f(j) = rho(j + 1) * (v_mod(j + 1) - v_mod(j));
    for (int i = j + 1; i < nSU - 1; ++i) {
      double prod_theta = 1.0;
      for (int k = i; k < nSU - 1; ++k) {
        prod_theta *= theta(k + 1);
      }
      f(j) += prod_theta * rho(i + 1) * (v_mod(i + 1) - v_mod(j));
    }
  }
  f(nSU - 1) = rho(nSU - 1) * (v_mod(nSU - 1) - v_mod(nSU - 2));

  if (printDebug)
    std::cout << "f : " << f.transpose() << std::endl;


  i_branch(nSU - 1) = (-Inew - f.sum()) / cn; // -Inew because discharge is positive in SLIDE.

  for (int j = nSU - 2; j >= 0; --j)
    i_branch(j) = theta(j + 1) * i_branch(j + 1) + rho(j + 1) * (v_mod(j + 1) - v_mod(j));


  auto StatusNow = Status::Success;

  Eigen::VectorXd Vb(nSU);
  std::vector<double> Vcontant_included(nSU);

  for (int i{}; i < nSU; i++) {
    auto cp = dynamic_cast<Cell_ECM<1> *>(SUs[i].get());
    StatusNow = std::max(StatusNow, cp->setCurrent(-i_branch(i)));

    Vb(i) = cp->V();
  }

  getVall(Vcontant_included, false); // Get voltages with contact resistances included.

  if (printDebug) {
    std::cout << "I:\n"
              << -i_branch << "\n";

    std::cout << "Vb:\n"
              << Vb << "\n"; // These should be equal if no contact resistance.

    std::cout << "Vcontant_included:\n";
    for (auto V_ : Vcontant_included)
      std::cout << V_ << '\n';

    std::cout << "\n";
  }


  return StatusNow;
}

Status Module_p::setCurrent_previous_impl(double Inew, bool checkV, bool print)
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
    Ib(i) = SUs[i]->I();
    Vb(i) = SUs[i]->V();
  }


  double Icumulative{}, error{};
  for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
    Icumulative += Ib[i];

    if (i != 0)
      b(i) = Vb(i) - Vb(i - 1) - Icumulative * Rcontact[i];
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
      Icumulative += Ib(i);

      if (i != 0)
        b(i) = Vb(i) - Vb(i - 1) - Icumulative * Rcontact[i];
      else
        b(i) = Inew - Icumulative; // #TODO????

      error = std::max(std::abs(b(i)), error);
    }

    if (error < tolerance) break;

    LU = getLU();

    deltaI = LU.solve(b.matrix()).array();


    for (size_t i = SUs.size() - 1; i < SUs.size(); i--) {
      double Vb_old = Vb(i);
      StatusNow = std::max(StatusNow, SUs[i]->setCurrent(Ib[i] - deltaI[i]));
      Ib(i) = SUs[i]->I();
      Vb(i) = SUs[i]->V();


      if (std::abs(deltaI[i]) > 1e-6 || iter == 0) {
        double r_est_i = (Vb(i) - Vb_old) / deltaI[i];

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

} // namespace slide