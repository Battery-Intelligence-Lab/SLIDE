/*
 * parallel_model_sln.cpp
 *
 * Example for integration test of the analytical parallel model solution.
 *
 * Created on: 01 Aug 2024
 * Author(s): Nilsu Atlan, Volkan Kumtepeli
 */


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <slide.hpp>
#include <chrono>
#include <boost/numeric/odeint.hpp>

// using namespace std; // bu yerine fonksiyonun içine ekle
// using namespace Eigen;

// Global variable to store i_branch values
Eigen::MatrixXd i_branch_store;

using namespace boost::numeric;

typedef Eigen::VectorXd state_type;

void harmonic_oscillator(const state_type &x, state_type &dxdt, double t)
{
  dxdt[0] = x[1];
  dxdt[1] = -x[0];
}

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

  std::cout << "Error: " << error << std::endl;
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
  i_branch_store = i_branch; // Store i_branch globally
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
  i_branch_store = i_branch; // Store i_branch globally
  Eigen::VectorXd xdot = A11 * x + A12 * (-m * phi);
  return xdot;
}

void parallel_model_dae5(double t, const Eigen::VectorXd &x, double I, const Eigen::MatrixXd &A11, const Eigen::MatrixXd &A12, const Eigen::MatrixXd &A21, const Eigen::MatrixXd &A22, const std::vector<double> &ocv_coefs, const Eigen::VectorXd &R, const Eigen::VectorXd &r, Eigen::VectorXd &xdot)
{
  std::cout << "x : " << x << std::endl;


  int n = x.size() / 2;
  Eigen::VectorXd SOC(n), OCV_branch(n), w(n), v_mod(n);

  for (int i = 0; i < n; ++i) {
    SOC(i) = x(2 * i);
    OCV_branch(i) = ocv_eval(ocv_coefs, SOC(i));
    w(i) = x(2 * i + 1);
    v_mod(i) = OCV_branch(i) + w(i);
  }

  std::cout << "SOC : " << SOC << std::endl;
  std::cout << "v_mod : " << v_mod << std::endl;
  std::cout << "w : " << w << std::endl;
  std::cout << "OCV_branch : " << OCV_branch << std::endl;


  Eigen::VectorXd i_branch(n), theta(n), rho(n);
  theta(0) = 0;
  rho(0) = 0;

  for (int i = 1; i < n; ++i) {
    theta(i) = (R(i) + r(i)) / r(i - 1);
    rho(i) = 1 / r(i - 1);
  }

  std::cout << "theta : " << theta << std::endl;
  std::cout << "rho : " << rho << std::endl;


  double cn = 1.0;
  for (int j = 1; j < n; ++j) {
    double prod_theta = 1.0;
    for (int k = j; k < n; ++k) {
      prod_theta *= theta(k);
    }
    cn += prod_theta;
  }

  std::cout << "cn : " << cn << std::endl;


  Eigen::VectorXd f = Eigen::VectorXd::Zero(n);

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
  std::cout << "f : " << f << std::endl;

  i_branch(n - 1) = (I - f.sum()) / cn;

  for (int j = n - 2; j >= 0; --j) {
    i_branch(j) = theta(j + 1) * i_branch(j + 1) + rho(j + 1) * (v_mod(j + 1) - v_mod(j));
  }

  i_branch_store = i_branch; // Store i_branch globally
  xdot = A11 * x + A12 * i_branch;

  std::cout << "i_branch : " << i_branch << std::endl;


  std::cout << "xdot : " << xdot << std::endl;
}


int main()
{
  int bat_select_num = 3;
  // Eigen::VectorXd R(n_par), C(n_par), Q(n_par), tau(n_par), r(n_par);

  // R << 0, 0, 0;
  // C << 634.0, 634.0, 634.0;
  // Q << 9000.0, 9000.0, 9000.0;
  // tau << 0.04, 0.04, 0.04;
  // r << 0.029, 0.029, 0.029;

  std::vector<double> capacitance_Ah_all = { 2.50, 3.0, 2.60 };
  std::vector<double> R1_all = { 0.0021, 0.0019, 0.0394 };
  std::vector<double> C1_all = { 23558, 23340, 634 };
  std::vector<double> r_all = { 0.0084, 0.0037, 0.0291 };
  std::vector<double> ocv_coefs = { -5147.78793933142, 28772.2793299200, -68881.8548465757 };
  // std::vector<double> ocv_coefs = {-5147.78793933142,	28772.2793299200,	-68881.8548465757}; // Example OCV coefficients

  // Open circuit voltage
  int nz = 100;
  std::vector<double> z_space(nz), OCV_plot(nz);
  for (int i = 0; i < nz; ++i) {
    z_space[i] = static_cast<double>(i) / (nz - 1);
    OCV_plot[i] = ocv_eval(ocv_coefs, z_space[i]);
  }

  // Model parameters
  // int n_par = 100; // The number of cells in parallel
  int n_par = 3; // for testing

  double capacitance_Ah = capacitance_Ah_all[bat_select_num - 1];
  double capacitance = 3600 * capacitance_Ah;

  double C_rate = 0.01;
  double I = n_par * C_rate * capacitance_Ah;

  double R_nom = 0 * 1e-5;

  double r_nom = r_all[bat_select_num - 1];
  double C_nom = C1_all[bat_select_num - 1];

  Eigen::VectorXd Q = Eigen::VectorXd::Constant(n_par, capacitance);
  Eigen::VectorXd R = Eigen::VectorXd::Constant(n_par, R_nom);
  Eigen::VectorXd C = Eigen::VectorXd::Constant(n_par, C_nom);
  Eigen::VectorXd r = Eigen::VectorXd::Constant(n_par, r_nom);

  // Add noise to each of the parameters to simulate aging
  double sd_vars = 1 * 1e-3;
  double var_Q = sd_vars * capacitance;
  double var_R = sd_vars * R_nom;
  double var_r = sd_vars * r_nom;
  double var_C = sd_vars * C_nom;

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);

  for (int j = 0; j < n_par; ++j) {
    Q(j) += var_Q * distribution(generator);
    R(j) += var_R * distribution(generator);
    r(j) += var_r * distribution(generator);
    C(j) += var_C * distribution(generator);
  }

  R_nom = R1_all[bat_select_num - 1];
  Eigen::VectorXd F = Eigen::VectorXd::Constant(n_par, R_nom);
  Eigen::VectorXd tau = F.cwiseProduct(C).cwiseInverse();

  Eigen::MatrixXd A11, A12, A21, A22, m;
  new_compute_A11_A12_A21_A22(n_par, R, C, Q, tau, r, A11, A12, A21, A22, m);

  // std::cout << "A11:\n" << A11 << std::endl;
  // std::cout << "A12:\n" << A12 << std::endl;
  // std::cout << "A21:\n" << A21 << std::endl;
  std::cout << "A22:\n"
            << A22 << std::endl;
  std::cout << "m:\n"
            << m << std::endl;
  // std::cout << "tau:\n" << tau << std::endl;
  // std::cout << "I:\n" << I << std::endl;
  // std::cout << "x0:\n" << x0 << std::endl;
  // std::cout << "tspan:\n" << tspan << std::endl;

  Eigen::VectorXd x0 = Eigen::VectorXd::Constant(2 * n_par, 0.2);
  for (int j = 0; j < n_par; ++j) {
    x0(2 * j + 1) = 0;
  }

  double t_f = 0.6 * 3600 / C_rate;
  int num_points = 100; // need to explicitly specify the number of points because there is no default value (MATLAB default = 100)
  Eigen::VectorXd tspan = Eigen::VectorXd::LinSpaced(num_points, 0, t_f);

  auto start = std::chrono::high_resolution_clock::now();
  // Call the functions
  // Eigen::VectorXd xdot_inv = parallel_model(t, x, I, A11, A12, A21, A22, ocv_coefs);
  // Eigen::VectorXd xdot_m = parallel_model_no_mat_inverse(t, x, I, A11, A12, A21, A22, ocv_coefs, m);
  // Eigen::VectorXd xdot_dae = parallel_model_dae5(t, x, I, A11, A12, A21, A22, ocv_coefs, m);

  // std::cout << "xdot_inv:\n" << xdot1 << std::endl;
  // std::cout << "xdot_m:\n" << xdot2 << std::endl;
  // std::cout << "xdot_dae:\n" << xdot2 << std::endl;

  // Use the Runge-Kutta Cash-Karp method (similar to MATLAB's ode45)
  odeint::runge_kutta_cash_karp54<state_type> stepper;

  // Integrate the ODE from t=0 to t=10 with step size 0.1
  double t_start = 0.0, t_end = 10.0, dt = 0.1;
  auto myODE = [&](auto &x, auto &dxdt, double t) {
    parallel_model_dae5(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r, dxdt);
  };

  odeint::integrate_adaptive(odeint::make_controlled(1.0e-12, 1.0e-12, stepper), myODE, x0, t_start, t_end, dt);


  // ...........
  // ODE SOLVER
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_dae = end - start;
  std::cout << "DAE5 computation took: " << time_dae.count() << " seconds" << std::endl;

  Eigen::MatrixXd i_branch_values = i_branch_store;
  std::cout << "i_branch_values:\n"
            << i_branch_values << std::endl;

  // x_ode = x_ode_2;
  // t_ode = t_ode_2;

  /*int size_t = x_ode.rows();
  Eigen::MatrixXd z_ode(size_t, n_par), w_ode(size_t, n_par);
  for (int i = 0; i < n_par; ++i) {
      z_ode.col(i) = x_ode.col(2 * i);
      w_ode.col(i) = x_ode.col(2 * i + 1);
  }

  int nT = tspan.size();
  Eigen::MatrixXd i_branch(n_par, nT);
  Eigen::VectorXd phi(n_par), phi_init(n_par - 1);

  Eigen::MatrixXd OCV_ode(nT, n_par), volts(nT, n_par);
  Eigen::VectorXd current_sum(nT);

  for (int j = 0; j < nT; ++j) {
      for (int i = 0; i < n_par; ++i) {
          OCV_ode(j, i) = ocv_eval(ocv_coefs, z_ode(j, i));
          volts(j, i) = OCV_ode(j, i) + w_ode(j, i);
      }

      for (int i = 1; i < n_par; ++i) {
          phi_init(i - 1) = volts(j, i) - volts(j, 0);
      }
      phi << I, phi_init;

      i_branch.col(j) = A22.inverse() * phi; // Get the current going into each parallel branch
      current_sum(j) = i_branch.col(j).sum() / I;
  }
  */
  // the commented out part is for battery behaviour analysis might be unnecessary for SLIDE

  return 0;
}
