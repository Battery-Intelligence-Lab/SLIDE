/*
 * main.cpp
 *
 * This file implements the main-function.
 * Here, the user has to choose what to simulate by commenting things in (by placing \\ at the start of the line) and out (by removing the \\)
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENSE for more information.
 */

//!< Include header files
#include "slide.hpp"
#include "../benchmark/benchmarks.hpp"
#include "../examples/examples.hpp"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <range/v3/all.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <numbers>
#include <ctime>
#include <thread>
#include <array>
#include <tuple>
#include <random>
#include <cmath>
#include <iomanip>
#include <complex>

/**
 * @brief Computes the matrices A, B, C, and D for the particle diffusion state-space model of the single particle model for each electrode.
 * @tparam Nch The number of Chebyshev nodes used to discretize the diffusion PDE in each particle.
 */
template <int Nch>
struct Model
{
  constexpr static auto N = Nch + 1;
  constexpr static auto M = 2 * N;
  int zero{};
  Eigen::Matrix<double, Nch, Nch> Vn, Vp;
  Eigen::Matrix<double, M + 1, M + 1> Q;
  Eigen::Matrix<double, N, Nch> Cn, Cp;
  Eigen::Vector<double, Nch> Ap, An, Bn, Bp, Dn, Dp;
  Eigen::Vector<double, N> Cc;
};


/**
 * @brief Computes the Chebyshev integration matrix.
 * @param N The number of Chebyshev points.
 * @return The matrix that maps function values at N Chebyshev points to values of the integral of the interpolating polynomial at those points.
 */
template <int N>
Eigen::Matrix<double, N + 1, N + 1> cumsummat()
{
  using std::numbers::pi;
  Eigen::Vector<double, 2 * N> arr = Eigen::Vector<double, 2 * N>::LinSpaced(2 * N, 0, 2 * N - 1);

  // Matrix mapping coeffs -> values.
  Eigen::Matrix<double, N + 1, N + 1> T = (((pi / N) * arr.head(N + 1) * arr.transpose().head(N + 1).rowwise().reverse())).array().cos().matrix().transpose();
  Eigen::Matrix<double, N + 1, 2 * N> F_real = (((pi / N) * arr.head(N + 1) * arr.transpose())).array().cos().matrix();

  std::ofstream abc{ "testFreal.txt", std::ios::out };

  abc << "Freal: \n"
      << F_real << '\n';


  Eigen::Matrix<double, N + 1, N + 1> Tinv; // Matrix mapping values -> coeffs.

  Tinv.leftCols(1) = F_real.col(N) / N;
  Tinv.rightCols(1) = F_real.col(0) / N;
  Tinv.middleCols(1, N - 1) = (F_real.middleCols(1, N - 1).rowwise().reverse() + F_real.middleCols(N + 1, N - 1)) / N;

  Tinv.row(0) /= 2.0;
  Tinv.row(N) /= 2.0;

  abc << "\nT: \n"
      << T << '\n';

  abc << "\nTinv: \n"
      << Tinv << '\n';

  // Matrix mapping coeffs -> integral coeffs. Note that the highest order term is truncated.
  Eigen::Matrix<double, N + 1, N + 1> B = Eigen::Matrix<double, N + 1, N + 1>::Zero();
  for (int i = 1; i < N; ++i) // - diag(1./k2,1)
  {
    B(i, i + 1) = -1.0 / (2 * i);
    if (i % 2 == 0) // neg.
      B(0, i + 1) -= B(i, i + 1);
    else
      B(0, i + 1) += B(i, i + 1);
  }

  for (int i = 1; i <= N; ++i) // diag(1./(2*k),-1)
  {
    B(i, i - 1) = 1.0 / (2 * i);

    if (i % 2 == 0) // positive.
      B(0, i - 1) -= B(i, i - 1);
    else
      B(0, i - 1) += B(i, i - 1);
  }

  B.col(0) *= 2;

  abc << "\nBBB: \n"
      << B << '\n';

  Eigen::Matrix<double, N + 1, N + 1> Q = T * B * Tinv;

  Q.row(0).setZero();

  abc << "\nQ: \n"
      << Q << '\n';

  return Q;
}


template <int Nch>
Model<Nch> get_model_vk_slide_impl()
{
  using std::numbers::pi;
  Model<Nch> model;

  constexpr double Rn = 12.5e-6; // Anode solid particles' radius [m]
  constexpr double Rp = 8.5e-6;  // Cathode solid particles' radius [m]

  constexpr int N = Nch + 1;
  constexpr int M = 2 * N;
  constexpr int Ncheb = M + 1;
  constexpr double dtheta = pi / (Ncheb - 1);

  // Computational coordinates (Chebyshev nodes)
  Eigen::Vector<double, Ncheb> xm;
  for (int i = 0; i < Ncheb; ++i)
    xm(i) = std::sin((Ncheb - 1 - 2 * i) * dtheta / 2);

  Eigen::Vector<double, N - 1> xr = xm.segment<N - 1>(1);
  Eigen::Vector<double, N - 1> xp = xr * Rp;
  Eigen::Vector<double, N - 1> xn = xr * Rn;

  // Computing the Chebyshev differentiation matrices
  Eigen::Matrix<double, Ncheb, Ncheb> D_vk = Eigen::Matrix<double, Ncheb, Ncheb>::Identity();
  for (int i = 0; i < Ncheb; ++i) {
    double row_sum = 0;
    for (int j = 0; j < Ncheb; ++j) {
      if (i == j) continue;

      const double DX = std::cos(dtheta * i) - std::cos(dtheta * j);
      double C_vk = 1;

      if (i == 0 || i == Ncheb - 1) C_vk *= 2;
      if (j == 0 || j == Ncheb - 1) C_vk /= 2;
      if ((i + j) % 2 == 1) C_vk = -C_vk;

      D_vk(i, j) = C_vk * D_vk(i, i) / DX;
      row_sum -= D_vk(i, j);
    }
    D_vk(i, i) = row_sum;
  }

  Eigen::RowVector<double, Ncheb> DM1 = D_vk.row(0);

  constexpr int order = 2;
  for (int i = 0; i < Ncheb; ++i) {
    double row_sum = 0;
    for (int j = 0; j < Ncheb; ++j) {
      if (i == j) continue;

      const double DX = std::cos(dtheta * i) - std::cos(dtheta * j);
      double C_vk = 1;

      if (i == 0 || i == Ncheb - 1) C_vk *= 2;
      if (j == 0 || j == Ncheb - 1) C_vk /= 2;
      if ((i + j) % 2 == 1) C_vk = -C_vk;

      D_vk(i, j) = order * (C_vk * D_vk(i, i) - D_vk(i, j)) / DX;
      row_sum -= D_vk(i, j);
    }
    D_vk(i, i) = row_sum;
  }

  Eigen::Matrix<double, N, Ncheb> DM2 = D_vk.topRows(N);
  Eigen::Matrix<double, N, N> DN2 = DM2.leftCols(N) - DM2.rightCols(N).rowwise().reverse();
  Eigen::RowVector<double, N> DN1 = DM1.leftCols(N) - DM1.rightCols(N).rowwise().reverse();

  const double temp = (1 - DN1(0, 0));

  Eigen::Matrix<double, N - 1, N - 1> A = DN2.block<N - 1, N - 1>(1, 1) + DN2.block<N - 1, 1>(1, 0) * DN1.block<1, N - 1>(0, 1) / temp;
  Eigen::Matrix<double, N - 1, 1> B = DN2.block<N - 1, 1>(1, 0) / temp;
  Eigen::Matrix<double, 1, N - 1> C = DN1.block<1, N - 1>(0, 1) / temp;
  double D = 1.0 / temp;

  Eigen::Matrix<double, N - 1, N - 1> A1 = A / (Rn * Rn);
  Eigen::Matrix<double, N - 1, 1> B1 = B;
  Eigen::Matrix<double, N, Nch> C1;

  C1.row(0) = C / Rn;
  C1.bottomRows(Nch) = xn.array().inverse().matrix().asDiagonal();

  Eigen::Vector<double, Nch> D1 = Eigen::Vector<double, Nch>::Zero();
  D1(0) = Rn * D;

  Eigen::Matrix<double, N - 1, N - 1> A3 = A / (Rp * Rp);
  Eigen::Matrix<double, N - 1, 1> B3 = B;
  Eigen::Matrix<double, N, Nch> C3;

  C3.row(0) = C / Rp;
  C3.bottomRows(Nch) = xp.array().inverse().matrix().asDiagonal();

  Eigen::Vector<double, Nch> D3 = Eigen::Vector<double, Nch>::Zero();
  D3(0) = Rp * D;

  Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es1(A1);
  model.Vn = es1.eigenvectors().real();
  model.An = es1.eigenvalues().real();
  model.Bn = model.Vn.lu().solve(B1);
  model.Cn = C1 * model.Vn;
  model.Dn = D1;

  Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es3(A3);
  model.Vp = es3.eigenvectors().real();
  model.Ap = es3.eigenvalues().real();
  model.Bp = model.Vp.lu().solve(B3);
  model.Cp = C3 * model.Vp;
  model.Dp = D3;

  model.Vp = model.Vp.inverse();
  model.Vn = model.Vn.inverse();

  Eigen::Index minIndex;
  model.Ap.array().abs().minCoeff(&minIndex);

  model.Ap(minIndex) = 0.0;
  model.An(minIndex) = 0.0;

  model.zero = minIndex;

  model.Cc = DM1.leftCols(N) + DM1.rightCols(N).rowwise().reverse();
  model.Q = cumsummat<M>();

  return model;
}


int main()
{

  constexpr int Nch = 5; // Example value
  Model<Nch> model = get_model_vk_slide_impl<Nch>();

  // std::cout << "Anode A matrix:\n"
  //           << model.A1 << "\n";
  // std::cout << "Anode B matrix:\n"
  //           << model.B1 << "\n";
  // std::cout << "Anode C matrix:\n"
  //           << model.C1 << "\n";
  // std::cout << "Anode D matrix:\n"
  //           << model.D1 << "\n";

  // std::cout << "Cathode A matrix:\n"
  //           << model.A3 << "\n";
  // std::cout << "Cathode B matrix:\n"
  //           << model.B3 << "\n";
  // std::cout << "Cathode C matrix:\n"
  //           << model.C3 << "\n";
  // std::cout << "Cathode D matrix:\n"
  //           << model.D3 << "\n";


  /*
   * This function decides what will be simulated.
   * In the code below, you have to uncomment the line which you want to execute (to uncomment, remove the //!<in front of the line)
   * and comment all the other lines in the blocks with FUNCTION CALLS (to comment in, add //!<in front of the lines)
   */

  // estimateCharacterisation();


  //!< print that you start simulations
  //!< slide::tests::unit::test_all();
  //!< Make a clock to measure how long the simulation takes
  slide::Clock clk;
  //!< Read the values for the matrices for the spatial discretisation of the solid diffusion PDE.
  //!< these values are calculated by the MATLAB-script 'modelSetup.m', which writes them to csv files.
  //!< Here, we only need to read these csv files

  //!< slide::Model_SPM M; //!<Initialised at constructor.
  //!<  Look for extra settings at Constants.hpp (such as verbosity)

  //!< Choose the prefix, a string with which the names of all subfolders with the outputfiles will start (the folders will be called prefix_xxxx).
  //!< This is typically a number, and the reason why we use it is that results of previous simulations are not overwritten.
  //!< i.e. for one simulation, use pref "0", then the next simulation could have pref "1".
  //!< this ensures the results of the second simulation are written in different files and folders, rather than overwriting the files from the first simulation
  const std::string pref = "0"; //!< todo
                                //!<  will be used to name a folder with the results, so must comply with the naming convention on your operating system
                                //!<  letters, numbers, underscores and dots are fine
                                //!<  don't use spaces
                                //!<  don't use forbidden characters, e.g. for Windows the following are not allowed < > : " \ | / ? *
                                //!<  https://docs.microsoft.com/en-us/windows/desktop/fileio/naming-a-file

  //!< Choose which cell should be used for the simulations
  //!< const int cellType = slide::cellType::KokamNMC;
  //!< which cell to use for the simulation.
  //!< 0 	high power Kokam NMC cell (18650)
  //!< 1 	high energy LG Chem NMC cell (18650)
  //!< 2 	user cell (template class provided for where the user can define his own parameters)

  //!< Choose which degradation models to use for the simulations
  //!< The models are identified by a number, see DEG_ID in DEG_ID.hpp and below
  /* SEI_id 		identification for which model(s) to use for SEI growth
   * 					0	no SEI growth
   * 					1 	kinetic model only (Tafel kinetics)
   * 							ref: Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
   * 					2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
   * 							ref: Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
   * 					3	Christensen and Newman model: diffusion and kinetics model
   * 							ref: Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
   * SEI_porosity	integer deciding whether we reduce the volume fraction of active electrode material due to SEI growth
   *  				0	don't subtract anything
   * 					1	use correlation from Ashwin
   * 							ref: Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
   * CS_id	identification for which model(s) to use for surface crack growth
   *  				0 	no surface cracking
   * 					1 	Laresgoiti's stress + crack growth model
   * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 					2 	Dai stress model + Laresgoiti crack growth
   * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 								 Dai, Cai, White, Journal of Power sources 247, 2014
   * 					3 	model based on Deshpande and Bernardi
   * 							ref: Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
   * 					4 	model from Barai
   * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
   * 					5 	model from Ekstrom
   * 							ref: Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
   * CS_diffusion	integer deciding whether to reduce the negative diffusion constant due to surface cracks
   *  				0 	don't decrease diffusion
   * 					1	decrease according to Barai
   * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
   * LAM_id  identification for which model(s) to use for loss of active material
   *  				0 	no LAM
   * 					1	Dai's stress model and Laresgoiti's correlation to relate the stress to LAM
   * 							ref: Dai, Cai, White, Journal of Power sources 247, 2014
   * 								 Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 					2	delacourt's correlation between abs(current density) and porosity
   * 							ref: Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
   * 					3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
   * 							ref: Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
   * 					4 	Narayanrao's correlation between effective surface area and current density
   * 							ref: Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
   * pl_id 	identification for which model to use for li-plating
   *  				0 	no plating
   * 					1	Yang et al's model for thermodynamic plating (Tafel kinetics)
   * 							ref: Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
   */
  //!< For now, assume we want to include 'Pinson&Bazant'-SEI growth, 'Delacourt'-LAM, 'Kindermann'-LAM and 'Yang'-lithium plating
  slide::DEG_ID deg;
  deg.SEI_id.add_model(2); //!< Pinson & Bazant SEI growth
  deg.SEI_porosity = 0;    //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  deg.CS_id.add_model(0); //!< no surface cracks
  deg.CS_diffusion = 0;   //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  //!< there are 2 LAM models (Delacourt and Kindermann)
  deg.LAM_id.add_model(2); //!< Delacourt LAM
  deg.LAM_id.add_model(3); //!< Kindermann LAM

  // deg.pl_id = 1; //!< Yang lithium plating

  // std::cout << "Used ageing model: " << deg.print() << '\n';

  //!< Then the user has to choose what is simulated.
  //!< In the code below, uncomment the line which calls the function you want to execute (uncommenting means removing the //!<in front of the line)
  //!< and comment all the other lines (commenting means putting //!<in front of the line)

  auto numThreads = std::thread::hardware_concurrency(); //!<
  std::cout << "Available number of threads : " << numThreads << '\n';


  // Please see examples for using SLIDE. For previous version refer to SLIDE_v2 branch.
  // using namespace slide;

  // auto c = make<Cell_SPM>("MyCell", deg, 1, 1, 1, 1);
  // auto &st = c->getStateObj();

  // std::cout << "Voltage: " << c->V() << " SOC: " << 100 * st.SOC() << " %.\n";
  // c->setCurrent(16);
  // std::cout << "Voltage: " << c->V() << " SOC: " << 100 * st.SOC() << " %.\n";
  // c->timeStep_CC(1, 60);
  // std::cout << "Voltage: " << c->V() << " SOC: " << 100 * st.SOC() << " %.\n";

  // Cycler cyc(c, "Cycler1");

  // ThroughputData th{};
  // cyc.CC(-16, 4.2, 3600, 1, 1, th);
  // std::cout << "\nAfter CC charge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";

  // auto status = cyc.CV(4.2, 10e-3, 7200, 0.1, 1, th);
  // std::cout << "\nAfter CV charge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";
  // std::cout << getStatusMessage(status) << "\n";

  // cyc.CCCV(16, 2.7, 5e-3, 1, 1, th);
  // std::cout << "\nAfter CCCV discharge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";

  // cyc.writeData();
  // // Cycler:

  // Cycler cyc(c, "Cycler2");
  // ThroughputData th{};
  // cyc.CC(-16, 4.2, 3600, 1, 1, th);

  // std::cout << "\nAfter CC charge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";
  // std::cout << "Ah: " << th.Ah() << " Wh: " << th.Wh() << " time: " << th.time() << '\n';

  // auto status = cyc.CV(4.2, 10e-3, 7200, 1, 1, th);

  // std::cout << "Why did CV stopped is: " << getStatusMessage(status) << '\n';
  // std::cout << "\nAfter CV charge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";
  // std::cout << "Ah: " << th.Ah() << " Wh: " << th.Wh() << " time: " << th.time() << '\n';


  // auto status2 = cyc.CCCV(16, 2.7, 5e-3, 0.5, 1, th);
  // std::cout << "Why did CCCV stopped is: " << getStatusMessage(status2) << '\n';
  // std::cout << "\nAfter CCCV discharge:\n";
  // std::cout << "Voltage: " << c->V() << " I: " << 1000 * c->I() << "mA SOC: " << 100 * st.SOC() << " %.\n";
  // std::cout << "Ah: " << th.Ah() << " Wh: " << th.Wh() << " time: " << th.time() << '\n';

  // cyc.writeData();


  //!< Examples:
  //  slide::examples::drive_cycle_artemis();
  // slide::examples::GITT_test();
  //!< Benchmarks:

  // slide::benchmarks::run_Cell_Bucket();
  // slide::benchmarks::run_Cell_ECM();
  // slide::benchmarks::run_Cell_SPM_1(1);
  // slide::benchmarks::run_Cell_SPM_2(1);
  // slide::benchmarks::run_LP_case_SmallPack();
  // slide::benchmarks::run_LP_case_MediumPack();
  // slide::benchmarks::run_LP_case_LargePack();

  // MATLAB ECM benchmarks:
  // slide::benchmarks::run_Cell_Bucket_single_default_pulse();
  // slide::benchmarks::run_Cell_Bucket_single_default_CCCV();

  // slide::benchmarks::run_Cell_ECM_single_default_pulse();
  // slide::benchmarks::run_Cell_ECM_single_default_CCCV();

  // slide::benchmarks::run_Cell_ECM_2_RC_single_default_pulse();
  // slide::benchmarks::run_Cell_ECM_2_RC_single_default_CCCV();

  // slide::benchmarks::run_Cell_ECM_parallel_3_default_pulse();

  // slide::benchmarks::run_Cell_ECM_parallel_3_default_CCCV();

  // slide::benchmarks::run_Cell_ECM_parallel_3_withRcontact_CCCV();

  // slide::benchmarks::run_Cell_ECM_series_3_withRcontact_CCCV();

  //!<*********************************************** END ********************************************************
  //!< Now all the simulations have finished. Print this message, as well as how long it took to do the simulations
  std::cout << "finished all simulations in " << clk << ".\n";
}
