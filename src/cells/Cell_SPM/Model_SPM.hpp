/**
 * @file Model_SPM.hpp
 * @brief  Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 2021
 */

#pragma once

#include "../../settings/settings.hpp"

#include <Eigen/Dense>
#include <Eigen/LU>

#include <array>
#include <numbers>
#include <fstream>

namespace slide {
//!< Define a structure with the matrices of the spatial discretisation of the solid diffusion PDE
//!< See the matlab script modelSetup.m
//!< This class has 326 elements. So it should be allocated in heap once.

/**
 * @brief Computes the Chebyshev integration matrix.
 * @param N The number of Chebyshev points.
 * @return The matrix that maps function values at N Chebyshev points to values of the integral of the interpolating polynomial at those points.
 */
template <int N>
Eigen::Matrix<double, N + 1, N + 1> cumsummat();

template <int nch = settings::nch>
struct Model_SPM
{
  constexpr static auto N = nch + 1;
  constexpr static auto M = 2 * N;
  constexpr static auto Ncheb = M + 1;
  constexpr static double dtheta = std::numbers::pi / (Ncheb - 1);

  double Rn{ 12.5e-6 }; // Anode solid particles' radius [m]
  double Rp{ 8.5e-6 };  // Cathode solid particles' radius [m]

  Eigen::Index zero{};            // Location of the zero.
  Eigen::Vector<double, nch> xch; //!< location of the Chebyshev nodes in the positive domain EXCLUDING centre and surface

  //!< state space model
  //	dzpos/dt = Ap*zpos + Bp*jp		time derivative of (transformed) concentration at the inner nodes
  //	dzneg/dt = An*zneg + Bn*jn
  //	cp = Cp*zpos + Dp*jp			actual concentration [mol m-3] of the nodes (surface, inner)
  //!< 	cn = Cn*zpos + Dn*jn

  Eigen::Vector<double, nch> Ap, An; //!< only main diagonal is non-zero, so only store those values
  Eigen::Vector<double, nch> Bp, Bn;
  Eigen::Matrix<double, N, nch> Cp, Cn;
  Eigen::Vector<double, N> Cc, Dp, Dn;    //!< matrix to get the concentration at the centre node
  Eigen::Matrix<double, nch, nch> Vp, Vn; //!< inverse of the eigenvectors for the positive/negative electrode
  Eigen::Matrix<double, M + 1, M + 1> Q;  //!< Matrix for Chebyshev integration

  Model_SPM()
  {
    // Computational coordinates (Chebyshev nodes)
    Eigen::Vector<double, Ncheb> xm;
    for (int i = 0; i < Ncheb; ++i)
      xm(i) = std::sin((Ncheb - 1 - 2 * i) * dtheta / 2);

    xch = xm.template segment<N - 1>(1);
    const Eigen::Vector<double, N - 1> xp = xch * Rp;
    const Eigen::Vector<double, N - 1> xn = xch * Rn;

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

    const Eigen::Matrix<double, N, Ncheb> DM2 = D_vk.topRows(N);
    const Eigen::Matrix<double, N, N> DN2 = DM2.leftCols(N) - DM2.rightCols(N).rowwise().reverse();
    const Eigen::RowVector<double, N> DN1 = DM1.leftCols(N) - DM1.rightCols(N).rowwise().reverse();

    const double temp = (1 - DN1(0, 0));

    const Eigen::Matrix<double, N - 1, N - 1> A = DN2.template block<N - 1, N - 1>(1, 1) + DN2.template block<N - 1, 1>(1, 0) * DN1.template block<1, N - 1>(0, 1) / temp;
    const Eigen::Matrix<double, N - 1, 1> B = DN2.template block<N - 1, 1>(1, 0) / temp;
    const Eigen::Matrix<double, 1, N - 1> C = DN1.template block<1, N - 1>(0, 1) / temp;
    const double D = 1.0 / temp;

    const Eigen::Matrix<double, N - 1, N - 1> A1 = A / (Rn * Rn);
    const Eigen::Matrix<double, N - 1, 1> B1 = B;
    Eigen::Matrix<double, N, nch> C1;

    C1.row(0) = C / Rn;
    C1.bottomRows(nch) = xn.array().inverse().matrix().asDiagonal();

    Eigen::Vector<double, N> D1 = Eigen::Vector<double, N>::Zero();
    D1(0) = Rn * D;

    const Eigen::Matrix<double, N - 1, N - 1> A3 = A / (Rp * Rp);
    const Eigen::Matrix<double, N - 1, 1> B3 = B;
    Eigen::Matrix<double, N, nch> C3;

    C3.row(0) = C / Rp;
    C3.bottomRows(nch) = xp.array().inverse().matrix().asDiagonal();

    Eigen::Vector<double, N> D3 = Eigen::Vector<double, N>::Zero();
    D3(0) = Rp * D;

    Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es1(A1);
    Vn = es1.eigenvectors().real();
    An = es1.eigenvalues().real();
    Bn = Vn.lu().solve(B1);
    Cn = C1 * Vn;
    Dn = D1;

    Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es3(A3);
    Vp = es3.eigenvectors().real();
    Ap = es3.eigenvalues().real();
    Bp = Vp.lu().solve(B3);
    Cp = C3 * Vp;
    Dp = D3;

    Vp = Vp.inverse();
    Vn = Vn.inverse();

    Ap.array().abs().minCoeff(&zero);

    Ap(zero) = An(zero) = 0.0;

    Cc = DM1.leftCols(N) + DM1.rightCols(N).rowwise().reverse();
    Q = cumsummat<M>();
  }

  static Model_SPM *makeModel() //!< #TODO make other type of models possible.
  {
    static Model_SPM<> model;
    return &model;
  }
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
  Eigen::Matrix<double, N + 1, N + 1> Tinv; // Matrix mapping values -> coeffs.

  Tinv.leftCols(1) = F_real.col(N) / N;
  Tinv.rightCols(1) = F_real.col(0) / N;
  Tinv.middleCols(1, N - 1) = (F_real.middleCols(1, N - 1).rowwise().reverse() + F_real.middleCols(N + 1, N - 1)) / N;

  Tinv.row(0) /= 2.0;
  Tinv.row(N) /= 2.0;

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

  Eigen::Matrix<double, N + 1, N + 1> Q = T * B * Tinv;

  Q.row(0).setZero();

  return Q;
}

template <int nch>
void testwriteChebyshevModel()
{
  std::ofstream file{ "testChebyshev.txt", std::ios::out };

  auto &model = *Model_SPM<nch>::makeModel();


  file << "xch: \n"
       << model.xch << "\n\n"

       << "An: \n"
       << model.An << "\n\n"

       << "Bn: \n"
       << model.Bn << "\n\n"

       << "Cn: \n"
       << model.Cn << "\n\n"

       << "Dn: \n"
       << model.Dn << "\n\n"

       << "Ap: \n"
       << model.Ap << "\n\n"

       << "Bp: \n"
       << model.Bp << "\n\n"

       << "Cp: \n"
       << model.Cp << "\n\n"

       << "Dp: \n"
       << model.Dp << "\n\n"

       << "Vn: \n"
       << model.Vn << "\n\n"

       << "Vp: \n"
       << model.Vp << "\n\n"

       << "Q: \n"
       << model.Q << "\n\n"

       << "Cc: \n"
       << model.Cc << "\n\n";
}


} // namespace slide