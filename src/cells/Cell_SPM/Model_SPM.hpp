/**
 * @file Model_SPM.hpp
 * @brief  Defines a struct to store the matrices for the spatial discretisation of the solid diffusion PDE
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 2021
 */

#pragma once

#include "../../settings/settings.hpp"
#include "enum_definitions.hpp"

#include <Eigen/Dense>
#include <Eigen/LU>

#include <array>
#include <cassert>
#include <iostream>
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

  std::array<Eigen::Vector<double, nch>, 2> A, B; //!< only main diagonal is non-zero, so only store those values
  std::array<Eigen::Matrix<double, N, nch>, 2> C;
  std::array<Eigen::Vector<double, N>, 2> D;
  std::array<Eigen::Matrix<double, nch, nch>, 2> V; //!< inverse of the eigenvectors for the positive/negative electrode

  Eigen::Vector<double, N> Cc;           //!< matrix to get the concentration at the centre node
  double cc_coeff{};                     //!< coefficient for centre concentration: c_centre = cc_coeff * (Cc.c + flux*R/D)
  Eigen::Matrix<double, M + 1, M + 1> Q; //!< Matrix for Chebyshev integration

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

    const Eigen::Matrix<double, N - 1, N - 1> A_ = DN2.template block<N - 1, N - 1>(1, 1) + DN2.template block<N - 1, 1>(1, 0) * DN1.template block<1, N - 1>(0, 1) / temp;
    const Eigen::Matrix<double, N - 1, 1> B_ = DN2.template block<N - 1, 1>(1, 0) / temp;
    const Eigen::Matrix<double, 1, N - 1> C_ = DN1.template block<1, N - 1>(0, 1) / temp;
    const double D_ = 1.0 / temp;

    const Eigen::Matrix<double, N - 1, N - 1> A1 = A_ / (Rn * Rn);
    const Eigen::Matrix<double, N - 1, 1> B1 = B_;
    Eigen::Matrix<double, N, nch> C1;

    C1.row(0) = C_ / Rn;
    C1.bottomRows(nch) = xn.array().inverse().matrix().asDiagonal();

    Eigen::Vector<double, N> D1 = Eigen::Vector<double, N>::Zero();
    D1(0) = Rn * D_;

    const Eigen::Matrix<double, N - 1, N - 1> A3 = A_ / (Rp * Rp);
    const Eigen::Matrix<double, N - 1, 1> B3 = B_;
    Eigen::Matrix<double, N, nch> C3;

    C3.row(0) = C_ / Rp;
    C3.bottomRows(nch) = xp.array().inverse().matrix().asDiagonal();

    Eigen::Vector<double, N> D3 = Eigen::Vector<double, N>::Zero();
    D3(0) = Rp * D_;

    // REVIEW-MARK(2026-07-09, Fable): v4 build() hardening targets (full audit:
    // .claude/reports/chebyshev-math-audit-2026-07-09.md). (C1) es1/es3.info() unchecked — failed
    // eigensolve consumed silently; (C2) imag-part contamination only warns to stderr then truncates
    // via .real() — v4 build() must Status-fail (cold path); (C3) zero-mode index selection + the
    // Release-silent assert below — superseded by checking eigenvalues against the ANALYTIC roots of
    // tan(mu)=mu (audit §3; mass mode exact, modes matched by value not solver order). Legacy file
    // behaviour unchanged per PLAN §5.1 — these land in the v4 per-batch model build.
    Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es1(A1);
    {
      const double imag_max = es1.eigenvalues().imag().array().abs().maxCoeff();
      const double real_max = es1.eigenvalues().real().array().abs().maxCoeff();
      if (real_max > 0 && imag_max > 1e-10 * real_max)
        std::cerr << "WARNING: Model_SPM eigenvalues (neg) have significant imaginary parts "
                  << "(ratio=" << imag_max / real_max << ") for nch=" << nch << "\n";
    }
    V[neg] = es1.eigenvectors().real();
    A[neg] = es1.eigenvalues().real();
    B[neg] = V[neg].lu().solve(B1);
    C[neg] = C1 * V[neg];
    D[neg] = D1;

    Eigen::EigenSolver<Eigen::Matrix<double, N - 1, N - 1>> es3(A3);
    {
      const double imag_max = es3.eigenvalues().imag().array().abs().maxCoeff();
      const double real_max = es3.eigenvalues().real().array().abs().maxCoeff();
      if (real_max > 0 && imag_max > 1e-10 * real_max)
        std::cerr << "WARNING: Model_SPM eigenvalues (pos) have significant imaginary parts "
                  << "(ratio=" << imag_max / real_max << ") for nch=" << nch << "\n";
    }
    V[pos] = es3.eigenvectors().real();
    A[pos] = es3.eigenvalues().real();
    B[pos] = V[pos].lu().solve(B3);
    C[pos] = C3 * V[pos];
    D[pos] = D3;

    V[pos] = V[pos].inverse().eval(); //!< .eval() required to avoid Eigen aliasing for small matrices (nch <= 4)
    V[neg] = V[neg].inverse().eval();

    //!< Find the zero eigenvalue using relative threshold (matching MATLAB reference).
    //!< Normalise eigenvalues by max absolute value, then find the one closest to zero.
    //!< Do this independently for each electrode to avoid index mismatch from EigenSolver ordering.
    auto findZeroEigenvalue = [](const auto &eigenvalues) -> Eigen::Index {
      const double max_abs = eigenvalues.array().abs().maxCoeff();
      Eigen::Index idx{};
      (eigenvalues.array().abs() / max_abs).minCoeff(&idx);
      return idx;
    };

    const Eigen::Index zero_pos = findZeroEigenvalue(A[pos]);
    const Eigen::Index zero_neg = findZeroEigenvalue(A[neg]);

    A[pos](zero_pos) = 0.0;
    A[neg](zero_neg) = 0.0;

    assert(zero_pos == zero_neg && "Zero eigenvalue index mismatch between electrodes");
    zero = zero_pos;

    Cc = DM1.leftCols(N) + DM1.rightCols(N).rowwise().reverse();
    cc_coeff = -1.0 / DM1(N); //!< DM1(N) = 2*(-1)^N, so cc_coeff = -0.5 for even N, +0.5 for odd N
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

       << "A[neg]: \n"
       << model.A[neg] << "\n\n"

       << "B[neg]: \n"
       << model.B[neg] << "\n\n"

       << "C[neg]: \n"
       << model.C[neg] << "\n\n"

       << "D[neg]: \n"
       << model.D[neg] << "\n\n"

       << "A[pos]: \n"
       << model.A[pos] << "\n\n"

       << "B[pos]: \n"
       << model.B[pos] << "\n\n"

       << "C[pos]: \n"
       << model.C[pos] << "\n\n"

       << "D[pos]: \n"
       << model.D[pos] << "\n\n"

       << "V[neg]: \n"
       << model.V[neg] << "\n\n"

       << "V[pos]: \n"
       << model.V[pos] << "\n\n"

       << "Q: \n"
       << model.Q << "\n\n"

       << "Cc: \n"
       << model.Cc << "\n\n";
}


} // namespace slide