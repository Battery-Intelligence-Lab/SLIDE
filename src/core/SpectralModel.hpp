/**
 * @file SpectralModel.hpp
 * @brief Validated cold compiler for the spherical Chebyshev diffusion model.
 * @surface internal
 */

#pragma once

#include "CellDesign.hpp"
#include "Numeric.hpp"
#include "../types/Status.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numbers>

namespace slide::core {

template <int NCH>
struct CompiledSpectralModel
{
  static constexpr int output_nodes = NCH + 1;
  static constexpr int full_nodes = 2 * NCH + 3;

  PerDomain<std::array<real_t, NCH>> A{};
  PerDomain<std::array<real_t, NCH>> B{};
  PerDomain<std::array<std::array<real_t, NCH>, output_nodes>> C{};
  PerDomain<std::array<real_t, output_nodes>> D{};
  PerDomain<std::array<std::array<real_t, NCH>, NCH>> state_transform{};
  PerDomain<int> zero_mode{};
  std::array<real_t, NCH> x_inner{};
  std::array<real_t, output_nodes> Cc{};
  real_t cc_coeff{};
  std::array<std::array<real_t, full_nodes>, full_nodes> integration{};
};

namespace detail {

  template <int NCH>
  slide::Status validateCompiledSpectralModelFiniteness(
    const CompiledSpectralModel<NCH> &model)
  {
    // Pass by reference so fast-math cannot attach a finite-value assumption to
    // the predicate parameter before this explicit cold-path validation.
    const auto finite = [](const double &value) { return is_finite(value); };
    bool valid = true;
    for (const Domain domain : domains) {
      const auto d = domain_index(domain);
      valid = valid
              && std::all_of(model.A[d].begin(), model.A[d].end(), finite)
              && std::all_of(model.B[d].begin(), model.B[d].end(), finite)
              && std::all_of(model.D[d].begin(), model.D[d].end(), finite);
      for (const auto &row : model.C[d])
        valid = valid && std::all_of(row.begin(), row.end(), finite);
      for (const auto &row : model.state_transform[d])
        valid = valid && std::all_of(row.begin(), row.end(), finite);
    }
    valid = valid
            && std::all_of(model.x_inner.begin(), model.x_inner.end(), finite)
            && std::all_of(model.Cc.begin(), model.Cc.end(), finite)
            && is_finite(model.cc_coeff);
    for (const auto &row : model.integration)
      valid = valid && std::all_of(row.begin(), row.end(), finite);
    return valid ? slide::Status::Success : slide::Status::Numerical_failure;
  }

  template <int N>
  Eigen::Matrix<double, N + 1, N + 1> cumulativeIntegrationMatrix()
  {
    using std::numbers::pi;
    const Eigen::Vector<double, 2 * N> index = Eigen::Vector<double, 2 * N>::LinSpaced(2 * N, 0, 2 * N - 1);
    const Eigen::Matrix<double, N + 1, N + 1> values_to_coefficients_source = (((pi / N) * index.head(N + 1)
                                                                                * index.transpose().head(N + 1).rowwise().reverse()))
                                                                                .array()
                                                                                .cos()
                                                                                .matrix()
                                                                                .transpose();
    const Eigen::Matrix<double, N + 1, 2 * N> cosine = (((pi / N) * index.head(N + 1) * index.transpose()))
                                                         .array()
                                                         .cos()
                                                         .matrix();
    Eigen::Matrix<double, N + 1, N + 1> values_to_coefficients;
    values_to_coefficients.leftCols(1) = cosine.col(N) / N;
    values_to_coefficients.rightCols(1) = cosine.col(0) / N;
    values_to_coefficients.middleCols(1, N - 1) = (cosine.middleCols(1, N - 1).rowwise().reverse()
                                                   + cosine.middleCols(N + 1, N - 1))
                                                  / N;
    values_to_coefficients.row(0) /= 2.0;
    values_to_coefficients.row(N) /= 2.0;

    Eigen::Matrix<double, N + 1, N + 1> integrate = Eigen::Matrix<double, N + 1, N + 1>::Zero();
    for (int i = 1; i < N; ++i) {
      integrate(i, i + 1) = -1.0 / (2 * i);
      integrate(0, i + 1) += i % 2 == 0 ? -integrate(i, i + 1)
                                        : integrate(i, i + 1);
    }
    for (int i = 1; i <= N; ++i) {
      integrate(i, i - 1) = 1.0 / (2 * i);
      integrate(0, i - 1) += i % 2 == 0 ? -integrate(i, i - 1)
                                        : integrate(i, i - 1);
    }
    integrate.col(0) *= 2.0;
    Eigen::Matrix<double, N + 1, N + 1> result = values_to_coefficients_source * integrate * values_to_coefficients;
    result.row(0).setZero();
    return result;
  }

  inline double sphericalEigenRoot(int mode)
  {
    const double pi = std::numbers::pi;
    double left = static_cast<double>(mode) * pi;
    double right = (static_cast<double>(mode) + 0.5) * pi;
    left = std::nextafter(left, right);
    right = std::nextafter(right, left);
    for (int iteration = 0; iteration < 100; ++iteration) {
      const double middle = 0.5 * (left + right);
      if (std::tan(middle) - middle > 0.0)
        right = middle;
      else
        left = middle;
    }
    return 0.5 * (left + right);
  }

  template <int NCH>
  constexpr double fundamentalRootTolerance()
  {
    if constexpr (NCH == 5)
      return 5e-5;
    else if constexpr (NCH == 8)
      return 1e-9;
    else
      return 1e-12;
  }

  template <int NCH>
  constexpr int resolvedModes()
  {
    if constexpr (NCH == 5)
      return 1;
    else if constexpr (NCH == 8)
      return 3;
    else
      return 6;
  }

  template <int NCH, class Eigenvalues>
  slide::Status validateSpectrum(const Eigenvalues &eigenvalues,
                                 double radius,
                                 int &zero_mode)
  {
    std::array<double, NCH - 1> numerical_roots{};
    Eigen::Index zero{};
    eigenvalues.array().abs().minCoeff(&zero);
    zero_mode = static_cast<int>(zero);
    const double mass_eigenvalue = eigenvalues(zero_mode) * radius * radius;
    if (!is_finite(mass_eigenvalue) || std::abs(mass_eigenvalue) > 1e-8)
      return slide::Status::Numerical_failure;

    int cursor = 0;
    for (int i = 0; i < NCH; ++i) {
      if (i == zero_mode)
        continue;
      const double dimensionless = eigenvalues(i) * radius * radius;
      if (!is_finite(dimensionless) || !(dimensionless < 0.0))
        return slide::Status::Numerical_failure;
      numerical_roots[static_cast<std::size_t>(cursor++)] = std::sqrt(-dimensionless);
    }
    std::sort(numerical_roots.begin(), numerical_roots.end());
    for (std::size_t mode = 1; mode < numerical_roots.size(); ++mode)
      if (!(numerical_roots[mode] > numerical_roots[mode - 1]))
        return slide::Status::Numerical_failure;
    for (int mode = 0; mode < resolvedModes<NCH>(); ++mode) {
      const double exact = sphericalEigenRoot(mode + 1);
      const double relative = std::abs(numerical_roots[static_cast<std::size_t>(mode)] - exact)
                              / exact;
      const double tolerance = mode == 0 ? fundamentalRootTolerance<NCH>() : 1e-3;
      if (!is_finite(relative) || relative > tolerance)
        return slide::Status::Numerical_failure;
    }
    return slide::Status::Success;
  }

} // namespace detail

/**
 * Compile one of the explicitly registered Chebyshev orders. The result is assigned only after
 * every eigensolver, complex-part, invertibility, analytic-spectrum, and finiteness gate passes.
 */
template <int NCH>
[[nodiscard]] slide::Status compileSpectralModel(
  const PerDomain<real_t> &particle_radius,
  CompiledSpectralModel<NCH> &output)
{
  static_assert(NCH == 5 || NCH == 8 || NCH == 12,
                "NCH must be an explicitly validated registry order");
  output = {};
  for (const Domain domain : domains)
    if (!(is_finite(particle_radius[domain_index(domain)])
          && particle_radius[domain_index(domain)] > 0.0))
      return slide::Status::Invalid_parameters;

  constexpr int N = NCH + 1;
  constexpr int M = 2 * N;
  constexpr int Ncheb = M + 1;
  constexpr double dtheta = std::numbers::pi / (Ncheb - 1);
  CompiledSpectralModel<NCH> candidate;

  Eigen::Vector<double, Ncheb> xm;
  for (int i = 0; i < Ncheb; ++i)
    xm(i) = std::sin((Ncheb - 1 - 2 * i) * dtheta / 2);
  const Eigen::Vector<double, NCH> x_inner = xm.template segment<NCH>(1);
  for (int i = 0; i < NCH; ++i)
    candidate.x_inner[static_cast<std::size_t>(i)] = x_inner(i);

  Eigen::Matrix<double, Ncheb, Ncheb> differentiation = Eigen::Matrix<double, Ncheb, Ncheb>::Identity();
  for (int i = 0; i < Ncheb; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < Ncheb; ++j) {
      if (i == j)
        continue;
      const double dx = std::cos(dtheta * i) - std::cos(dtheta * j);
      double coefficient = 1.0;
      if (i == 0 || i == Ncheb - 1)
        coefficient *= 2.0;
      if (j == 0 || j == Ncheb - 1)
        coefficient /= 2.0;
      if ((i + j) % 2 == 1)
        coefficient = -coefficient;
      differentiation(i, j) = coefficient * differentiation(i, i) / dx;
      row_sum -= differentiation(i, j);
    }
    differentiation(i, i) = row_sum;
  }
  const Eigen::RowVector<double, Ncheb> first_derivative = differentiation.row(0);

  constexpr int derivative_order = 2;
  for (int i = 0; i < Ncheb; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < Ncheb; ++j) {
      if (i == j)
        continue;
      const double dx = std::cos(dtheta * i) - std::cos(dtheta * j);
      double coefficient = 1.0;
      if (i == 0 || i == Ncheb - 1)
        coefficient *= 2.0;
      if (j == 0 || j == Ncheb - 1)
        coefficient /= 2.0;
      if ((i + j) % 2 == 1)
        coefficient = -coefficient;
      differentiation(i, j) = derivative_order
                              * (coefficient * differentiation(i, i) - differentiation(i, j)) / dx;
      row_sum -= differentiation(i, j);
    }
    differentiation(i, i) = row_sum;
  }

  const Eigen::Matrix<double, N, Ncheb> second_derivative = differentiation.topRows(N);
  const Eigen::Matrix<double, N, N> folded_second = second_derivative.leftCols(N)
                                                    - second_derivative.rightCols(N).rowwise().reverse();
  const Eigen::RowVector<double, N> folded_first = first_derivative.leftCols(N)
                                                   - first_derivative.rightCols(N).rowwise().reverse();
  const double boundary_denominator = 1.0 - folded_first(0, 0);
  // This coefficient depends only on the fixed Chebyshev nodes for the three
  // registered orders, not on caller input. Registry-order compilation tests
  // exercise the invariant directly.
  assert(is_finite(boundary_denominator) && boundary_denominator != 0.0);

  const Eigen::Matrix<double, NCH, NCH> dimensionless_A = folded_second.template block<NCH, NCH>(1, 1)
                                                          + folded_second.template block<NCH, 1>(1, 0)
                                                              * folded_first.template block<1, NCH>(0, 1) / boundary_denominator;
  const Eigen::Matrix<double, NCH, 1> input = folded_second.template block<NCH, 1>(1, 0) / boundary_denominator;
  const Eigen::Matrix<double, 1, NCH> surface_output = folded_first.template block<1, NCH>(0, 1) / boundary_denominator;
  const double direct_output = 1.0 / boundary_denominator;

  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    const double radius = particle_radius[d];
    const Eigen::Matrix<double, NCH, NCH> physical_A = dimensionless_A / (radius * radius);
    Eigen::Matrix<double, N, NCH> physical_C;
    physical_C.row(0) = surface_output / radius;
    physical_C.bottomRows(NCH) = (x_inner * radius).array().inverse().matrix().asDiagonal();
    Eigen::Vector<double, N> physical_D = Eigen::Vector<double, N>::Zero();
    physical_D(0) = radius * direct_output;

    Eigen::EigenSolver<Eigen::Matrix<double, NCH, NCH>> solver(physical_A);
    if (solver.info() != Eigen::Success)
      return slide::Status::Numerical_failure;
    const double real_scale = solver.eigenvalues().real().array().abs().maxCoeff();
    const double eigenvalue_imag = solver.eigenvalues().imag().array().abs().maxCoeff();
    const double eigenvector_imag = solver.eigenvectors().imag().array().abs().maxCoeff();
    if (!(is_finite(real_scale) && real_scale > 0.0
          && is_finite(eigenvalue_imag) && is_finite(eigenvector_imag))
        || eigenvalue_imag > 1e-12 * real_scale
        || eigenvector_imag > 1e-10)
      return slide::Status::Numerical_failure;

    const Eigen::Vector<double, NCH> eigenvalues = solver.eigenvalues().real();
    int zero_mode{};
    const auto spectrum_status = detail::validateSpectrum<NCH>(eigenvalues, radius, zero_mode);
    if (spectrum_status != slide::Status::Success)
      return spectrum_status;
    const Eigen::Matrix<double, NCH, NCH> eigenvectors = solver.eigenvectors().real();
    Eigen::FullPivLU<Eigen::Matrix<double, NCH, NCH>> invertibility(eigenvectors);
    // validateSpectrum proves all accepted real eigenvalues are distinct;
    // nevertheless, retain the finite-precision decomposition check because a
    // backend may return a numerically rank-deficient approximate basis.
    if (!invertibility.isInvertible())
      return slide::Status::Numerical_failure;

    const Eigen::Vector<double, NCH> modal_input = eigenvectors.lu().solve(input);
    const Eigen::Matrix<double, N, NCH> modal_output = physical_C * eigenvectors;
    const Eigen::Matrix<double, NCH, NCH> state_transform = eigenvectors.inverse().eval();
    for (int mode = 0; mode < NCH; ++mode) {
      candidate.A[d][static_cast<std::size_t>(mode)] = mode == zero_mode ? 0.0 : eigenvalues(mode);
      candidate.B[d][static_cast<std::size_t>(mode)] = modal_input(mode);
      for (int column = 0; column < NCH; ++column)
        candidate.state_transform[d][static_cast<std::size_t>(mode)]
                                 [static_cast<std::size_t>(column)] = state_transform(mode, column);
    }
    candidate.zero_mode[d] = zero_mode;
    for (int row = 0; row < N; ++row) {
      candidate.D[d][static_cast<std::size_t>(row)] = physical_D(row);
      for (int mode = 0; mode < NCH; ++mode)
        candidate.C[d][static_cast<std::size_t>(row)]
                   [static_cast<std::size_t>(mode)] = modal_output(row, mode);
    }
  }

  const Eigen::Vector<double, N> centre_output = first_derivative.leftCols(N)
                                                 + first_derivative.rightCols(N).rowwise().reverse();
  candidate.cc_coeff = -1.0 / first_derivative(N);
  for (int i = 0; i < N; ++i)
    candidate.Cc[static_cast<std::size_t>(i)] = centre_output(i);
  const auto integration = detail::cumulativeIntegrationMatrix<M>();
  for (int row = 0; row < M + 1; ++row)
    for (int column = 0; column < M + 1; ++column)
      candidate.integration[static_cast<std::size_t>(row)]
                           [static_cast<std::size_t>(column)] = integration(row, column);

  const auto finite_status =
    detail::validateCompiledSpectralModelFiniteness(candidate);
  if (finite_status != slide::Status::Success)
    return finite_status;

  output = candidate;
  return slide::Status::Success;
}

} // namespace slide::core
