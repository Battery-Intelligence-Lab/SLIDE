/*
 * MCMC.hpp
 *
 * Markov Chain Monte Carlo
 *
 *  Created on: 12 Nov 2021
 *  Author(s): Volkan Kumtepeli
 *
 *
 * See : Aitio, Antti, Scott G. Marquis, Pedro Ascencio, and David Howey.
 * "Bayesian parameter estimation applied to the Li-ion battery single particle model with electrolyte dynamics."
 * IFAC-PapersOnLine 53, no. 2 (2020): 12497-12504.
 *
 * Copyright (c) 2021, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <cstdlib>
#include <array>

namespace slide {
template <size_t N>
auto S_update(slide::Matrix<double, N, N> &S, std::array<double, N> w, double alpha, double n)
{
  constexpr double gamma = 0.5; //!< a parameter determining the speed of adaptation of the proposal density covariance matrix.
  constexpr double alpha_star = 0.234;

  const double m = std::pow(n, -gamma) * (alpha - alpha_star) / slide::norm_sq<2>(w);
  const double m_root = std::sqrt(std::abs(m));

  for (auto &w_i : w)
    w_i *= m_root;

  auto I = slide::eye<N>(1.0);

  slide::cholUpdate<N>(I, w, m < 0);

  auto S_new = slide::zeros<N, N>();

  for (size_t i = 0; i < w.size(); i++)
    for (size_t j = 0; j <= i; j++)
      for (size_t o = 0; o <= i; o++)
        S_new[i][j] += S[i][o] * I[o][j];

  return S_new; //!< S_new;
}

template <size_t N_PARAM, size_t N_OUTPUT = 1>
struct MCMC_Input
{
  std::array<double, N_PARAM> param;
  std::array<double, N_OUTPUT> logR;

  auto scale(const std::array<double, N_PARAM> &param_scalar) const
  {
    auto param_scaled = param;

    for (size_t i = 0; i < param_scaled.size(); i++)
      param_scaled[i] *= param_scalar[i];

    return param_scaled;
  }

  auto &operator[](size_t i)
  {
    if (i < N_PARAM)
      return param[i];
    else if (i == N_PARAM)
      return logR[i - N_PARAM];
    else {
      std::cout << "THIS SHOULD NOT HAPPEN out of bounds for MCMC_Input.\n";
      throw 123511;
    }
  }

  auto operator[](size_t i) const
  {
    if (i < N_PARAM)
      return param[i];
    else if (i == N_PARAM)
      return logR[i - N_PARAM];
    else {
      std::cout << "THIS SHOULD NOT HAPPEN out of bounds for MCMC_Input.\n";
      throw 123511;
    }
  }

  auto input_size() const { return N_PARAM; }

  auto output_size() const { return N_OUTPUT; }

  auto size() const { return N_PARAM + N_OUTPUT; }

  auto get_R() const { return std::exp(logR[0]); };
};

template <size_t N_PARAM, size_t N_OUTPUT = 1>
struct MCMCSettings
{
  MCMC_Input<N_PARAM, N_OUTPUT> theta_init{};
  std::array<double, N_PARAM> theta_scalar{};
  Matrix<double, N_PARAM, 2> distParam{};

  size_t n_iter = 1000;
};

template <typename Tdistribution>
struct Prior
{
  std::vector<Tdistribution> distributionVec;

  Prior(std::vector<Tdistribution> distributionVec) : distributionVec(distributionVec) {}

  template <typename Tparam>
  auto operator()(const Tparam &b)
  {
    for (const auto b_i : b)
      if (b_i <= 0)
        return std::pair(-1e25, false); //!< -Inf

    double sum{ 0 };
    for (size_t i = 0; i < b.size(); i++)
      sum += std::log(distributionVec[i](b[i]));

    return std::pair(sum, true);
  }
};

using LL_Result = std::array<double, 5>;
//!< struct LL_Result
//!< {
//!<     double LLout{0}, LLout1{0}, LLout2{0}, P{}, err_sqr{};
//!< };

template <typename ThetaType, typename PriorType, typename CostType, typename SettingsType>
auto LL(const ThetaType &theta_cand, PriorType &prior, const CostType &costFun, const size_t Nobservations, const SettingsType &mcmcSettings)
{
  auto [P, flag] = prior(theta_cand.param);

  if (!flag)
    return std::pair(LL_Result{}, false); //!< std::pair(P, flag);
  else {
    auto err_sqr = costFun(theta_cand.scale(mcmcSettings.theta_scalar));

    const auto R_volt_exp = theta_cand.get_R();

    const double LLout1 = 0.5 * std::log(2.0 * PhyConst::pi * R_volt_exp) * Nobservations;
    const double LLout2 = 0.5 * err_sqr / R_volt_exp;
    double LLout = -P + LLout1 + LLout2;

    return std::pair(LL_Result{ LLout, LLout1, LLout2, P, err_sqr }, flag); //!<   std::pair(LLout, true);
  }
}

template <size_t N_PARAM, size_t N_OUTPUT = 1>
struct MCMC_Record
{
  MCMC_Input<N_PARAM, N_OUTPUT> input;
  LL_Result LLt_struct;

  auto operator[](size_t i) const
  {
    if (i < input.size())
      return input[i];
    else
      return LLt_struct[i - input.size()];
  }

  auto back() const
  {
    return LLt_struct.back();
  }

  auto size() const
  {
    return input.size() + LLt_struct.size();
  }
};

} // namespace slide
