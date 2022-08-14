/*
 * util.hpp
 *
 * Some utility functions. So the code is less verbose.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

// Include other util files.
#pragma once

#include <string>
#include <iostream>
#include <array>
#include <thread>
#include <vector>
#include <cmath>
#include <ctime>

//#include "util_error.hpp"
#include "../settings/settings.hpp"
#include "slide_aux.hpp"
#include "util_debug.hpp"

namespace slide {
struct Clock
{
  std::clock_t tstart{ std::clock() };
  Clock() = default;
  auto now() const { return std::clock(); }
  auto start() const { return tstart; }
  double duration() const { return (now() - start()) / (double)CLOCKS_PER_SEC; }
};

inline std::ostream &operator<<(std::ostream &ofs, const Clock &clk)
{
  const auto duration = clk.duration();
  ofs << std::floor(duration / 60) << ":"
      << duration - std::floor(duration / 60) * 60
      << " min:sec";

  return ofs;
}

constexpr auto abs_sqrt(auto x)
{
  return std::sqrt(std::abs(x));
}

template <typename T, typename FilePath>
void output_printer(std::vector<T> vec, const FilePath &save_path)
{
  std::ofstream out_file(save_path, std::ios_base::out);

  for (const auto &state : vec) {
    for (size_t i{ 0 }; i < (state.size() - 1); i++)
      out_file << state[i] << ',';

    out_file << state.back() << '\n';
  }

  out_file.close();
}

template <typename Tfun>
void run(Tfun task_indv, int i_end, unsigned int numMaxParallelWorkers = settings::numMaxParallelWorkers)
{

  auto task_par = [&](int i_begin, int i_end, int Nth) {
    while (i_begin < i_end) {
      task_indv(i_begin);
      i_begin += Nth;
    }
  };

  if constexpr (settings::isParallel) {
    if (numMaxParallelWorkers == 1)
      task_par(0, i_end, 1);
    else {
      if (numMaxParallelWorkers < 1)
        numMaxParallelWorkers = std::thread::hardware_concurrency();

      const unsigned int N_th_max = std::min(numMaxParallelWorkers, std::thread::hardware_concurrency());

      std::vector<std::thread> threads;
      threads.reserve(N_th_max);

      for (unsigned int i_begin = 0; i_begin < N_th_max; i_begin++) // indices for the threads
      {
        // Multi threaded simul:

        threads.emplace_back(task_par, i_begin, i_end, N_th_max);
      }

      for (auto &th : threads) {
        if (th.joinable())
          th.join();
      }
    }
  } else {
    task_par(0, i_end, 1);
  }
}
} // namespace slide

namespace slide::util {
// template <int pLevel>
// inline void print(const std::string& message) // Print messages.
//  {
//  if constexpr (settings::verbose >= pLevel)
//   {
//        std::cout << message << '\n';
//     }
//  }

// struct Counter // Counts how many times it is called.
// {

//     // To use counter create a static Counter object in the function you would like to measure
//     // slide::util::Counter myCounterName(printEveryXiter);  Then just call it
//     // myCounterName();
//     static count{0};
//     long long printEveryXiter{500000};

//     Counter() = default;
//     Counter(unsigned long long printEveryXiter) : printEveryXiter(printEveryXiter){};

//     void operator()()
//     {
//         count++;

//         if (count % printEveryXiter == 0)
//             std::cout << "Counter: " << count << '\n';
//     }
// };

} // namespace slide::util

namespace slide {
std::vector<double> linstep(double x_min, double x_step, int Nstep);
std::vector<double> logstep(double x_min, double x_step, int Nstep);

FixedData<double> linstep_fix(double x_min, double x_step, int Nstep);
FixedData<double> logstep_fix(double x_min, double x_step, int Nstep);

FixedData<double> range_fix(double x_min, double x_max, double dx);

// FixedData<int> range(int stop); #TODO -> FixedData is not good since it has function
// FixedData<int> range(int start, int stop, int step = 1);

std::vector<double> linspace(double x1, double x2, int N);
FixedData<double> linspace_fix(double x1, double x2, int N);

template <size_t N>
constexpr std::array<double, N> linspace(double x1, double x2)
{
  std::array<double, N> out;

  double dx{ 1 };
  if (N < 1)
    return out;
  else if (N > 1)
    dx = (x2 - x1) / static_cast<double>(N - 1);

  for (size_t n{ 0 }; n < N; n++) {
    out[N - 1 - n] = x2 - n * dx;
  }

  return out;
}

constexpr inline double C_to_Kelvin(double Celsius) { return PhyConst::Kelvin + Celsius; }
constexpr inline double K_to_Celsius(double Kelvin) { return Kelvin - PhyConst::Kelvin; }

} // namespace slide

namespace slide {

inline slide::FixedData<double> range_fix(double x_min, double x_max, double x_step)
{
  int Nstep = static_cast<int>((x_max - x_min) / x_step) + 1;
  return slide::FixedData<double>(x_min, x_step, Nstep);
}

inline slide::FixedData<double> linstep_fix(double x_min, double x_step, int Nstep)
{
  return slide::FixedData<double>(x_min, x_step, Nstep);
}

inline slide::FixedData<double> logstep_fix(double x_min, double x_step, int Nstep)
{
  auto fun = [](double x_min, double x_step, int i) { return x_min * std::pow(x_step, i); };

  return slide::FixedData<double>(x_min, x_step, Nstep, fun);
}

inline slide::FixedData<double> linspace_fix(double x1, double x2, int N)
{
  if (N < 1)
    return FixedData(x2, 0.0, 0);
  else if (N == 1)
    return FixedData(x2, 0.0, 1);
  else
    return FixedData(x1, (x2 - x1) / static_cast<double>(N - 1), N);
}

inline std::vector<double> linspace(double x1, double x2, int N)
{
  N = std::max(0, N);
  std::vector<double> out(N);

  double dx{ 1 };
  if (N < 1)
    return out;
  else if (N > 1)
    dx = (x2 - x1) / static_cast<double>(N - 1);

  for (int n{ 0 }; n < N; n++) {
    out[N - 1 - n] = x2 - n * dx;
  }

  return out;
}
} // namespace slide