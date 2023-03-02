/*
 * parallelisation.hpp
 *
 * Some utility functions for parallelisation.

 *  Created on: 16 Oct 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include "../settings/settings.hpp"

#include <thread>
#include <vector>

namespace slide {

template <typename Tfun> // #TODO change with parallel algorithms.
void run(Tfun task_indv, int i_end, unsigned int numMaxParallelWorkers = settings::numMaxParallelWorkers)
{

  auto task_par = [&](int i_begin, int i_end_, int Nth) {
    while (i_begin < i_end_) {
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

      for (unsigned int i_begin = 0; i_begin < N_th_max; i_begin++) //!< indices for the threads
      {
        //!< Multi threaded simul:

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