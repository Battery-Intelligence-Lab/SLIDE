/**
 * @file parallelisation.hpp
 * @brief Some utility functions for parallelisation.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 16 Oct 2022
 */

#pragma once

#include "../settings/settings.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <thread>
#include <vector>

namespace slide {

namespace legacy_parallel_detail {

  inline unsigned workerCount(unsigned requested,
                              unsigned hardware,
                              unsigned tasks) noexcept
  {
    if (tasks == 0)
      return 0;
    hardware = std::max(1U, hardware);
    if (requested == 0)
      requested = hardware;
    return std::max(1U, std::min({ requested, hardware, tasks }));
  }

} // namespace legacy_parallel_detail

template <typename Tfun> // #TODO change with parallel algorithms.
void run(Tfun task_indv, int i_end, unsigned int numMaxParallelWorkers = settings::numMaxParallelWorkers)
{
  if (i_end <= 0)
    return;

  if constexpr (settings::isParallel) {
    const unsigned workers = legacy_parallel_detail::workerCount(
      numMaxParallelWorkers,
      std::thread::hardware_concurrency(),
      static_cast<unsigned>(i_end));
    if (workers == 1) {
      std::exception_ptr first;
      for (int index = 0; index < i_end; ++index)
        try {
          task_indv(index);
        } catch (...) {
          if (first == nullptr)
            first = std::current_exception();
        }
      if (first != nullptr)
        std::rethrow_exception(first);
      return;
    }

    struct Failure
    {
      int index{ std::numeric_limits<int>::max() };
      std::exception_ptr exception{};
    };
    std::vector<Failure> failures(workers);
    std::vector<std::thread> threads;
    threads.reserve(workers);
    const auto task = [&](unsigned worker) {
      const auto end = static_cast<unsigned>(i_end);
      for (unsigned index = worker; index < end; index += workers) {
        try {
          task_indv(static_cast<int>(index));
        } catch (...) {
          if (failures[worker].exception == nullptr) {
            failures[worker].index = static_cast<int>(index);
            failures[worker].exception = std::current_exception();
          }
        }
      }
    };
    try {
      for (unsigned worker = 0; worker < workers; ++worker)
        threads.emplace_back(task, worker);
    } catch (...) {
      for (auto &th : threads) {
        if (th.joinable())
          th.join();
      }
      throw;
    }
    for (auto &thread : threads)
      if (thread.joinable())
        thread.join();

    const Failure *first{};
    for (const auto &failure : failures)
      if (failure.exception != nullptr
          && (first == nullptr || failure.index < first->index))
        first = &failure;
    if (first != nullptr)
      std::rethrow_exception(first->exception);
  } else {
    for (int index = 0; index < i_end; ++index)
      task_indv(index);
  }
}
} // namespace slide
