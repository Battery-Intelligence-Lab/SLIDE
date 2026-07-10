/**
 * @file ThreadPool.hpp
 * @brief Persistent dependency-free batch task runtime and deterministic reductions.
 */

#pragma once

#include "Numeric.hpp"
#include "../types/Status.hpp"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace slide::core {

class ThreadPool
{
public:
  /** Zero selects hardware_concurrency (with a standards-compliant fallback to one). */
  explicit ThreadPool(unsigned workers = 0);
  ~ThreadPool();
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;
  ThreadPool(ThreadPool &&) = delete;
  ThreadPool &operator=(ThreadPool &&) = delete;

  unsigned workerCount() const { return static_cast<unsigned>(workers_.size()); }

  /**
   * Execute every index in [0,count) exactly once. The callback may return void or Status.
   * All tasks finish even after the first failure so result ownership is deterministic.
   */
  template <class Function>
  [[nodiscard]] slide::Status parallelFor(std::size_t count, Function &&function)
  {
    using Callable = std::remove_reference_t<Function>;
    Callable *callable = std::addressof(function);
    const auto invoke = [](void *context, std::size_t index) -> slide::Status {
      auto &target = *static_cast<Callable *>(context);
      if constexpr (std::is_same_v<std::invoke_result_t<Callable &, std::size_t>,
                                   slide::Status>)
        return target(index);
      else {
        static_assert(std::is_void_v<
                        std::invoke_result_t<Callable &, std::size_t>>,
                      "ThreadPool callbacks return void or slide::Status");
        target(index);
        return slide::Status::Success;
      }
    };
    return submit(count, callable, invoke);
  }

private:
  using Invoke = slide::Status (*)(void *, std::size_t);

  [[nodiscard]] slide::Status submit(std::size_t count,
                                     void *context,
                                     Invoke invoke);
  void workerLoop();

  std::vector<std::thread> workers_{};
  std::mutex mutex_{};
  std::condition_variable start_{};
  std::condition_variable done_{};
  std::atomic<std::size_t> next_{ 0 };
  std::atomic<int> first_status_{ static_cast<int>(slide::Status::Success) };
  std::size_t count_{};
  std::size_t workers_pending_{};
  std::uint64_t generation_{};
  void *context_{};
  Invoke invoke_{};
  bool running_{};
  bool stopping_{};
};

/** Fixed-order sum: result does not depend on pool scheduling or worker count. */
real_t fixedOrderSum(std::span<const real_t> values);

} // namespace slide::core

namespace slide::test {

struct ParallelisationDiagnostic
{
  unsigned logical_cores{};
  unsigned workers{};
  std::string backend{};
  double serial_seconds{};
  double parallel_seconds{};
  double speedup{};
  double checksum{};
  slide::Status status{ slide::Status::Success };
};

/** Measured smoke diagnostic; speedup is intentionally informational, never a gate. */
ParallelisationDiagnostic parallelisation(unsigned workers = 0,
                                          std::size_t values = 1U << 20U);

} // namespace slide::test
