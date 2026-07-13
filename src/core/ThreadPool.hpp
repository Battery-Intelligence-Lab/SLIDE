/**
 * @file ThreadPool.hpp
 * @brief Persistent dependency-free batch task runtime and deterministic reductions.
 * @surface api
 */

#pragma once

#include "Numeric.hpp"
#include "../types/Status.hpp"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace slide::core {

namespace detail {

  struct BatchExecutorTestAccess;

  template <class Callable>
  slide::Status invokeIndexed(Callable &target, std::size_t index)
  {
    using Result = std::invoke_result_t<Callable &, std::size_t>;
    if constexpr (std::is_same_v<std::remove_cvref_t<Result>, slide::Status>)
      return std::invoke(target, index);
    else {
      static_assert(std::is_void_v<Result>,
                    "indexed callbacks return void or slide::Status");
      std::invoke(target, index);
      return slide::Status::Success;
    }
  }

} // namespace detail

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
   * All tasks finish after a failure; the status from the lowest failing index is returned.
   * A count whose worker stop sentinels would overflow size_t is rejected.
   */
  template <class Function>
  [[nodiscard]] slide::Status parallelFor(std::size_t count, Function &&function)
  {
    if (count == 0)
      return slide::Status::Success;
    using Callable = std::remove_reference_t<Function>;
    if constexpr (std::is_pointer_v<std::remove_cv_t<Callable>>)
      if (function == nullptr)
        return slide::Status::Invalid_parameters;
    struct Context
    {
      Callable *callable;
    } context{ std::addressof(function) };
    const auto invoke = [](void *context, std::size_t index) -> slide::Status {
      auto &target = *static_cast<Context *>(context)->callable;
      return detail::invokeIndexed(target, index);
    };
    return submit(count, std::addressof(context), invoke);
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
  std::size_t count_{};
  std::size_t failure_index_{};
  std::size_t workers_pending_{};
  std::uint64_t generation_{};
  void *context_{};
  Invoke invoke_{};
  slide::Status first_status_{ slide::Status::Success };
  bool running_{};
  bool stopping_{};
};

/**
 * Movable production execution context for independent archetype batches.
 * Configuration creates at most one persistent pool; a selected single worker
 * executes inline and is reported explicitly rather than silently falling back.
 */
class BatchExecutor
{
public:
  BatchExecutor() = default;
  BatchExecutor(BatchExecutor &&other) noexcept;
  BatchExecutor &operator=(BatchExecutor &&other) noexcept;
  BatchExecutor(const BatchExecutor &) = delete;
  BatchExecutor &operator=(const BatchExecutor &) = delete;

  [[nodiscard]] slide::Status configure(std::size_t batches,
                                        unsigned workers = 0);
  unsigned workerCount() const { return configured_ ? workers_ : 0; }

  template <class Function>
  [[nodiscard]] slide::Status parallelFor(std::size_t count,
                                          Function &&function)
  {
    if (count == 0)
      return slide::Status::Success;
    if (!configured_ || count > batches_)
      return slide::Status::Invalid_parameters;
    using Callable = std::remove_reference_t<Function>;
    if constexpr (std::is_pointer_v<std::remove_cv_t<Callable>>)
      if (function == nullptr)
        return slide::Status::Invalid_parameters;
    if (pool_ != nullptr)
      return pool_->parallelFor(count, function);

    slide::Status first = slide::Status::Success;
    for (std::size_t index = 0; index < count; ++index) {
      slide::Status status = slide::Status::Unknown_problem;
      try {
        status = detail::invokeIndexed(function, index);
      } catch (...) {
        status = slide::Status::Unknown_problem;
      }
      if (first == slide::Status::Success && status != slide::Status::Success)
        first = status;
    }
    return first;
  }

private:
  using PoolFactory = std::unique_ptr<ThreadPool> (*)(unsigned);
  static std::unique_ptr<ThreadPool> makePool(unsigned workers);

  std::unique_ptr<ThreadPool> pool_{};
  PoolFactory pool_factory_{ &BatchExecutor::makePool };
  std::size_t batches_{};
  unsigned workers_{};
  bool configured_{};

  friend struct detail::BatchExecutorTestAccess;
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
