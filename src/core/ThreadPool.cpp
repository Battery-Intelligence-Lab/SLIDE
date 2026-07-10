/**
 * @file ThreadPool.cpp
 * @brief Persistent std::thread runtime implementation.
 */

#include "ThreadPool.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace slide::core {

ThreadPool::ThreadPool(unsigned workers)
{
  if (workers == 0)
    workers = std::thread::hardware_concurrency();
  if (workers == 0)
    workers = 1;
  workers_.reserve(workers);
  try {
    for (unsigned i = 0; i < workers; ++i)
      workers_.emplace_back([this] { workerLoop(); });
  } catch (...) {
    {
      const std::lock_guard lock{ mutex_ };
      stopping_ = true;
    }
    start_.notify_all();
    for (auto &worker : workers_)
      if (worker.joinable())
        worker.join();
    throw;
  }
}

ThreadPool::~ThreadPool()
{
  {
    const std::lock_guard lock{ mutex_ };
    stopping_ = true;
  }
  start_.notify_all();
  for (auto &worker : workers_)
    if (worker.joinable())
      worker.join();
}

BatchExecutor::BatchExecutor(BatchExecutor &&other) noexcept
  : pool_{ std::move(other.pool_) },
    batches_{ std::exchange(other.batches_, 0) },
    workers_{ std::exchange(other.workers_, 0) },
    configured_{ std::exchange(other.configured_, false) }
{}

BatchExecutor &BatchExecutor::operator=(BatchExecutor &&other) noexcept
{
  if (this != &other) {
    pool_ = std::move(other.pool_);
    batches_ = std::exchange(other.batches_, 0);
    workers_ = std::exchange(other.workers_, 0);
    configured_ = std::exchange(other.configured_, false);
  }
  return *this;
}

slide::Status BatchExecutor::configure(std::size_t batches, unsigned workers)
{
  if (batches == 0)
    return slide::Status::Invalid_parameters;
  if (workers == 0)
    workers = std::thread::hardware_concurrency();
  if (workers == 0)
    workers = 1;
  workers = static_cast<unsigned>(
    std::min<std::size_t>(workers, batches));

  std::unique_ptr<ThreadPool> pool;
  try {
    if (workers > 1)
      pool = std::make_unique<ThreadPool>(workers);
  } catch (...) {
    return slide::Status::Numerical_failure;
  }
  pool_ = std::move(pool);
  batches_ = batches;
  workers_ = workers;
  configured_ = true;
  return slide::Status::Success;
}

slide::Status ThreadPool::submit(std::size_t count,
                                 void *context,
                                 Invoke invoke)
{
  if (count == 0)
    return slide::Status::Success;
  if (context == nullptr || invoke == nullptr || workers_.empty()
      || count > std::numeric_limits<std::size_t>::max() - workers_.size())
    return slide::Status::Invalid_parameters;

  std::unique_lock lock{ mutex_ };
  if (running_ || stopping_)
    return slide::Status::Invalid_states;
  running_ = true;
  count_ = count;
  context_ = context;
  invoke_ = invoke;
  next_.store(0, std::memory_order_relaxed);
  failure_index_ = count;
  first_status_ = slide::Status::Success;
  workers_pending_ = workers_.size();
  ++generation_;
  lock.unlock();
  start_.notify_all();

  lock.lock();
  done_.wait(lock, [this] { return !running_; });
  return first_status_;
}

void ThreadPool::workerLoop()
{
  std::uint64_t observed_generation{};
  while (true) {
    std::unique_lock lock{ mutex_ };
    start_.wait(lock, [&] {
      return stopping_ || generation_ != observed_generation;
    });
    if (stopping_)
      return;
    const auto generation = generation_;
    const auto count = count_;
    auto *context = context_;
    const auto invoke = invoke_;
    lock.unlock();

    for (std::size_t index = next_.fetch_add(1, std::memory_order_relaxed);
         index < count;
         index = next_.fetch_add(1, std::memory_order_relaxed)) {
      slide::Status status = slide::Status::Unknown_problem;
      try {
        status = invoke(context, index);
      } catch (...) {
        status = slide::Status::Unknown_problem;
      }
      if (status != slide::Status::Success) {
        const std::lock_guard failure_lock{ mutex_ };
        if (index < failure_index_) {
          failure_index_ = index;
          first_status_ = status;
        }
      }
    }

    lock.lock();
    observed_generation = generation;
    if (--workers_pending_ == 0) {
      running_ = false;
      context_ = nullptr;
      invoke_ = nullptr;
      done_.notify_one();
    }
  }
}

#if defined(_MSC_VER) && !defined(__clang__)
#pragma float_control(precise, on, push)
#endif
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-fast-math")))
#endif
real_t fixedOrderSum(std::span<const real_t> values)
{
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
  volatile real_t sum{};
  for (const real_t value : values) {
    volatile real_t updated = sum + value;
    sum = updated;
  }
  return sum;
}
#if defined(_MSC_VER) && !defined(__clang__)
#pragma float_control(pop)
#endif

} // namespace slide::core

namespace slide::test {
namespace {

  template <class Function>
  double elapsed(Function &&function)
  {
    const auto start = std::chrono::steady_clock::now();
    function();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
  }

} // namespace

ParallelisationDiagnostic parallelisation(unsigned workers, std::size_t values)
{
  ParallelisationDiagnostic diagnostic;
  diagnostic.logical_cores = std::thread::hardware_concurrency();
  diagnostic.backend = "std::thread persistent pool";
  if (values == 0) {
    diagnostic.status = slide::Status::Invalid_parameters;
    return diagnostic;
  }
  try {
    slide::core::ThreadPool pool{ workers };
    diagnostic.workers = pool.workerCount();
    std::vector<double> input(values), serial(values), parallel(values);
    for (std::size_t i = 0; i < values; ++i)
      input[i] = 1.0 + static_cast<double>(i % 1009) / 1009.0;
    const auto kernel = [&](std::span<double> output, std::size_t index) {
      double value = input[index];
      for (int iteration = 0; iteration < 24; ++iteration)
        value = std::sqrt(value + 0.5) + 0.125 * value;
      output[index] = value;
    };
    diagnostic.serial_seconds = elapsed([&] {
      for (std::size_t i = 0; i < values; ++i)
        kernel(serial, i);
    });
    diagnostic.parallel_seconds = elapsed([&] {
      diagnostic.status = pool.parallelFor(values, [&](std::size_t i) {
        kernel(parallel, i);
      });
    });
    if (diagnostic.status != slide::Status::Success || serial != parallel) {
      diagnostic.status = slide::Status::Numerical_failure;
      return diagnostic;
    }
    diagnostic.checksum = slide::core::fixedOrderSum(parallel);
    diagnostic.speedup = diagnostic.parallel_seconds > 0.0
                           ? diagnostic.serial_seconds
                               / diagnostic.parallel_seconds
                           : std::numeric_limits<double>::max();
  } catch (...) {
    diagnostic.status = slide::Status::Unknown_problem;
  }
  return diagnostic;
}

} // namespace slide::test
