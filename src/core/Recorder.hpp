/**
 * @file Recorder.hpp
 * @brief Allocation-free snapshot recorder and portable CSV/mmap sinks.
 * @surface api
 */

#pragma once

#include "SpmFactory.hpp"

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <memory>
#include <span>
#include <vector>

namespace slide::core {

namespace detail {

  struct RecorderTestAccess;

  struct BinaryRecordingLayout
  {
    std::uint64_t current_bytes{};
    std::uint64_t state_bytes{};
    std::uint64_t record_bytes{};
    std::uint64_t table_bytes{};
    std::uint64_t data_offset{};
    std::uint64_t file_size{};
  };

  /** Pure checked layout used by the mmap writer and boundary tests. */
  [[nodiscard]] slide::Status binaryRecordingLayout(
    std::uint64_t lanes,
    std::uint64_t state_values,
    std::uint64_t snapshots,
    BinaryRecordingLayout &output);

} // namespace detail

enum class BackpressurePolicy : unsigned char {
  stop, //!< return a failure when the preallocated ring is full
  thin  //!< keep simulating and count every omitted cadence point
};

struct RecorderConfig
{
  std::size_t cadence{ 1 };
  std::size_t capacity{ 1024 };
  BackpressurePolicy backpressure{ BackpressurePolicy::stop };
};

struct SnapshotView
{
  std::uint64_t accepted_step{};
  real_t time{};
  std::span<const real_t> current_density{};
  std::span<const real_t> state{};
};

/**
 * A fixed-capacity, deferred-I/O recorder. configure() allocates every byte; record()
 * only copies into the next ring slot and therefore performs no heap or file work.
 */
class Recorder
{
public:
  [[nodiscard]] slide::Status configure(SpmBatch &batch,
                                        RecorderConfig config = {});
  [[nodiscard]] slide::Status record(
    std::uint64_t accepted_step,
    std::span<const real_t> total_current_A);

  bool configured() const { return batch_ != nullptr; }
  std::size_t size() const { return count_; }
  std::size_t capacity() const { return config_.capacity; }
  std::uint64_t thinnedSnapshots() const { return thinned_; }
  int n_rows() const { return rows_; }
  int n_lanes() const { return lanes_; }
  int stride() const { return stride_; }
  std::size_t stateValues() const { return state_values_; }
  SnapshotView snapshot(std::size_t index) const;
  void clear();

  /** Uses the batch's production observable kernel against recorded state bytes. */
  [[nodiscard]] slide::Status terminalVoltage(
    std::size_t index,
    std::span<real_t> output);

  [[nodiscard]] slide::Status writeCsv(
    const std::filesystem::path &path);
  [[nodiscard]] slide::Status writeBinary(
    const std::filesystem::path &path) const;
  /** Optional Arrow/Parquet sink; returns NotImplementedYet when disabled. */
  [[nodiscard]] slide::Status writeParquet(
    const std::filesystem::path &path);

private:
  using MappingFlush = bool (*)(void *mapping);

  slide::Status writeCsvStream(std::ostream &output);
  static bool flushMapping(void *mapping);

  SpmBatch *batch_{};
  RecorderConfig config_{};
  int rows_{};
  int lanes_{};
  int stride_{};
  std::size_t state_values_{};
  std::size_t count_{};
  std::uint64_t thinned_{};
  std::uint64_t previous_step_{};
  bool has_previous_step_{};
  std::vector<std::uint64_t> accepted_steps_{};
  std::vector<real_t> times_{};
  std::vector<real_t> current_density_{};
  std::vector<real_t> states_{};
  MappingFlush mapping_flush_{ &Recorder::flushMapping };

  friend struct detail::RecorderTestAccess;
};

/** Read-only hardened view over a SLIDE recording mapped directly from disk. */
class BinaryRecording
{
public:
  BinaryRecording();
  ~BinaryRecording();
  BinaryRecording(const BinaryRecording &) = delete;
  BinaryRecording &operator=(const BinaryRecording &) = delete;
  BinaryRecording(BinaryRecording &&) noexcept;
  BinaryRecording &operator=(BinaryRecording &&) noexcept;

  [[nodiscard]] slide::Status open(const std::filesystem::path &path);
  void close();
  bool valid() const;
  std::size_t size() const { return snapshots_; }
  int n_rows() const { return rows_; }
  int n_lanes() const { return lanes_; }
  int stride() const { return stride_; }
  SnapshotView snapshot(std::size_t index) const;

  struct Mapping;

private:
  std::unique_ptr<Mapping> mapping_{};
  std::size_t snapshots_{};
  int rows_{};
  int lanes_{};
  int stride_{};
  std::uint64_t table_offset_{};
  std::size_t state_values_{};
};

} // namespace slide::core
