/**
 * @file Recorder.cpp
 * @brief The snapshot recorder and its text/columnar sinks.
 *
 * Owns: `Recorder`'s configuration, snapshot capture, lazy voltage, and the CSV and Parquet
 * sinks. Implements PLAN.md §3.7. The binary recording format -- header, CRC, mmap, and the
 * `BinaryRecording` reader -- lives in `RecordingFormat.cpp`.
 * @surface internal
 */

#include "Recorder.hpp"
#include "detail/CheckedArithmetic.hpp"
#include "detail/RecordingFormatCommon.hpp"
#include "detail/SnapshotIndexing.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(SLIDE_WITH_ARROW)
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

namespace slide::core {

using detail::allocationFailureStatus;
using detail::checkedAdd;
using detail::checkedMultiply;

namespace {

  [[nodiscard]] std::vector<std::string> recorderColumnNames(int rows,
                                                             int lanes)
  {
    std::vector<std::string> names{ "accepted_step", "time_s" };
    for (int lane = 0; lane < lanes; ++lane)
      names.push_back(
        "current_density_lane" + std::to_string(lane) + "_A_m2");
    for (int row = 0; row < rows; ++row)
      for (int lane = 0; lane < lanes; ++lane)
        names.push_back(
          "state_r" + std::to_string(row) + "_lane"
          + std::to_string(lane));
    for (int lane = 0; lane < lanes; ++lane)
      names.push_back(
        "terminal_voltage_lane" + std::to_string(lane) + "_V");
    return names;
  }

  [[nodiscard]] std::span<const real_t> snapshotRow(
    const SnapshotView &snapshot,
    int row,
    int rows,
    int stride,
    int lanes) noexcept
  {
    assert(0 <= row && row < rows && 0 < lanes && lanes <= stride);
    return snapshot.state.subspan(
      static_cast<std::size_t>(row) * static_cast<std::size_t>(stride),
      static_cast<std::size_t>(lanes));
  }

} // namespace

slide::Status Recorder::configure(SpmBatch &batch, RecorderConfig config)
{
  if (!batch.valid() || config.cadence == 0 || config.capacity == 0
      || !(config.backpressure == BackpressurePolicy::stop
           || config.backpressure == BackpressurePolicy::thin))
    return slide::Status::Invalid_parameters;
  const std::size_t state_values = batch.state().size();
  const auto lanes = static_cast<std::size_t>(batch.n_lanes());
  std::size_t state_storage{}, current_storage{};
  if (!checkedMultiply(config.capacity, state_values, state_storage)
      || !checkedMultiply(config.capacity, lanes, current_storage)
      || config.capacity > std::vector<std::uint64_t>{}.max_size()
      || config.capacity > std::vector<real_t>{}.max_size()
      || state_storage > std::vector<real_t>{}.max_size()
      || current_storage > std::vector<real_t>{}.max_size())
    return slide::Status::Invalid_parameters;
  try {
    std::vector<std::uint64_t> steps(config.capacity);
    std::vector<real_t> times(config.capacity);
    std::vector<real_t> currents(current_storage);
    std::vector<real_t> states(state_storage);
    batch_ = &batch;
    config_ = config;
    rows_ = batch.state().n_rows();
    lanes_ = batch.state().n_lanes();
    stride_ = batch.state().stride();
    state_values_ = state_values;
    count_ = 0;
    thinned_ = 0;
    previous_step_ = 0;
    has_previous_step_ = false;
    accepted_steps_ = std::move(steps);
    times_ = std::move(times);
    current_density_ = std::move(currents);
    states_ = std::move(states);
  } catch (const std::bad_alloc &) {
    return allocationFailureStatus();
  } catch (const std::length_error &) {
    return allocationFailureStatus();
  }
  return slide::Status::Success;
}

slide::Status Recorder::record(std::uint64_t accepted_step,
                               std::span<const real_t>
                                 total_current_A)
{
  if (!configured() || total_current_A.size() != static_cast<std::size_t>(lanes_))
    return slide::Status::Invalid_parameters;
  if (accepted_step % config_.cadence != 0)
    return slide::Status::Success;
  const real_t time = batch_->state().at(batch_->layout().elapsed_time, 0, 0);
  if (!is_finite(time))
    return slide::Status::Invalid_parameters;
  for (const real_t current : total_current_A) {
    if (!is_finite(current)
        || !is_finite(current / batch_->electrode_area()))
      return slide::Status::Invalid_parameters;
  }
  if (has_previous_step_ && accepted_step <= previous_step_)
    return slide::Status::Invalid_parameters;
  previous_step_ = accepted_step;
  has_previous_step_ = true;
  if (count_ == config_.capacity) {
    if (config_.backpressure == BackpressurePolicy::thin) {
      ++thinned_;
      return slide::Status::Success;
    }
    return slide::Status::Numerical_failure;
  }
  accepted_steps_[count_] = accepted_step;
  times_[count_] = time;
  const std::size_t current_offset = count_ * static_cast<std::size_t>(lanes_);
  for (int lane = 0; lane < lanes_; ++lane)
    current_density_[current_offset + static_cast<std::size_t>(lane)] = total_current_A[static_cast<std::size_t>(lane)] / batch_->electrode_area();
  std::memcpy(states_.data() + count_ * state_values_,
              batch_->state().raw().data(),
              state_values_ * sizeof(real_t));
  ++count_;
  return slide::Status::Success;
}

SnapshotView Recorder::snapshot(std::size_t index) const
{
  assert(index < count_);
  return detail::snapshotView(index,
                              static_cast<std::size_t>(lanes_),
                              state_values_,
                              accepted_steps_,
                              times_,
                              current_density_,
                              states_);
}

void Recorder::clear()
{
  count_ = 0;
  thinned_ = 0;
  previous_step_ = 0;
  has_previous_step_ = false;
}

slide::Status Recorder::terminalVoltage(std::size_t index,
                                        std::span<real_t>
                                          output)
{
  if (!configured() || index >= count_
      || output.size() != static_cast<std::size_t>(lanes_))
    return slide::Status::Invalid_parameters;
  const auto recorded = snapshot(index);
  return batch_->terminalVoltageAt(recorded.state, recorded.current_density, output);
}

slide::Status Recorder::writeCsv(const std::filesystem::path &path)
try {
  if (!configured())
    return slide::Status::Invalid_parameters;
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    return slide::Status::Invalid_parameters;
  return writeCsvStream(output);
} catch (const std::bad_alloc &) {
  return allocationFailureStatus();
} catch (const std::length_error &) {
  return allocationFailureStatus();
}

slide::Status Recorder::writeCsvStream(std::ostream &output)
{
  const auto names = recorderColumnNames(rows_, lanes_);
  for (std::size_t index = 0; index < names.size(); ++index)
    output << (index == 0 ? "" : ",") << names[index];
  output << '\n'
         << std::setprecision(std::numeric_limits<real_t>::max_digits10)
         << std::scientific;
  std::vector<real_t> voltage(static_cast<std::size_t>(lanes_));
  for (std::size_t index = 0; index < count_; ++index) {
    const auto recorded = snapshot(index);
    const auto status = terminalVoltage(index, voltage);
    if (status != slide::Status::Success)
      return status;
    output << recorded.accepted_step << ',' << recorded.time;
    for (const real_t current : recorded.current_density)
      output << ',' << current;
    for (int row = 0; row < rows_; ++row)
      for (const real_t value :
           snapshotRow(recorded, row, rows_, stride_, lanes_))
        output << ',' << value;
    for (const real_t value : voltage)
      output << ',' << value;
    output << '\n';
  }
  return output.good() ? slide::Status::Success
                       : slide::Status::Numerical_failure;
}

slide::Status Recorder::writeParquet(const std::filesystem::path &path)
{
#if !defined(SLIDE_WITH_ARROW)
  (void)path;
  return slide::Status::NotImplementedYet;
#else
  try {
    if (!configured())
      return slide::Status::Invalid_parameters;
    const auto names = recorderColumnNames(rows_, lanes_);
    std::size_t column{};
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> columns;
    auto append_uint64 = [&](const std::string &name, auto value) {
      arrow::UInt64Builder builder;
      for (std::size_t index = 0; index < count_; ++index)
        if (!builder.Append(value(index)).ok())
          return false;
      std::shared_ptr<arrow::Array> array;
      if (!builder.Finish(&array).ok())
        return false;
      fields.push_back(arrow::field(name, arrow::uint64()));
      columns.push_back(std::move(array));
      return true;
    };
    auto append_double = [&](const std::string &name, auto value) {
      arrow::DoubleBuilder builder;
      for (std::size_t index = 0; index < count_; ++index)
        if (!builder.Append(value(index)).ok())
          return false;
      std::shared_ptr<arrow::Array> array;
      if (!builder.Finish(&array).ok())
        return false;
      fields.push_back(arrow::field(name, arrow::float64()));
      columns.push_back(std::move(array));
      return true;
    };
    if (!append_uint64(names[column++], [&](std::size_t i) { return accepted_steps_[i]; })
        || !append_double(names[column++], [&](std::size_t i) { return times_[i]; }))
      return slide::Status::Numerical_failure;
    for (int lane = 0; lane < lanes_; ++lane) {
      if (!append_double(names[column++],
                         [&](std::size_t i) {
                           return snapshot(i).current_density[static_cast<std::size_t>(lane)];
                         }))
        return slide::Status::Numerical_failure;
    }
    for (int row = 0; row < rows_; ++row)
      for (int lane = 0; lane < lanes_; ++lane)
        if (!append_double(names[column++],
                           [&](std::size_t i) {
                             return snapshotRow(
                               snapshot(i), row, rows_, stride_, lanes_)
                               [static_cast<std::size_t>(lane)];
                           }))
          return slide::Status::Numerical_failure;
    std::vector<std::vector<real_t>> voltage(
      static_cast<std::size_t>(lanes_), std::vector<real_t>(count_));
    std::vector<real_t> one_voltage(static_cast<std::size_t>(lanes_));
    for (std::size_t index = 0; index < count_; ++index) {
      if (terminalVoltage(index, one_voltage) != slide::Status::Success)
        return slide::Status::Numerical_failure;
      for (int lane = 0; lane < lanes_; ++lane)
        voltage[static_cast<std::size_t>(lane)][index] = one_voltage[static_cast<std::size_t>(lane)];
    }
    for (int lane = 0; lane < lanes_; ++lane)
      if (!append_double(names[column++],
                         [&](std::size_t i) {
                           return voltage[static_cast<std::size_t>(lane)][i];
                         }))
        return slide::Status::Numerical_failure;

    assert(column == names.size());
    const auto table = arrow::Table::Make(arrow::schema(std::move(fields)),
                                          std::move(columns));
    auto output_result = arrow::io::FileOutputStream::Open(path.string());
    if (!output_result.ok())
      return slide::Status::Invalid_parameters;
    auto output = *output_result;
    const auto status = parquet::arrow::WriteTable(
      *table, arrow::default_memory_pool(), output, std::max<std::int64_t>(1, static_cast<std::int64_t>(count_)));
    if (!status.ok() || !output->Close().ok())
      return slide::Status::Numerical_failure;
    return slide::Status::Success;
  } catch (const std::bad_alloc &) {
    return allocationFailureStatus();
  } catch (const std::length_error &) {
    return allocationFailureStatus();
  }
#endif
}

} // namespace slide::core
