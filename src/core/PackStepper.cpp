/**
 * @file PackStepper.cpp
 * @brief Transactional staggered electrical/thermal advance for compiled SPM packs.
 */

#include "PackStepper.hpp"

#include <algorithm>
#include <cstring>
#include <type_traits>

namespace slide::core {

slide::Status PackStepper::configure(
  const CompiledPackTopology &topology,
  std::span<SpmBatch *const>
    batches)
{
  if (topology.cells.empty() || batches.size() != topology.batch_archetypes.size())
    return slide::Status::Invalid_parameters;

  std::vector<int> required_lanes(batches.size());
  std::vector<unsigned char> thermal_batch(batches.size());
  for (const auto &cell : topology.cells) {
    if (cell.location.batch >= batches.size())
      return slide::Status::Invalid_parameters;
    required_lanes[cell.location.batch] = std::max(
      required_lanes[cell.location.batch], static_cast<int>(cell.location.lane + 1));
    thermal_batch[cell.location.batch] |= static_cast<unsigned char>(cell.thermal);
  }

  // A thermally uncoupled series of identical parallel-width layers may preserve a repeated
  // lane brick exactly. Derive the candidate here, but do not mutate the caller-owned batch
  // until all Status-returning configure checks have passed.
  int trusted_lane_period{};
  if (batches.size() == 1 && batches[0] != nullptr
      && topology.thermal.edges.empty()
      && topology.electrical.series_parallel_ladder
      && topology.electrical.ladder_offsets.size() > 2) {
    const int width = static_cast<int>(topology.electrical.ladder_offsets[1]
                                       - topology.electrical.ladder_offsets[0]);
    bool repeated = width > 0;
    for (std::size_t layer = 1;
         layer + 1 < topology.electrical.ladder_offsets.size() && repeated;
         ++layer)
      repeated = topology.electrical.ladder_offsets[layer + 1]
                   - topology.electrical.ladder_offsets[layer]
                 == static_cast<std::uint32_t>(width);
    if (repeated)
      trusted_lane_period = width;
  }

  std::vector<TheveninBatchView> views;
  std::vector<EulerLegacy> steppers;
  std::vector<ExponentialModal> exponential_steppers;
  std::vector<std::vector<real_t>> current_density;
  std::vector<std::size_t> checkpoint_offsets(batches.size() + 1);
  views.reserve(batches.size());
  steppers.reserve(batches.size());
  exponential_steppers.reserve(batches.size());
  current_density.reserve(batches.size());
  for (std::size_t batch = 0; batch < batches.size(); ++batch) {
    const auto *ptr = batches[batch];
    if (ptr == nullptr || !ptr->valid() || ptr->n_lanes() != required_lanes[batch])
      return slide::Status::Invalid_parameters;
    const bool pipeline_thermal = ptr->composition() == SpmComposition::thermal
                                  || ptr->composition() == SpmComposition::thermal_ageing;
    if (pipeline_thermal != (thermal_batch[batch] != 0))
      return slide::Status::Invalid_parameters;
    views.push_back(TheveninBatchView::bind(*batches[batch], required_lanes[batch]));
    current_density.emplace_back(static_cast<std::size_t>(required_lanes[batch]));
    checkpoint_offsets[batch + 1] = checkpoint_offsets[batch] + ptr->state().size();
  }

  PackSolver solver;
  auto status = solver.configure(topology, views);
  if (status != slide::Status::Success)
    return status;

  for (auto *batch : batches) {
    steppers.emplace_back();
    status = steppers.back().configure(*batch);
    if (status != slide::Status::Success)
      return status;
    exponential_steppers.emplace_back();
    status = exponential_steppers.back().configure(*batch);
    if (status != slide::Status::Success)
      return status;
  }

  CompiledPackTopology candidate_topology = topology;
  std::vector<SpmBatch *> candidate_batches(batches.begin(), batches.end());
  PackSolution solver_checkpoint_solution;
  solver_checkpoint_solution.cell_current.assign(topology.cells.size(), 0.0);
  solver_checkpoint_solution.node_voltage.assign(
    topology.electrical.node_count, 0.0);
  std::vector<real_t> checkpoint(checkpoint_offsets.back());
  std::vector<real_t> cell_temperature(topology.cells.size());
  std::vector<real_t> cell_external_heat(topology.cells.size());
  std::vector<real_t> boundary_heat(topology.thermal.boundary_count);
  std::vector<real_t> cell_external_heat_checkpoint(topology.cells.size());
  std::vector<real_t> boundary_heat_checkpoint(topology.thermal.boundary_count);

  if (trusted_lane_period > 0) {
    const auto period_status = batches[0]->setTrustedLanePeriod(trusted_lane_period);
    if (period_status == slide::Status::Success) {
      // The first configure allocated full-lane capacity. Rebinding to the smaller trusted
      // period only reuses those buffers, so the external batch mutation is the final
      // potentially throwing boundary before no-throw member publication.
      status = steppers[0].configure(*batches[0]);
      if (status != slide::Status::Success)
        return status;
    } else if (period_status != slide::Status::Invalid_states) {
      return period_status;
    }
  }

  static_assert(std::is_nothrow_move_assignable_v<CompiledPackTopology>);
  static_assert(std::is_nothrow_move_assignable_v<PackSolver>);
  topology_ = std::move(candidate_topology);
  batches_ = std::move(candidate_batches);
  steppers_ = std::move(steppers);
  exponential_steppers_ = std::move(exponential_steppers);
  solver_ = std::move(solver);
  current_density_ = std::move(current_density);
  checkpoint_offsets_ = std::move(checkpoint_offsets);
  checkpoint_ = std::move(checkpoint);
  cell_temperature_ = std::move(cell_temperature);
  cell_external_heat_ = std::move(cell_external_heat);
  boundary_heat_ = std::move(boundary_heat);
  solver_checkpoint_solution_ = std::move(solver_checkpoint_solution);
  solver_checkpoint_diagnostics_ = {};
  cell_external_heat_checkpoint_ = std::move(cell_external_heat_checkpoint);
  boundary_heat_checkpoint_ = std::move(boundary_heat_checkpoint);
  solver_checkpoint_has_solution_ = false;
  configured_ = true;
  return slide::Status::Success;
}

void PackStepper::saveCheckpoint()
{
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    const auto state = batches_[batch]->state().raw();
    std::memcpy(checkpoint_.data() + checkpoint_offsets_[batch], state.data(), state.size_bytes());
  }
  solver_checkpoint_solution_.cell_current = solver_.solution_.cell_current;
  solver_checkpoint_solution_.node_voltage = solver_.solution_.node_voltage;
  solver_checkpoint_solution_.terminal_voltage = solver_.solution_.terminal_voltage;
  solver_checkpoint_diagnostics_ = solver_.diagnostics_;
  solver_checkpoint_has_solution_ = solver_.has_solution_;
  std::copy(cell_external_heat_.begin(),
            cell_external_heat_.end(),
            cell_external_heat_checkpoint_.begin());
  std::copy(boundary_heat_.begin(),
            boundary_heat_.end(),
            boundary_heat_checkpoint_.begin());
}

void PackStepper::restoreCheckpoint()
{
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    auto state = batches_[batch]->state().raw();
    std::memcpy(state.data(), checkpoint_.data() + checkpoint_offsets_[batch], state.size_bytes());
  }
  solver_.solution_.cell_current = solver_checkpoint_solution_.cell_current;
  solver_.solution_.node_voltage = solver_checkpoint_solution_.node_voltage;
  solver_.solution_.terminal_voltage = solver_checkpoint_solution_.terminal_voltage;
  solver_.diagnostics_ = solver_checkpoint_diagnostics_;
  solver_.has_solution_ = solver_checkpoint_has_solution_;
  solver_.workspace_.invalidate();
  std::copy(cell_external_heat_checkpoint_.begin(),
            cell_external_heat_checkpoint_.end(),
            cell_external_heat_.begin());
  std::copy(boundary_heat_checkpoint_.begin(),
            boundary_heat_checkpoint_.end(),
            boundary_heat_.begin());
}

slide::Status PackStepper::checkpoint(std::span<real_t> destination) const
{
  if (!configured_ || destination.size() != checkpoint_.size())
    return slide::Status::Invalid_parameters;
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    const auto state = batches_[batch]->state().raw();
    std::memcpy(destination.data() + checkpoint_offsets_[batch], state.data(), state.size_bytes());
  }
  return slide::Status::Success;
}

slide::Status PackStepper::restore(std::span<const real_t> source)
{
  if (!configured_ || source.size() != checkpoint_.size())
    return slide::Status::Invalid_parameters;
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    auto state = batches_[batch]->state().raw();
    std::memcpy(state.data(), source.data() + checkpoint_offsets_[batch], state.size_bytes());
  }
  solver_.invalidate();
  return slide::Status::Success;
}

slide::Status PackStepper::step(
  real_t applied_current,
  real_t time,
  real_t dt,
  std::span<const real_t>
    boundary_temperature,
  PackSolveMode mode,
  real_t current_tolerance,
  int substeps)
{
  return stepImpl(applied_current, time, dt, boundary_temperature, mode, current_tolerance, substeps, false);
}

slide::Status PackStepper::stepExponential(
  real_t applied_current,
  real_t time,
  real_t dt,
  std::span<const real_t>
    boundary_temperature,
  PackSolveMode mode,
  real_t current_tolerance)
{
  return stepImpl(applied_current, time, dt, boundary_temperature, mode, current_tolerance, 1, true);
}

slide::Status PackStepper::stepImpl(
  real_t applied_current,
  real_t time,
  real_t dt,
  std::span<const real_t>
    boundary_temperature,
  PackSolveMode mode,
  real_t current_tolerance,
  int substeps,
  bool exponential)
{
  if (!configured_ || !is_finite(applied_current) || !is_finite(time)
      || !is_finite(dt) || !(dt > 0.0) || substeps <= 0
      || boundary_temperature.size() != topology_.thermal.boundary_count)
    return slide::Status::Invalid_parameters;

  saveCheckpoint();
  auto fail = [&](slide::Status status) {
    restoreCheckpoint();
    return status;
  };

  auto status = solver_.solve(applied_current, mode, current_tolerance);
  if (status != slide::Status::Success)
    return fail(status);

  for (std::size_t cell = 0; cell < topology_.cells.size(); ++cell) {
    const auto location = topology_.cells[cell].location;
    auto &batch = *batches_[location.batch];
    current_density_[location.batch][location.lane] = solver_.solution().cell_current[cell]
                                                      / batch.electrode_area();
    cell_temperature_[cell] = batch.state().at(batch.layout().spm.temperature,
                                               0,
                                               static_cast<int>(location.lane));
  }

  status = topology_.thermal.assemble(cell_temperature_, boundary_temperature, cell_external_heat_, boundary_heat_);
  if (status != slide::Status::Success)
    return fail(status);
  for (std::size_t cell = 0; cell < topology_.cells.size(); ++cell) {
    if (!topology_.cells[cell].thermal)
      continue;
    const auto location = topology_.cells[cell].location;
    auto &batch = *batches_[location.batch];
    batch.state().at(batch.layout().thermal.external_heat_flow,
                     0,
                     static_cast<int>(location.lane)) = cell_external_heat_[cell];
  }

  for (int substep = 0; substep < substeps; ++substep)
    for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
      status = exponential
                 ? exponential_steppers_[batch].step(
                     *batches_[batch], current_density_[batch], time + static_cast<real_t>(substep) * dt, dt)
                 : steppers_[batch].step(
                     *batches_[batch], current_density_[batch], time + static_cast<real_t>(substep) * dt, dt);
      if (status != slide::Status::Success)
        return fail(status);
    }
  return slide::Status::Success;
}

} // namespace slide::core
