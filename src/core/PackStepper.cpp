/**
 * @file PackStepper.cpp
 * @brief Transactional staggered electrical/thermal advance for compiled SPM packs.
 */

#include "PackStepper.hpp"

#include <algorithm>
#include <cstring>

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

  // A thermally uncoupled series of identical parallel-width layers preserves a repeated
  // lane brick exactly. Validate before sizing the batch steppers so their rollback and
  // cumulative workspaces also use only the unique brick.
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
      (void)batches[0]->setTrustedLanePeriod(width);
  }

  std::vector<TheveninBatchView> views;
  std::vector<EulerLegacy> steppers;
  std::vector<std::vector<real_t>> current_density;
  std::vector<std::size_t> checkpoint_offsets(batches.size() + 1);
  views.reserve(batches.size());
  steppers.reserve(batches.size());
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
    steppers.emplace_back(*batches[batch]);
    current_density.emplace_back(static_cast<std::size_t>(required_lanes[batch]));
    checkpoint_offsets[batch + 1] = checkpoint_offsets[batch] + ptr->state().size();
  }

  PackSolver solver;
  auto status = solver.configure(topology, views);
  if (status != slide::Status::Success)
    return status;

  topology_ = topology;
  batches_.assign(batches.begin(), batches.end());
  steppers_ = std::move(steppers);
  solver_ = std::move(solver);
  current_density_ = std::move(current_density);
  checkpoint_offsets_ = std::move(checkpoint_offsets);
  checkpoint_.assign(checkpoint_offsets_.back(), 0.0);
  cell_temperature_.assign(topology.cells.size(), 0.0);
  cell_external_heat_.assign(topology.cells.size(), 0.0);
  boundary_heat_.assign(topology.thermal.boundary_count, 0.0);
  configured_ = true;
  return slide::Status::Success;
}

void PackStepper::saveCheckpoint()
{
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    const auto state = batches_[batch]->state().raw();
    std::memcpy(checkpoint_.data() + checkpoint_offsets_[batch], state.data(), state.size_bytes());
  }
}

void PackStepper::restoreCheckpoint()
{
  for (std::size_t batch = 0; batch < batches_.size(); ++batch) {
    auto state = batches_[batch]->state().raw();
    std::memcpy(state.data(), checkpoint_.data() + checkpoint_offsets_[batch], state.size_bytes());
  }
  solver_.invalidate();
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
      status = steppers_[batch].step(*batches_[batch], current_density_[batch], time + static_cast<real_t>(substep) * dt, dt);
      if (status != slide::Status::Success)
        return fail(status);
    }
  return slide::Status::Success;
}

} // namespace slide::core
