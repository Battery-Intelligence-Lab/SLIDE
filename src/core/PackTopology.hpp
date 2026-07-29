/**
 * @file PackTopology.hpp
 * @brief Value-semantic pack combinators and compiled electrical/thermal topology.
 * Owns D-19 pack authoring/compilation plus canonical electrical and D-21 thermal outputs.
 * Implements PLAN.md section 3.4 (D-19 and D-21).
 * Cold: description authoring and compilePackDescription().
 * Hot: CompiledThermalGraph::assemble(), once per PackStepper step attempt.
 * @surface api
 */

#pragma once

#include "Numeric.hpp"
#include "../types/Status.hpp"

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace slide::core {

struct PackCellSpec
{
  std::string archetype{ "cell" };
  bool thermal{};
};

struct PackLink
{
  real_t resistance{}; //!< inserted between series replicas / ahead of parallel replicas [ohm]
};

enum class PackNodeKind : unsigned char { cell,
                                          series,
                                          parallel };

struct PackNode
{
  PackNodeKind kind{ PackNodeKind::cell };
  PackCellSpec cell{};
  PackLink link{};
  std::vector<PackNode> children{};
};

[[nodiscard]] PackNode cell(PackCellSpec spec = {});
[[nodiscard]] PackNode series(int count, const PackNode &child, PackLink link = {});
[[nodiscard]] PackNode parallel(int count, const PackNode &child, PackLink link = {});
[[nodiscard]] PackNode series(std::vector<PackNode> children, PackLink link = {});
[[nodiscard]] PackNode parallel(std::vector<PackNode> children, PackLink link = {});

struct ThermalBoundarySpec
{
  std::string name{};
};

struct ThermalLinkSpec
{
  std::string endpoint_a{}; //!< compiled cell path or boundary name
  std::string endpoint_b{};
  real_t conductance{}; //!< W/K
};

struct PackDescription
{
  PackNode root{};
  std::vector<ThermalBoundarySpec> thermal_boundaries{};
  std::vector<ThermalLinkSpec> thermal_links{};
};

enum class ElectricalBranchKind : unsigned char { cell,
                                                  resistor };

struct CompiledElectricalBranch
{
  std::uint32_t node_positive{};
  std::uint32_t node_negative{};
  ElectricalBranchKind kind{ ElectricalBranchKind::cell };
  std::uint32_t cell{}; //!< valid only for cell branches
  real_t resistance{};  //!< valid only for resistor branches
};

struct BatchLaneLocation
{
  std::uint32_t batch{};
  std::uint32_t lane{};
};

struct CompiledCell
{
  std::string path{};
  std::string archetype{};
  BatchLaneLocation location{};
  bool thermal{};
};

struct CompiledElectricalNetlist
{
  std::uint32_t node_count{};
  std::uint32_t terminal_positive{};
  std::uint32_t terminal_negative{ 1 };
  std::vector<CompiledElectricalBranch> branches{};
  std::vector<std::pair<std::uint32_t, std::uint32_t>> nodal_sparsity{};
  std::vector<std::uint32_t> ladder_offsets{}; //!< layer -> range in ladder_cells
  std::vector<std::uint32_t> ladder_cells{};   //!< cell indices, positive-to-negative
  std::vector<std::uint32_t> ladder_nodes{};   //!< positive-to-negative path nodes
  bool connected{};
  bool index1_candidate{};
  bool series_parallel_ladder{};
};

struct ThermalEdge
{
  std::uint32_t low{};
  std::uint32_t high{};
  real_t conductance{};
};

struct ThermalIncident
{
  std::uint32_t edge{};
  std::int32_t sign{};
};

struct CompiledThermalGraph
{
  std::uint32_t cell_count{};
  std::uint32_t boundary_count{};
  std::vector<ThermalEdge> edges{};
  std::vector<std::uint32_t> offsets{};
  std::vector<ThermalIncident> incidents{};
  std::vector<real_t> edge_flux{};                   //!< last successfully assembled edge fluxes
  std::vector<real_t> trial_edge_flux{};             //!< preallocated transactional scratch
  std::vector<real_t> trial_endpoint_heat{};         //!< preallocated transactional scratch
  std::vector<unsigned char> trial_edge_incidence{}; //!< low/high membership bits

  [[nodiscard]] slide::Status assemble(std::span<const real_t> cell_temperature,
                                       std::span<const real_t> boundary_temperature,
                                       std::span<real_t> q_ext,
                                       std::span<real_t> boundary_heat) noexcept;
};

struct CompiledPackTopology
{
  std::vector<CompiledCell> cells{};
  std::vector<std::string> batch_archetypes{};
  CompiledElectricalNetlist electrical{};
  CompiledThermalGraph thermal{};
};

/** Cold compile is atomic: output is unchanged on validation failure. */
[[nodiscard]] slide::Status compilePackDescription(const PackDescription &description,
                                                   CompiledPackTopology &output);

} // namespace slide::core
