/**
 * @file CudaSpmData.hpp
 * @brief Private POD boundary between the dependency-light host facade and CUDA runtime.
 */

#pragma once

#include <cstddef>
#include <span>

namespace slide::core::cuda_detail {

inline constexpr int maximum_nch = 12;

struct CurveRef
{
  int x_offset{};
  int y_offset{};
  int count{};
};

struct ElectrodeData
{
  double cs_max{};
  double reaction_rate_ref{};
  double reaction_activation{};
};

/** All integers are row indices into the SoA state arena. */
struct KernelParams
{
  int nch{};
  int lanes{};
  int stride{};
  int rows{};
  int z_row[2]{};
  int temperature_row{};
  int sei_thickness_row{};
  int diffusion_row[2]{};
  int thickness_row[2]{};
  int specific_area_row[2]{};
  int specific_resistance_row[2]{};
  int collector_resistance_row{};
  int elapsed_time_row{};
  int charge_throughput_row{};
  int energy_throughput_row{};
  double A[2][maximum_nch]{};
  double B[2][maximum_nch]{};
  double surface_C[2][maximum_nch]{};
  double surface_D[2]{};
  double diffusion_activation[2]{};
  ElectrodeData electrode[2]{};
  CurveRef electrode_ocv[2]{};
  CurveRef total_entropic{};
  double F{ 96487.0 };
  double Rg{ 8.314 };
  double n{ 1.0 };
  double electrolyte_concentration{ 1000.0 };
  double reference_temperature{ 298.15 };
  double electrode_area{};
  double sei_resistivity_area{};
};

struct Runtime;

enum class RuntimeResult : int {
  success,
  invalid_parameter,
  cuda_failure,
  invalid_state
};

int deviceCount() noexcept;
RuntimeResult create(const KernelParams &params,
                     std::span<const double> curves,
                     std::span<const double> state,
                     Runtime *&runtime) noexcept;
void destroy(Runtime *runtime) noexcept;
RuntimeResult uploadState(Runtime *runtime,
                          std::span<const double> state) noexcept;
RuntimeResult downloadState(Runtime *runtime,
                            std::span<double> state) noexcept;
RuntimeResult launchStep(Runtime *runtime,
                         std::span<const double> current_density,
                         double time,
                         double dt) noexcept;
RuntimeResult synchronize(Runtime *runtime,
                          std::span<double> terminal_voltage) noexcept;
RuntimeResult checkpoint(Runtime *runtime) noexcept;
RuntimeResult restore(Runtime *runtime) noexcept;
std::size_t deviceBytes(const Runtime *runtime) noexcept;
std::size_t deviceAllocations(const Runtime *runtime) noexcept;
std::size_t deviceWideSynchronizations(const Runtime *runtime) noexcept;

} // namespace slide::core::cuda_detail
