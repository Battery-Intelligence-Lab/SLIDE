/**
 * @file AgeingKernel.hpp
 * @brief Shared zero-overhead scaffolding for scalar-generic ageing kernels.
 *
 * M0.6 / 9C-2 keeps SEI, surface-crack, LAM, and lithium-plating equations in
 * their named mechanism headers. This file owns only their common model-mask,
 * checked field-major scratch, clearing, and traversal idiom. Traversal order is
 * part of the numerical contract: enabled models are visited in ascending order,
 * and every model visits lanes in ascending order.
 */

#pragma once

#include "../types/Status.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

// Besides removing hot-path call overhead, forced inlining preserves the caller's
// opaque-reference finite checks under Clang fast-math. A merely-inline traversal
// caused the surface-crack overflow regression to return Success in Release.
#if defined(_MSC_VER) && !defined(__clang__)
#define SLIDE_AGEING_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SLIDE_AGEING_FORCE_INLINE inline __attribute__((always_inline))
#else
#define SLIDE_AGEING_FORCE_INLINE inline
#endif

namespace slide::core {

template <unsigned ModelCount>
constexpr std::uint8_t ageing_model_bit(unsigned model) noexcept
{
  static_assert(ModelCount > 0 && ModelCount <= 8);
  return model >= 1 && model <= ModelCount
           ? static_cast<std::uint8_t>(std::uint16_t{ 1 } << (model - 1))
           : std::uint8_t{};
}

template <unsigned ModelCount>
constexpr std::uint8_t ageing_model_mask() noexcept
{
  static_assert(ModelCount > 0 && ModelCount <= 8);
  return static_cast<std::uint8_t>((std::uint16_t{ 1 } << ModelCount) - 1);
}

template <unsigned ModelCount>
constexpr bool valid_ageing_model_mask(std::uint8_t mask) noexcept
{
  return mask != 0
         && (mask & static_cast<std::uint8_t>(~ageing_model_mask<ModelCount>())) == 0;
}

template <unsigned ModelCount>
constexpr bool valid_optional_ageing_model_mask(std::uint8_t mask) noexcept
{
  return (mask & static_cast<std::uint8_t>(~ageing_model_mask<ModelCount>())) == 0;
}

namespace detail {

  template <class Real, std::size_t FieldCount>
  class AgeingScratchStorage
  {
  public:
    static_assert(FieldCount > 0);

    explicit AgeingScratchStorage(int n_lanes)
      : n_lanes_{ checked_lane_count(n_lanes) },
        storage_(required_elements(n_lanes))
    {}

    std::span<Real> field(std::size_t index)
    {
      assert(index < FieldCount);
      const auto lanes = static_cast<std::size_t>(n_lanes_);
      return std::span<Real>{ storage_ }.subspan(index * lanes, lanes);
    }

    int n_lanes() const noexcept { return n_lanes_; }
    static constexpr std::size_t field_count = FieldCount;

  private:
    static int checked_lane_count(int n_lanes)
    {
      if (n_lanes <= 0)
        throw std::invalid_argument{ "ageing scratch requires at least one lane" };
      return n_lanes;
    }

    static std::size_t required_elements(int n_lanes)
    {
      if (n_lanes <= 0)
        throw std::invalid_argument{ "ageing scratch requires at least one lane" };
      const auto lanes = static_cast<std::size_t>(n_lanes);
      if (lanes > std::numeric_limits<std::size_t>::max() / FieldCount)
        throw std::length_error{ "ageing scratch extent is not representable" };
      return lanes * FieldCount;
    }

    int n_lanes_{};
    std::vector<Real> storage_{};
  };

  template <class Real, std::size_t FieldCount>
  SLIDE_AGEING_FORCE_INLINE void clear_ageing_fields(
    int n_lanes,
    const std::array<std::span<Real>, FieldCount> &fields)
  {
    assert(n_lanes > 0);
    [[maybe_unused]] const auto lanes = static_cast<std::size_t>(n_lanes);
    for (const auto field : fields) {
      assert(field.size() == lanes);
      std::fill(field.begin(), field.end(), Real{});
    }
  }

  template <class Body>
  SLIDE_AGEING_FORCE_INLINE void for_each_ageing_lane(int n_lanes, Body &&body)
  {
    assert(n_lanes > 0);
    for (int lane = 0; lane < n_lanes; ++lane)
      body(lane);
  }

  template <class Body>
  [[nodiscard]] SLIDE_AGEING_FORCE_INLINE slide::Status for_each_ageing_lane_while_success(
    int n_lanes,
    Body &&body)
  {
    assert(n_lanes > 0);
    for (int lane = 0; lane < n_lanes; ++lane) {
      const auto status = body(lane);
      if (status != slide::Status::Success)
        return status;
    }
    return slide::Status::Success;
  }

  template <unsigned ModelCount, class Body>
  [[nodiscard]] SLIDE_AGEING_FORCE_INLINE slide::Status for_each_enabled_ageing_model_lane(
    std::uint8_t mask,
    int n_lanes,
    Body &&body)
  {
    assert(valid_ageing_model_mask<ModelCount>(mask) && n_lanes > 0);
    for (unsigned model = 1; model <= ModelCount; ++model) {
      if ((mask & ageing_model_bit<ModelCount>(model)) == 0)
        continue;
      for (int lane = 0; lane < n_lanes; ++lane) {
        const auto status = body(model, lane);
        if (status != slide::Status::Success)
          return status;
      }
    }
    return slide::Status::Success;
  }

  /** The common pipeline transaction: disabled -> no-op; compute -> check -> apply. */
  template <class Output, class Compute, class Apply>
  [[nodiscard]] SLIDE_AGEING_FORCE_INLINE slide::Status evaluate_ageing_stage(
    bool enabled,
    Output output,
    Compute &&compute,
    Apply &&apply)
  {
    if (!enabled)
      return slide::Status::Success;
    const auto status = compute(output);
    if (status != slide::Status::Success)
      return status;
    apply(output);
    return slide::Status::Success;
  }

} // namespace detail
} // namespace slide::core

#undef SLIDE_AGEING_FORCE_INLINE
