/*
 * Histogram.hpp
 *
 * Histogram code to hold histogram data.
 *
 *  Created on: 14 Dec 2021
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <span>
#include <fstream>

namespace slide {

enum class HistogramType {
  slidepack = 0,
  equidistant = 1

};

template <HistogramType histogramType = HistogramType::equidistant, typename Tdata = double>
class Histogram
{
  std::vector<size_t> bins;
  Tdata x_min{}, x_max{}, dx{ 1 };
  int Nbins;

public:
  Histogram() = default;

  Histogram(Tdata x_min_, Tdata x_max_, int Nbins_ = settings::DATASTORE_NHIST)
    : bins(Nbins_ + 2), x_min(x_min_), x_max(x_max_), dx((x_max_ - x_min_) / Nbins_), Nbins(Nbins_)
  {
  }

  constexpr void add(Tdata x) noexcept
  {
    /*
     * increase the correct bin counter
     * data is in bin i if edge[i-1] <= data < edge[i]
     * bins[0] is number of elements less than x_min
     * bins[end] is number of elements more than x_max
     * */
    auto i = static_cast<int>(1 + (x - x_min) / dx);

    const auto ssize = static_cast<int>(bins.size());
    i = std::max(0, std::min(i, ssize - 1)); // #TODO what happens if container is empty?
    bins[i]++;
  }

  std::span<const size_t> viewBinValues() const noexcept { return bins; }

  [[nodiscard]] constexpr auto begin() noexcept { return bins.begin(); }
  [[nodiscard]] constexpr auto end() noexcept { return bins.end(); }

  auto getEdgeValues() const noexcept { return range_fix(x_min, x_max, dx); }

  auto size() const noexcept { return bins.size(); }

  friend std::ostream &operator<<(std::ostream &ofs, const Histogram<> &hist);
};

inline std::ostream &operator<<(std::ostream &ofs, const Histogram<> &hist)
{
  auto bins = hist.viewBinValues();
  auto edges = hist.getEdgeValues();

  ofs << "Edges:\n";

  for (auto edge : edges)
    ofs << edge << ',';

  ofs << '\n';
  ofs << "Bins:\n";

  for (auto bin : bins)
    ofs << bin << ',';

  ofs << '\n';

  return ofs;
}

static Histogram<> EmptyHistogram{};

} // namespace slide
