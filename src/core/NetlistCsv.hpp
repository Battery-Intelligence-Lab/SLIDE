/**
 * @file NetlistCsv.hpp
 * @brief Bounded liionpack-compatible CSV absorption into a compiled topology.
 * @surface api
 */

#pragma once

#include "PackTopology.hpp"

#include <filesystem>
#include <string>
#include <string_view>

namespace slide::core {

struct NetlistCsvOptions
{
  std::string archetype{ "cell" };
  bool thermal{};
};

struct NetlistCsvDiagnostic
{
  std::size_t row{};    //!< one-based logical CSV row, zero before row parsing
  std::size_t offset{}; //!< source cursor at/immediately after a reader-detected
                        //!< syntax/limit failure; zero when unavailable
  std::string message{};
};

/**
 * Parse liionpack's `desc,node1,node2,value` table (additional coordinate
 * columns are ignored). `V*` rows become cells, every `R*` row including
 * `Ri*` becomes a literal fixed link, and the single `I*` row identifies pack
 * terminals. V/I magnitudes are validated but runtime models/experiments own
 * the actual cell voltage and applied current. This is graph-schema
 * compatibility, not behavioural parity or a lossless round-trip format.
 *
 * The operation is atomic: `output` is unchanged on every non-Success result.
 * Input is bounded to 4 MiB, 100,000 rows, 32 columns, and 65,536 bytes/field.
 */
[[nodiscard]] slide::Status parseLiionpackNetlistCsv(
  std::string_view csv,
  CompiledPackTopology &output,
  NetlistCsvDiagnostic &diagnostic,
  const NetlistCsvOptions &options = {});

/** Bounded exact-read file wrapper around parseLiionpackNetlistCsv(). */
[[nodiscard]] slide::Status loadLiionpackNetlistCsv(
  const std::filesystem::path &path,
  CompiledPackTopology &output,
  NetlistCsvDiagnostic &diagnostic,
  const NetlistCsvOptions &options = {});

} // namespace slide::core
