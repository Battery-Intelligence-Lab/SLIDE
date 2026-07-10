/** Raw-byte libFuzzer driver for liionpack CSV absorption. */

#include "FuzzSupport.hpp"

#include <cstdint>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t *data,
                                      std::size_t size)
try {
  const char *bytes = size == 0 ? "" : reinterpret_cast<const char *>(data);
  const std::string_view source{ bytes, size };
  const auto poison = slide::fuzz::poisonTopology();
  auto first = poison;
  auto second = poison;
  slide::core::NetlistCsvDiagnostic first_diagnostic, second_diagnostic;
  const auto first_status = slide::core::parseLiionpackNetlistCsv(
    source, first, first_diagnostic);
  const auto second_status = slide::core::parseLiionpackNetlistCsv(
    source, second, second_diagnostic);
  slide::fuzz::require(first_status == second_status);
  slide::fuzz::require(first_diagnostic.row == second_diagnostic.row);
  slide::fuzz::require(first_diagnostic.offset == second_diagnostic.offset);
  slide::fuzz::require(first_diagnostic.message == second_diagnostic.message);
  slide::fuzz::require(slide::fuzz::sameTopology(first, second));
  if (first_status == slide::Status::Success)
    slide::fuzz::require(slide::fuzz::validImportedTopology(first));
  else
    slide::fuzz::require(slide::fuzz::sameTopology(first, poison));
  return 0;
} catch (...) {
  slide::fuzz::trap();
}
