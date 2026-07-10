/** Raw-byte libFuzzer driver for BPX JSON absorption. */

#include "FuzzSupport.hpp"

#include <cstdint>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t *data,
                                      std::size_t size)
try {
  const char *bytes = size == 0 ? "" : reinterpret_cast<const char *>(data);
  const std::string_view source{ bytes, size };
  const auto poison = slide::fuzz::poisonParameters();
  auto first = poison;
  auto second = poison;
  std::string first_diagnostic, second_diagnostic;
  const auto first_status = slide::core::ParameterSet::fromBpxJson(
    source, first, first_diagnostic);
  const auto second_status = slide::core::ParameterSet::fromBpxJson(
    source, second, second_diagnostic);
  slide::fuzz::require(first_status == second_status);
  slide::fuzz::require(first_diagnostic == second_diagnostic);
  slide::fuzz::require(slide::fuzz::sameParameters(first, second));
  if (first_status == slide::Status::Success) {
    slide::core::SpmFactoryInput input;
    slide::fuzz::require(!slide::fuzz::hasPoisonParameters(first));
    slide::fuzz::require(first.toSpmInput(input) == slide::Status::Success);
  } else {
    slide::fuzz::require(slide::fuzz::sameParameters(first, poison));
  }
  return 0;
} catch (...) {
  slide::fuzz::trap();
}
