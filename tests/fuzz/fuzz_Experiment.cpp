/** Raw-byte libFuzzer driver for the Experiment grammar. */

#include "FuzzSupport.hpp"

#include <cstdint>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t *data,
                                      std::size_t size)
try {
  std::vector<std::string> steps;
  if (size != 0) {
    const std::string_view source{ reinterpret_cast<const char *>(data), size };
    std::size_t begin{};
    while (true) {
      const auto end = source.find('\n', begin);
      steps.emplace_back(source.substr(
        begin, end == std::string_view::npos ? source.size() - begin : end - begin));
      if (end == std::string_view::npos)
        break;
      begin = end + 1;
      if (begin == source.size())
        break;
    }
  }

  const auto poison = slide::fuzz::poisonExperiment();
  auto first = poison;
  auto second = poison;
  slide::core::ParseDiagnostic first_diagnostic, second_diagnostic;
  const auto first_status = slide::core::Experiment::parse(
    steps, first, first_diagnostic);
  const auto second_status = slide::core::Experiment::parse(
    steps, second, second_diagnostic);
  slide::fuzz::require(first_status == second_status);
  slide::fuzz::require(first_diagnostic.step == second_diagnostic.step);
  slide::fuzz::require(first_diagnostic.offset == second_diagnostic.offset);
  slide::fuzz::require(first_diagnostic.message == second_diagnostic.message);
  slide::fuzz::require(slide::fuzz::sameExperiment(first, second));
  if (first_status == slide::Status::Success)
    slide::fuzz::require(slide::fuzz::validExperiment(first));
  else
    slide::fuzz::require(slide::fuzz::sameExperiment(first, poison));
  return 0;
} catch (...) {
  slide::fuzz::trap();
}
