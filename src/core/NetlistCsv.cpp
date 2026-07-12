/**
 * @file NetlistCsv.cpp
 * @brief Bounded liionpack-compatible CSV absorption.
 */

#include "NetlistCsv.hpp"
#include "BoundedFileReader.hpp"
#include "PackTopologyInternal.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <fstream>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <stdexcept>
#include <vector>

namespace slide::core {
namespace {

  constexpr std::size_t max_csv_bytes = 4U * 1024U * 1024U;
  constexpr std::size_t max_csv_rows = 100'000;
  constexpr std::size_t max_csv_columns = 32;
  constexpr std::size_t max_csv_field_bytes = 65'536;
  constexpr std::size_t max_descriptor_bytes = 127;

  void assignDiagnosticNoThrow(std::string &target,
                               std::string_view message) noexcept
  {
    try {
      target.assign(message);
    } catch (...) {
      target.clear();
    }
  }

  slide::Status allocationFailure(NetlistCsvDiagnostic &diagnostic,
                                  std::string_view message) noexcept
  {
    diagnostic.row = 0;
    diagnostic.offset = 0;
    assignDiagnosticNoThrow(diagnostic.message, message);
    return slide::Status::Numerical_failure;
  }

  class CsvReader
  {
  public:
    CsvReader(std::string_view source, NetlistCsvDiagnostic &diagnostic)
      : source_{ source }, diagnostic_{ diagnostic }
    {}

    bool next(std::vector<std::string> &fields, bool &has_row)
    {
      has_row = false;
      fields.clear();
      if (cursor_ == source_.size())
        return true;
      while (true) {
        if (fields.size() >= max_csv_columns)
          return fail("CSV row exceeds 32 columns");
        std::string field;
        if (source_[cursor_] == '"') {
          ++cursor_;
          bool closed{};
          while (cursor_ < source_.size()) {
            const char value = source_[cursor_++];
            if (value == '"') {
              if (cursor_ < source_.size() && source_[cursor_] == '"') {
                ++cursor_;
                field.push_back('"');
              } else {
                closed = true;
                break;
              }
            } else {
              field.push_back(value);
            }
            if (field.size() > max_csv_field_bytes)
              return fail("CSV field exceeds 65536 bytes");
          }
          if (!closed)
            return fail("unterminated quoted CSV field");
          if (cursor_ < source_.size() && source_[cursor_] != ','
              && source_[cursor_] != '\r' && source_[cursor_] != '\n')
            return fail("unexpected byte after quoted CSV field");
        } else {
          while (cursor_ < source_.size() && source_[cursor_] != ','
                 && source_[cursor_] != '\r' && source_[cursor_] != '\n') {
            if (source_[cursor_] == '"')
              return fail("quote inside unquoted CSV field");
            field.push_back(source_[cursor_++]);
            if (field.size() > max_csv_field_bytes)
              return fail("CSV field exceeds 65536 bytes");
          }
        }
        fields.push_back(std::move(field));
        if (cursor_ == source_.size()) {
          has_row = true;
          ++row_;
          return true;
        }
        if (source_[cursor_] == ',') {
          ++cursor_;
          if (cursor_ == source_.size()) {
            fields.emplace_back();
            has_row = true;
            ++row_;
            return true;
          }
          continue;
        }
        if (source_[cursor_] == '\r') {
          if (cursor_ + 1 >= source_.size() || source_[cursor_ + 1] != '\n')
            return fail("bare carriage return in CSV");
          cursor_ += 2;
        } else {
          ++cursor_;
        }
        has_row = true;
        ++row_;
        return true;
      }
    }

    std::size_t row() const noexcept { return row_; }

  private:
    bool fail(std::string_view message)
    {
      diagnostic_.row = row_;
      diagnostic_.offset = cursor_;
      diagnostic_.message.assign(message);
      return false;
    }

    std::string_view source_{};
    NetlistCsvDiagnostic &diagnostic_;
    std::size_t cursor_{};
    std::size_t row_{ 1 };
  };

  bool parseNode(std::string_view text, std::uint32_t &node)
  {
    if (text.empty())
      return false;
    const auto parsed = std::from_chars(
      text.data(), text.data() + text.size(), node, 10);
    return parsed.ec == std::errc{}
           && parsed.ptr == text.data() + text.size();
  }

  bool parseValue(std::string_view text, real_t &value)
  {
    if (text.empty())
      return false;
    std::size_t cursor{};
    if (text[cursor] == '-' && ++cursor == text.size())
      return false;
    if (text[cursor] == '0') {
      ++cursor;
      if (cursor < text.size() && text[cursor] >= '0' && text[cursor] <= '9')
        return false;
    } else if (text[cursor] >= '1' && text[cursor] <= '9') {
      while (cursor < text.size() && text[cursor] >= '0'
             && text[cursor] <= '9')
        ++cursor;
    } else {
      return false;
    }
    if (cursor < text.size() && text[cursor] == '.') {
      const std::size_t fraction = ++cursor;
      while (cursor < text.size() && text[cursor] >= '0'
             && text[cursor] <= '9')
        ++cursor;
      if (cursor == fraction)
        return false;
    }
    if (cursor < text.size()
        && (text[cursor] == 'e' || text[cursor] == 'E')) {
      ++cursor;
      if (cursor < text.size()
          && (text[cursor] == '+' || text[cursor] == '-'))
        ++cursor;
      const std::size_t exponent = cursor;
      while (cursor < text.size() && text[cursor] >= '0'
             && text[cursor] <= '9')
        ++cursor;
      if (cursor == exponent)
        return false;
    }
    if (cursor != text.size())
      return false;
    const auto parsed = std::from_chars(
      text.data(), text.data() + text.size(), value, std::chars_format::general);
    return parsed.ec == std::errc{}
           && parsed.ptr == text.data() + text.size() && is_finite(value);
  }

  bool validDescriptor(std::string_view descriptor)
  {
    if (descriptor.empty() || descriptor.size() > max_descriptor_bytes)
      return false;
    for (const unsigned char value : descriptor)
      if (!((value >= 'a' && value <= 'z')
            || (value >= 'A' && value <= 'Z')
            || (value >= '0' && value <= '9') || value == '_'
            || value == '-' || value == '.'))
        return false;
    return true;
  }

  bool failSemantic(NetlistCsvDiagnostic &diagnostic, std::size_t row,
                    std::string_view message)
  {
    diagnostic.row = row;
    diagnostic.message.assign(message);
    return false;
  }

  struct ParsedElement
  {
    std::string descriptor{};
    std::uint32_t node_positive{};
    std::uint32_t node_negative{};
    real_t value{};
    char kind{};
    std::size_t row{};
  };

  class DisjointSet
  {
  public:
    explicit DisjointSet(std::size_t size) : parent_(size)
    {
      for (std::size_t i = 0; i < size; ++i)
        parent_[i] = i;
    }

    std::size_t find(std::size_t node)
    {
      std::size_t root = node;
      while (parent_[root] != root)
        root = parent_[root];
      while (parent_[node] != node) {
        const auto next = parent_[node];
        parent_[node] = root;
        node = next;
      }
      return root;
    }

    void unite(std::size_t left, std::size_t right)
    {
      left = find(left);
      right = find(right);
      if (left == right)
        return;
      if (right < left)
        std::swap(left, right);
      parent_[right] = left;
    }

  private:
    std::vector<std::size_t> parent_{};
  };

} // namespace

slide::Status parseLiionpackNetlistCsv(
  std::string_view csv,
  CompiledPackTopology &output,
  NetlistCsvDiagnostic &diagnostic,
  const NetlistCsvOptions &options)
try {
  diagnostic = {};
  if (csv.empty()) {
    diagnostic.message = "liionpack CSV is empty";
    return slide::Status::Invalid_parameters;
  }
  if (csv.size() > max_csv_bytes) {
    diagnostic.message = "liionpack CSV exceeds 4194304 bytes";
    return slide::Status::Invalid_parameters;
  }
  if (csv.find('\0') != std::string_view::npos) {
    diagnostic.message = "NUL byte in liionpack CSV";
    return slide::Status::Invalid_parameters;
  }
  if (options.archetype.empty() || options.archetype.size() > 1024
      || options.archetype.find('\0') != std::string::npos) {
    diagnostic.message = "invalid imported cell archetype";
    return slide::Status::Invalid_parameters;
  }

  CsvReader reader{ csv, diagnostic };
  std::vector<std::string> fields;
  bool has_row{};
  if (!reader.next(fields, has_row) || !has_row) {
    if (diagnostic.message.empty())
      diagnostic.message = "liionpack CSV needs a header";
    return slide::Status::Invalid_parameters;
  }
  if (!fields.empty() && fields[0].size() >= 3
      && static_cast<unsigned char>(fields[0][0]) == 0xefU
      && static_cast<unsigned char>(fields[0][1]) == 0xbbU
      && static_cast<unsigned char>(fields[0][2]) == 0xbfU)
    fields[0].erase(0, 3);

  std::map<std::string, std::size_t, std::less<>> columns;
  for (std::size_t column = 0; column < fields.size(); ++column)
    if (!columns.emplace(fields[column], column).second) {
      diagnostic.row = 1;
      diagnostic.message = "duplicate liionpack CSV header";
      return slide::Status::Invalid_parameters;
    }
  constexpr std::array required_names{
    std::string_view{ "desc" }, std::string_view{ "node1" }, std::string_view{ "node2" }, std::string_view{ "value" }
  };
  std::array<std::size_t, required_names.size()> required{};
  for (std::size_t i = 0; i < required_names.size(); ++i) {
    const auto found = columns.find(required_names[i]);
    if (found == columns.end()) {
      diagnostic.row = 1;
      diagnostic.message = "liionpack CSV needs desc,node1,node2,value columns";
      return slide::Status::Invalid_parameters;
    }
    required[i] = found->second;
  }
  const std::size_t column_count = fields.size();

  std::vector<ParsedElement> elements;
  std::set<std::string, std::less<>> descriptors;
  std::vector<std::uint32_t> labels;
  std::size_t data_rows{};
  std::size_t cell_count{};
  std::size_t current_count{};
  while (true) {
    if (!reader.next(fields, has_row))
      return slide::Status::Invalid_parameters;
    if (!has_row)
      break;
    ++data_rows;
    if (data_rows > max_csv_rows) {
      diagnostic.row = reader.row() - 1;
      diagnostic.message = "liionpack CSV exceeds 100000 rows";
      return slide::Status::Invalid_parameters;
    }
    const std::size_t row = reader.row() - 1;
    if (fields.size() != column_count) {
      failSemantic(diagnostic, row, "liionpack CSV row width differs from header");
      return slide::Status::Invalid_parameters;
    }
    const auto &descriptor = fields[required[0]];
    ParsedElement element{ .descriptor = descriptor, .row = row };
    if (!validDescriptor(descriptor)
        || !(descriptor[0] == 'V' || descriptor[0] == 'R'
             || descriptor[0] == 'I')) {
      failSemantic(diagnostic, row, "unsupported liionpack descriptor");
      return slide::Status::Invalid_parameters;
    }
    if (!descriptors.emplace(descriptor).second) {
      failSemantic(diagnostic, row, "duplicate liionpack descriptor");
      return slide::Status::Invalid_parameters;
    }
    if (!parseNode(fields[required[1]], element.node_positive)
        || !parseNode(fields[required[2]], element.node_negative)) {
      failSemantic(diagnostic, row, "invalid liionpack node label");
      return slide::Status::Invalid_parameters;
    }
    if (element.node_positive == element.node_negative) {
      failSemantic(diagnostic, row, "liionpack element has identical endpoints");
      return slide::Status::Invalid_parameters;
    }
    if (!parseValue(fields[required[3]], element.value)) {
      failSemantic(diagnostic, row, "invalid liionpack element value");
      return slide::Status::Invalid_parameters;
    }
    element.kind = descriptor[0];
    if (element.kind == 'R' && element.value < 0.0) {
      failSemantic(diagnostic, row, "liionpack resistance must be non-negative");
      return slide::Status::Invalid_parameters;
    }
    cell_count += element.kind == 'V';
    current_count += element.kind == 'I';
    labels.push_back(element.node_positive);
    labels.push_back(element.node_negative);
    elements.push_back(std::move(element));
  }
  if (data_rows == 0 || cell_count == 0 || current_count != 1) {
    diagnostic.row = data_rows + 1;
    diagnostic.message = "liionpack CSV needs cells and exactly one current source";
    return slide::Status::Invalid_parameters;
  }
  if (cell_count > max_csv_bytes / options.archetype.size()) {
    diagnostic.message = "imported cell archetype exceeds retained-text budget";
    return slide::Status::Invalid_parameters;
  }

  std::sort(labels.begin(), labels.end());
  labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
  DisjointSet sets{ labels.size() };
  const auto labelIndex = [&](std::uint32_t label) {
    return static_cast<std::size_t>(
      std::lower_bound(labels.begin(), labels.end(), label) - labels.begin());
  };
  for (const auto &element : elements)
    if (element.kind == 'R' && element.value == 0.0)
      sets.unite(labelIndex(element.node_positive),
                 labelIndex(element.node_negative));

  std::vector<std::size_t> roots;
  roots.reserve(labels.size());
  for (std::size_t node = 0; node < labels.size(); ++node)
    roots.push_back(sets.find(node));
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
  if (roots.size() > std::numeric_limits<std::uint32_t>::max()) {
    diagnostic.message = "liionpack CSV has too many nodes";
    return slide::Status::Invalid_parameters;
  }
  std::vector<std::uint32_t> root_to_dense(
    labels.size(), std::numeric_limits<std::uint32_t>::max());
  for (std::size_t dense = 0; dense < roots.size(); ++dense)
    root_to_dense[roots[dense]] = static_cast<std::uint32_t>(dense);
  const auto denseNode = [&](std::uint32_t label) {
    return root_to_dense[sets.find(labelIndex(label))];
  };

  CompiledPackTopology candidate;
  candidate.cells.reserve(cell_count);
  candidate.electrical.branches.reserve(elements.size() - current_count);
  for (const auto &element : elements) {
    const auto positive = denseNode(element.node_positive);
    const auto negative = denseNode(element.node_negative);
    if (element.kind == 'R' && element.value == 0.0)
      continue;
    if (positive == negative) {
      failSemantic(diagnostic, element.row, "liionpack ideal wire shorts an element or terminal source");
      return slide::Status::Invalid_parameters;
    }
    if (element.kind == 'V') {
      const auto cell = static_cast<std::uint32_t>(candidate.cells.size());
      candidate.cells.push_back({ .path = element.descriptor,
                                  .archetype = options.archetype,
                                  .thermal = options.thermal });
      candidate.electrical.branches.push_back(
        { .node_positive = positive,
          .node_negative = negative,
          .kind = ElectricalBranchKind::cell,
          .cell = cell });
    } else if (element.kind == 'R') {
      candidate.electrical.branches.push_back(
        { .node_positive = positive,
          .node_negative = negative,
          .kind = ElectricalBranchKind::resistor,
          .resistance = element.value });
    } else {
      candidate.electrical.terminal_positive = positive;
      candidate.electrical.terminal_negative = negative;
    }
  }

  const auto status = detail::finalizeImportedPackTopology(
    candidate, static_cast<std::uint32_t>(roots.size()));
  if (status != slide::Status::Success) {
    diagnostic.message = "liionpack CSV graph is disconnected or inconsistent";
    return status;
  }
  diagnostic = {};
  output = std::move(candidate);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return allocationFailure(diagnostic, "liionpack CSV allocation failed");
} catch (const std::length_error &) {
  return allocationFailure(
    diagnostic, "liionpack CSV size is not representable");
}

slide::Status loadLiionpackNetlistCsv(
  const std::filesystem::path &path,
  CompiledPackTopology &output,
  NetlistCsvDiagnostic &diagnostic,
  const NetlistCsvOptions &options)
try {
  std::string contents;
  const auto read = detail::readBoundedFile(path, max_csv_bytes, contents);
  if (read != detail::BoundedFileRead::success) {
    if (read == detail::BoundedFileRead::open_failed)
      diagnostic = { .message = "could not open liionpack CSV file" };
    else if (read == detail::BoundedFileRead::too_large)
      diagnostic = { .message = "liionpack CSV exceeds 4194304 bytes" };
    else
      diagnostic = { .message = "could not read complete liionpack CSV file" };
    return detail::boundedFileStatus(read);
  }
  return parseLiionpackNetlistCsv(contents, output, diagnostic, options);
} catch (const std::bad_alloc &) {
  return allocationFailure(
    diagnostic, "liionpack CSV file allocation failed");
} catch (const std::length_error &) {
  return allocationFailure(
    diagnostic, "liionpack CSV file size is not representable");
}

} // namespace slide::core
