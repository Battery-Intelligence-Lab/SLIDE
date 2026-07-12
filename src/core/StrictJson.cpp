/**
 * @file StrictJson.cpp
 * @brief Resource-bounded strict JSON grammar and UTF-8 scanner.
 */

#include "detail/StrictJson.hpp"

#include "Numeric.hpp"

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace slide::core::detail {
namespace {

  constexpr std::size_t max_json_values = 65'536;
  using JsonValue = StrictJsonValue;

  class JsonParser
  {
  public:
    JsonParser(std::string_view source, std::string &diagnostic)
      : source_{ source }, diagnostic_{ diagnostic }
    {}

    bool parse(JsonValue &output)
    {
      diagnostic_.clear();
      skipSpace();
      if (!parseValue(output, 0))
        return false;
      skipSpace();
      if (cursor_ != source_.size())
        return fail("unexpected trailing JSON data");
      return true;
    }

  private:
    bool parseValue(JsonValue &output, int depth)
    {
      if (depth > 64)
        return fail("JSON nesting exceeds 64 levels");
      if (values_ >= max_json_values)
        return fail("JSON value count exceeds 65536");
      ++values_;
      skipSpace();
      if (cursor_ == source_.size())
        return fail("unexpected end of JSON");
      const char token = source_[cursor_];
      if (token == '{') return parseObject(output, depth + 1);
      if (token == '[') return parseArray(output, depth + 1);
      if (token == '"') {
        output.kind = JsonValue::Kind::string;
        return parseString(output.string);
      }
      if (token == '-' || (token >= '0' && token <= '9'))
        return parseNumber(output);
      if (consume("true")) {
        output.kind = JsonValue::Kind::boolean;
        output.boolean = true;
        return true;
      }
      if (consume("false")) {
        output.kind = JsonValue::Kind::boolean;
        output.boolean = false;
        return true;
      }
      if (consume("null")) {
        output.kind = JsonValue::Kind::null_value;
        return true;
      }
      return fail("invalid JSON value");
    }

    bool parseObject(JsonValue &output, int depth)
    {
      ++cursor_;
      output.kind = JsonValue::Kind::object;
      skipSpace();
      if (take('}')) return true;
      while (true) {
        std::string key;
        if (!parseString(key))
          return false;
        skipSpace();
        if (!take(':'))
          return fail("expected ':' after object key");
        JsonValue value;
        if (!parseValue(value, depth))
          return false;
        if (!output.object.emplace(std::move(key), std::move(value)).second)
          return fail("duplicate JSON object key");
        skipSpace();
        if (take('}')) return true;
        if (!take(','))
          return fail("expected ',' or '}' in object");
        skipSpace();
      }
    }

    bool parseArray(JsonValue &output, int depth)
    {
      ++cursor_;
      output.kind = JsonValue::Kind::array;
      skipSpace();
      if (take(']')) return true;
      while (true) {
        JsonValue value;
        if (!parseValue(value, depth))
          return false;
        output.array.push_back(std::move(value));
        skipSpace();
        if (take(']')) return true;
        if (!take(','))
          return fail("expected ',' or ']' in array");
        skipSpace();
      }
    }

    bool parseString(std::string &output)
    {
      skipSpace();
      if (!take('"'))
        return fail("expected JSON string");
      output.clear();
      while (cursor_ < source_.size()) {
        const unsigned char c = static_cast<unsigned char>(source_[cursor_++]);
        if (c == '"') return true;
        if (c < 0x20)
          return fail("control byte in JSON string");
        if (c != '\\') {
          if (c < 0x80U)
            output.push_back(static_cast<char>(c));
          else if (!appendRawUtf8(output, c))
            return false;
          continue;
        }
        if (cursor_ == source_.size())
          return fail("unfinished JSON escape");
        const char escape = source_[cursor_++];
        switch (escape) {
        case '"':
          output.push_back('"');
          break;
        case '\\':
          output.push_back('\\');
          break;
        case '/':
          output.push_back('/');
          break;
        case 'b':
          output.push_back('\b');
          break;
        case 'f':
          output.push_back('\f');
          break;
        case 'n':
          output.push_back('\n');
          break;
        case 'r':
          output.push_back('\r');
          break;
        case 't':
          output.push_back('\t');
          break;
        case 'u': {
          std::uint32_t codepoint{};
          if (!parseHex4(codepoint))
            return false;
          if (codepoint >= 0xd800U && codepoint <= 0xdbffU) {
            if (cursor_ + 2 > source_.size() || source_[cursor_] != '\\'
                || source_[cursor_ + 1] != 'u')
              return fail("high surrogate without low surrogate");
            cursor_ += 2;
            std::uint32_t low{};
            if (!parseHex4(low))
              return false;
            if (low < 0xdc00U || low > 0xdfffU)
              return fail("invalid low surrogate");
            codepoint = 0x10000U + ((codepoint - 0xd800U) << 10U)
                        + (low - 0xdc00U);
          } else if (codepoint >= 0xdc00U && codepoint <= 0xdfffU) {
            return fail("unpaired low surrogate");
          }
          appendUtf8(output, codepoint);
          break;
        }
        default:
          return fail("unsupported JSON escape (use UTF-8 directly)");
        }
      }
      return fail("unterminated JSON string");
    }

    bool parseNumber(JsonValue &output)
    {
      const std::size_t begin = cursor_;
      if (take('-') && cursor_ == source_.size())
        return fail("invalid JSON number");
      if (cursor_ == source_.size())
        return fail("invalid JSON number");
      if (source_[cursor_] == '0') {
        ++cursor_;
        if (cursor_ < source_.size() && source_[cursor_] >= '0'
            && source_[cursor_] <= '9')
          return fail("leading zero in JSON number");
      } else if (source_[cursor_] >= '1' && source_[cursor_] <= '9') {
        while (cursor_ < source_.size() && source_[cursor_] >= '0'
               && source_[cursor_] <= '9')
          ++cursor_;
      } else {
        return fail("invalid JSON number");
      }
      if (cursor_ < source_.size() && source_[cursor_] == '.') {
        ++cursor_;
        const std::size_t fraction = cursor_;
        while (cursor_ < source_.size() && source_[cursor_] >= '0'
               && source_[cursor_] <= '9')
          ++cursor_;
        if (cursor_ == fraction)
          return fail("JSON fraction needs a digit");
      }
      if (cursor_ < source_.size()
          && (source_[cursor_] == 'e' || source_[cursor_] == 'E')) {
        ++cursor_;
        if (cursor_ < source_.size()
            && (source_[cursor_] == '+' || source_[cursor_] == '-'))
          ++cursor_;
        const std::size_t exponent = cursor_;
        while (cursor_ < source_.size() && source_[cursor_] >= '0'
               && source_[cursor_] <= '9')
          ++cursor_;
        if (cursor_ == exponent)
          return fail("JSON exponent needs a digit");
      }

      const char *first = source_.data() + begin;
      const char *last = source_.data() + cursor_;
      real_t value{};
      const auto parsed =
        std::from_chars(first, last, value, std::chars_format::general);
      if (parsed.ec != std::errc{} || parsed.ptr != last || !is_finite(value))
        return fail("invalid JSON number");
      output.kind = JsonValue::Kind::number;
      output.number = value;
      return true;
    }

    bool appendRawUtf8(std::string &output, unsigned char lead)
    {
      const std::size_t begin = cursor_ - 1;
      std::size_t continuations{};
      std::uint32_t codepoint{};
      std::uint32_t minimum{};
      if (lead >= 0xc2U && lead <= 0xdfU) {
        continuations = 1;
        codepoint = lead & 0x1fU;
        minimum = 0x80U;
      } else if (lead >= 0xe0U && lead <= 0xefU) {
        continuations = 2;
        codepoint = lead & 0x0fU;
        minimum = 0x800U;
      } else if (lead >= 0xf0U && lead <= 0xf4U) {
        continuations = 3;
        codepoint = lead & 0x07U;
        minimum = 0x10000U;
      } else {
        return fail("invalid UTF-8 lead byte in JSON string");
      }
      if (continuations > source_.size() - cursor_)
        return fail("unfinished UTF-8 sequence in JSON string");
      for (std::size_t i = 0; i < continuations; ++i) {
        const auto continuation = static_cast<unsigned char>(source_[cursor_++]);
        if ((continuation & 0xc0U) != 0x80U)
          return fail("invalid UTF-8 continuation in JSON string");
        codepoint = (codepoint << 6U) | (continuation & 0x3fU);
      }
      if (codepoint < minimum || codepoint > 0x10ffffU
          || (codepoint >= 0xd800U && codepoint <= 0xdfffU))
        return fail("invalid UTF-8 scalar in JSON string");
      output.append(source_.substr(begin, continuations + 1));
      return true;
    }

    bool parseHex4(std::uint32_t &value)
    {
      if (cursor_ + 4 > source_.size())
        return fail("unfinished Unicode escape");
      value = 0;
      for (int i = 0; i < 4; ++i) {
        const char c = source_[cursor_++];
        unsigned digit{};
        if (c >= '0' && c <= '9')
          digit = static_cast<unsigned>(c - '0');
        else if (c >= 'a' && c <= 'f')
          digit = 10U + static_cast<unsigned>(c - 'a');
        else if (c >= 'A' && c <= 'F')
          digit = 10U + static_cast<unsigned>(c - 'A');
        else
          return fail("invalid hexadecimal Unicode escape");
        value = (value << 4U) | digit;
      }
      return true;
    }

    static void appendUtf8(std::string &output, std::uint32_t codepoint)
    {
      if (codepoint <= 0x7fU) {
        output.push_back(static_cast<char>(codepoint));
      } else if (codepoint <= 0x7ffU) {
        output.push_back(static_cast<char>(0xc0U | (codepoint >> 6U)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      } else if (codepoint <= 0xffffU) {
        output.push_back(static_cast<char>(0xe0U | (codepoint >> 12U)));
        output.push_back(
          static_cast<char>(0x80U | ((codepoint >> 6U) & 0x3fU)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      } else {
        output.push_back(static_cast<char>(0xf0U | (codepoint >> 18U)));
        output.push_back(
          static_cast<char>(0x80U | ((codepoint >> 12U) & 0x3fU)));
        output.push_back(
          static_cast<char>(0x80U | ((codepoint >> 6U) & 0x3fU)));
        output.push_back(static_cast<char>(0x80U | (codepoint & 0x3fU)));
      }
    }

    void skipSpace()
    {
      while (cursor_ < source_.size()
             && (source_[cursor_] == ' ' || source_[cursor_] == '\t'
                 || source_[cursor_] == '\r' || source_[cursor_] == '\n'))
        ++cursor_;
    }

    bool take(char expected)
    {
      if (cursor_ < source_.size() && source_[cursor_] == expected) {
        ++cursor_;
        return true;
      }
      return false;
    }

    bool consume(std::string_view text)
    {
      if (!source_.substr(cursor_).starts_with(text))
        return false;
      cursor_ += text.size();
      return true;
    }

    bool fail(std::string_view message)
    {
      diagnostic_ = std::string{ message } + " at byte "
                    + std::to_string(cursor_);
      return false;
    }

    std::string_view source_{};
    std::string &diagnostic_;
    std::size_t cursor_{};
    std::size_t values_{};
  };

} // namespace

bool parseStrictJson(std::string_view source,
                     StrictJsonValue &output,
                     std::string &diagnostic)
{
  static_assert(std::is_nothrow_move_assignable_v<StrictJsonValue>);
  StrictJsonValue candidate;
  JsonParser parser{ source, diagnostic };
  if (!parser.parse(candidate))
    return false;
  output = std::move(candidate);
  return true;
}

} // namespace slide::core::detail
