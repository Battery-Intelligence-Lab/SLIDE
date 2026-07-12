/**
 * @file BpxExpression.cpp
 * @brief Strict AST compilation and evaluation for cold BPX expressions.
 */

#include "detail/BpxExpression.hpp"

#include "Numeric.hpp"
#include "detail/ParameterCurve.hpp"

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace slide::core::detail {
namespace {

  /** Strict evaluator for BPX's documented one-variable expression language. */
  class BpxExpression
  {
  public:
    bool compile(std::string_view source, std::string &diagnostic)
    {
      source_ = source;
      cursor_ = 0;
      nodes_.clear();
      diagnostic.clear();
      if (source.empty() || source.size() > 65'536)
        return fail(diagnostic, "BPX expression is empty or too long");
      root_ = parseExpression(diagnostic, 0);
      skipSpace();
      if (root_ < 0 || cursor_ != source_.size()) {
        if (diagnostic.empty())
          fail(diagnostic, "unexpected BPX expression token");
        return false;
      }
      return true;
    }

    bool evaluate(real_t x, real_t &output) const
    {
      return root_ >= 0 && evaluateNode(root_, x, output, 0)
             && is_finite(output);
    }

  private:
    enum class Kind : unsigned char {
      literal,
      variable,
      add,
      subtract,
      multiply,
      divide,
      power,
      negate,
      exponential,
      hyperbolic_tangent,
      hyperbolic_cosine
    };

    struct Node
    {
      Kind kind{ Kind::literal };
      real_t value{};
      int left{ -1 };
      int right{ -1 };
    };

    int append(Node node, std::string &diagnostic)
    {
      if (nodes_.size() >= 1024) {
        fail(diagnostic, "BPX expression exceeds 1024 operations");
        return -1;
      }
      nodes_.push_back(node);
      return static_cast<int>(nodes_.size() - 1);
    }

    int parseExpression(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      int left = parseTerm(diagnostic, depth + 1);
      while (left >= 0) {
        skipSpace();
        Kind kind;
        if (take('+'))
          kind = Kind::add;
        else if (take('-'))
          kind = Kind::subtract;
        else
          break;
        const int right = parseTerm(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = kind, .left = left, .right = right },
                      diagnostic);
      }
      return left;
    }

    int parseTerm(std::string &diagnostic, int depth)
    {
      int left = parseUnary(diagnostic, depth + 1);
      while (left >= 0) {
        skipSpace();
        Kind kind;
        if (source_.substr(cursor_).starts_with("**"))
          break;
        if (take('*'))
          kind = Kind::multiply;
        else if (take('/'))
          kind = Kind::divide;
        else
          break;
        const int right = parseUnary(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = kind, .left = left, .right = right },
                      diagnostic);
      }
      return left;
    }

    int parseUnary(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      skipSpace();
      if (take('+'))
        return parseUnary(diagnostic, depth + 1);
      if (take('-')) {
        const int child = parseUnary(diagnostic, depth + 1);
        return child < 0 ? -1
                         : append({ .kind = Kind::negate, .left = child },
                                  diagnostic);
      }
      return parsePower(diagnostic, depth + 1);
    }

    int parsePower(std::string &diagnostic, int depth)
    {
      int left = parsePrimary(diagnostic, depth + 1);
      skipSpace();
      if (left >= 0 && consume("**")) {
        const int right = parseUnary(diagnostic, depth + 1);
        if (right < 0)
          return -1;
        left = append({ .kind = Kind::power, .left = left, .right = right },
                      diagnostic);
      }
      return left;
    }

    int parsePrimary(std::string &diagnostic, int depth)
    {
      if (depth > 128) {
        fail(diagnostic, "BPX expression nesting exceeds 128 levels");
        return -1;
      }
      skipSpace();
      if (take('(')) {
        const int result = parseExpression(diagnostic, depth + 1);
        skipSpace();
        if (result < 0 || !take(')')) {
          fail(diagnostic, "expected ')' in BPX expression");
          return -1;
        }
        return result;
      }
      if (cursor_ < source_.size()
          && (std::isdigit(static_cast<unsigned char>(source_[cursor_]))
              || source_[cursor_] == '.')) {
        const char *first = source_.data() + cursor_;
        const char *last = source_.data() + source_.size();
        real_t value{};
        const auto parsed =
          std::from_chars(first, last, value, std::chars_format::general);
        if (parsed.ec != std::errc{} || parsed.ptr == first
            || !is_finite(value)) {
          fail(diagnostic, "invalid number in BPX expression");
          return -1;
        }
        cursor_ = static_cast<std::size_t>(parsed.ptr - source_.data());
        return append({ .kind = Kind::literal, .value = value }, diagnostic);
      }
      if (cursor_ < source_.size()
          && std::isalpha(static_cast<unsigned char>(source_[cursor_]))) {
        const std::size_t begin = cursor_++;
        while (cursor_ < source_.size()
               && std::isalnum(static_cast<unsigned char>(source_[cursor_])))
          ++cursor_;
        const auto identifier = source_.substr(begin, cursor_ - begin);
        if (identifier == "x")
          return append({ .kind = Kind::variable }, diagnostic);
        Kind kind;
        if (identifier == "exp")
          kind = Kind::exponential;
        else if (identifier == "tanh")
          kind = Kind::hyperbolic_tangent;
        else if (identifier == "cosh")
          kind = Kind::hyperbolic_cosine;
        else {
          fail(diagnostic, "unsupported BPX expression identifier");
          return -1;
        }
        skipSpace();
        if (!take('(')) {
          fail(diagnostic, "expected '(' after BPX function");
          return -1;
        }
        const int child = parseExpression(diagnostic, depth + 1);
        skipSpace();
        if (child < 0 || !take(')')) {
          fail(diagnostic, "expected ')' after BPX function argument");
          return -1;
        }
        return append({ .kind = kind, .left = child }, diagnostic);
      }
      fail(diagnostic, "expected value in BPX expression");
      return -1;
    }

    bool evaluateNode(int index, real_t x, real_t &output, int depth) const
    {
      if (depth > 128 || index < 0
          || static_cast<std::size_t>(index) >= nodes_.size())
        return false;
      const auto &node = nodes_[static_cast<std::size_t>(index)];
      if (node.kind == Kind::literal) {
        output = node.value;
        return true;
      }
      if (node.kind == Kind::variable) {
        output = x;
        return is_finite(output);
      }
      real_t left{};
      if (!evaluateNode(node.left, x, left, depth + 1))
        return false;
      if (node.kind == Kind::negate)
        output = -left;
      else if (node.kind == Kind::exponential)
        output = std::exp(left);
      else if (node.kind == Kind::hyperbolic_tangent)
        output = std::tanh(left);
      else if (node.kind == Kind::hyperbolic_cosine)
        output = std::cosh(left);
      else {
        real_t right{};
        if (!evaluateNode(node.right, x, right, depth + 1))
          return false;
        switch (node.kind) {
        case Kind::add:
          output = left + right;
          break;
        case Kind::subtract:
          output = left - right;
          break;
        case Kind::multiply:
          output = left * right;
          break;
        case Kind::divide:
          if (right == 0.0)
            return false;
          output = left / right;
          break;
        case Kind::power:
          output = std::pow(left, right);
          break;
        default:
          return false;
        }
      }
      return is_finite(output);
    }

    void skipSpace()
    {
      while (cursor_ < source_.size()
             && std::isspace(static_cast<unsigned char>(source_[cursor_])))
        ++cursor_;
    }

    bool take(char token)
    {
      if (cursor_ < source_.size() && source_[cursor_] == token) {
        ++cursor_;
        return true;
      }
      return false;
    }

    bool consume(std::string_view token)
    {
      if (!source_.substr(cursor_).starts_with(token))
        return false;
      cursor_ += token.size();
      return true;
    }

    bool fail(std::string &diagnostic, std::string_view message)
    {
      if (diagnostic.empty())
        diagnostic = std::string{ message } + " at expression byte "
                     + std::to_string(cursor_);
      return false;
    }

    std::string_view source_{};
    std::size_t cursor_{};
    std::vector<Node> nodes_{};
    int root_{ -1 };
  };

} // namespace

bool evaluateBpxExpressionSamples(std::string_view source,
                                  std::span<const real_t>
                                    samples,
                                  std::span<real_t>
                                    output,
                                  std::string &diagnostic)
{
  if (samples.size() != output.size())
    return false;
  BpxExpression expression;
  if (!expression.compile(source, diagnostic))
    return false;
  for (std::size_t i = 0; i < samples.size(); ++i)
    if (!expression.evaluate(samples[i], output[i])) {
      diagnostic = "BPX scalar expression is non-finite";
      return false;
    }
  return true;
}

bool sampleBpxExpressionCurve(std::string_view source,
                              OCVCurve &output,
                              std::string &diagnostic)
{
  BpxExpression expression;
  if (!expression.compile(source, diagnostic))
    return false;
  auto candidate = sampleParameterCurve([&expression](real_t x) {
    real_t value_at_x{};
    return expression.evaluate(x, value_at_x)
             ? value_at_x
             : std::numeric_limits<real_t>::quiet_NaN();
  });
  if (!validParameterCurve(candidate)) {
    diagnostic = "BPX function is non-finite on stoichiometry [0,1]";
    return false;
  }
  output = std::move(candidate);
  return true;
}

} // namespace slide::core::detail
