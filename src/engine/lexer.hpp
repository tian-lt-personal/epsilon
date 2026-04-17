// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_ENGINE_INC_LEXER_HPP
#define EPSILON_ENGINE_INC_LEXER_HPP

#include <expected>
#include <string_view>
#include <variant>

namespace epx::script {

struct token_integer_literal {
  std::string_view raw;
};
struct token_real_literal {
  std::string_view raw;
};
struct token_id {
  std::string_view raw;
};
struct token_op_eq {};
struct token_op_dot {};
struct token_op_comma {};
struct token_op_plus {};
struct token_op_minus {};
struct token_op_mul {};
struct token_op_div {};
struct token_op_percent {};
struct token_lparen {};
struct token_rparen {};

using token = std::variant<token_integer_literal, token_real_literal, token_id, token_op_eq, token_op_dot,
                           token_op_comma, token_op_plus, token_op_minus, token_op_mul, token_op_div, token_op_percent,
                           token_lparen, token_rparen>;

enum struct token_ec { eof, bad_input };

class lexer {
 private:
  std::string_view input_;
  const char* cursor_;

 public:
  constexpr explicit lexer(std::string_view input) noexcept : input_(input), cursor_(input.data()) {}
  constexpr bool drained() noexcept { return cursor_ >= input_.data() + input_.length(); }
  std::expected<token, token_ec> operator()() noexcept;
};

}  // namespace epx::script

#endif  // EPSILON_ENGINE_INC_LEXER_HPP
