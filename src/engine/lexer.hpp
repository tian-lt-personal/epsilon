// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_ENGINE_INC_LEXER_HPP
#define EPSILON_ENGINE_INC_LEXER_HPP

#include <expected>
#include <string_view>
#include <variant>

namespace epx {

struct token_val_decimal {
  std::string_view raw;
};
struct token_op_add {};
struct token_op_sub {};
struct token_op_mul {};
struct token_op_div {};
struct token_op_percent {};
struct token_lparen {};
struct token_rparen {};

using token = std::variant<token_val_decimal, token_op_add, token_op_sub, token_op_mul, token_op_div, token_op_percent,
                           token_lparen, token_rparen>;

enum struct token_ec { eof, bad_input };
struct token_error {
  token_ec code;
  std::size_t line;
  std::size_t column;
};

class lexer {
 private:
  std::string_view input_;
  const char* cursor_;
  std::size_t line_ = 1;
  std::size_t column_ = 1;

 public:
  constexpr explicit lexer(std::string_view input) noexcept : input_(input), cursor_(input.data()) {}
  constexpr bool drained() noexcept { return cursor_ >= input_.data() + input_.length(); }
  std::expected<token, token_error> operator()() noexcept;
};

}  // namespace epx

#endif  // EPSILON_ENGINE_INC_LEXER_HPP
