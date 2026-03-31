// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "lexer.hpp"

namespace epxut {

TEST(lexer_tests, empty_input) {
  {
    epx::lexer lex{""};
    EXPECT_TRUE(lex.drained());
    auto token = lex();
    EXPECT_FALSE(token.has_value());
    EXPECT_EQ(token.error().code, epx::token_ec::eof);
  }
  {
    epx::lexer lex{"        \n\r\t  "};
    EXPECT_FALSE(lex.drained());
    auto token = lex();
    EXPECT_FALSE(token.has_value());
    EXPECT_EQ(token.error().code, epx::token_ec::eof);
    EXPECT_EQ(token.error().code, epx::token_ec::eof);  // drained state can be queried for any times.
  }
}

TEST(lexer_tests, operator_tokens) {
  epx::lexer lex{"+-*/%()"};
  EXPECT_TRUE(std::holds_alternative<epx::token_op_add>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_op_sub>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_op_mul>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_op_div>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_op_percent>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_lparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<epx::token_rparen>(lex().value()));
  EXPECT_TRUE(lex.drained());
}
TEST(lexer_tests, val_decimal_tokens) {
  constexpr auto verify_decimal = [](std::string_view expected_raw,
                                     std::expected<epx::token, epx::token_error> result) {
    EXPECT_TRUE(result.has_value());
    EXPECT_TRUE(std::holds_alternative<epx::token_val_decimal>(*result));
    EXPECT_EQ(expected_raw, std::get<epx::token_val_decimal>(*result).raw);
  };
  epx::lexer lex{"0 123 123.456 123. .123 +1 -2 +.123 -.321"};
  verify_decimal("0", lex());
  verify_decimal("123", lex());
  verify_decimal("123.456", lex());
  verify_decimal("123.", lex());
  verify_decimal(".123", lex());
  verify_decimal("+1", lex());
  verify_decimal("-2", lex());
  verify_decimal("+.123", lex());
  verify_decimal("-.321", lex());
  EXPECT_TRUE(lex.drained());
}

}  // namespace epxut
