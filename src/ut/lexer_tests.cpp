// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "lexer.hpp"

namespace script = epx::script;

namespace epxut {

TEST(lexer_tests, empty_input) {
  {
    script::lexer lex{""};
    EXPECT_TRUE(lex.drained());
    auto token = lex();
    EXPECT_FALSE(token.has_value());
    EXPECT_EQ(token.error(), script::token_ec::eof);
  }
  {
    script::lexer lex{"        \n\r\t  "};
    EXPECT_FALSE(lex.drained());
    auto token = lex();
    EXPECT_FALSE(token.has_value());
    EXPECT_EQ(token.error(), script::token_ec::eof);
    EXPECT_EQ(token.error(), script::token_ec::eof);  // drained state can be queried for any times.
  }
}

TEST(lexer_tests, operator_tokens) {
  script::lexer lex{"+-*/%()"};
  EXPECT_TRUE(std::holds_alternative<script::token_op_plus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_minus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_mul>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_div>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_percent>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_lparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_rparen>(lex().value()));
  EXPECT_TRUE(lex.drained());
}

TEST(lexer_tests, val_decimal_tokens) {
  //constexpr auto verify_decimal = [](std::string_view expected_raw,
  //                                   std::expected<script::token, script::token_ec> result) {
  //  EXPECT_TRUE(result.has_value());
  //  EXPECT_TRUE(std::holds_alternative<script::token_val_decimal>(*result));
  //  EXPECT_EQ(expected_raw, std::get<script::token_val_decimal>(*result).raw);
  //};
  //script::lexer lex{"0 123 123.456 123. .123 +1 -2 +.123 -.321"};
  //verify_decimal("0", lex());
  //verify_decimal("123", lex());
  //verify_decimal("123.456", lex());
  //verify_decimal("123.", lex());
  //verify_decimal(".123", lex());
  //verify_decimal("+1", lex());
  //verify_decimal("-2", lex());
  //verify_decimal("+.123", lex());
  //verify_decimal("-.321", lex());
  //EXPECT_TRUE(lex.drained());
}

}  // namespace epxut
