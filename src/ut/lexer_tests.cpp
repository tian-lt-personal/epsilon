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
  script::lexer lex{"+-*/%()."};
  EXPECT_TRUE(std::holds_alternative<script::token_op_plus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_minus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_mul>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_div>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_percent>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_lparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_rparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_dot>(lex().value()));
  EXPECT_TRUE(lex.drained());
}

TEST(lexer_tests, integer_literal) {
  constexpr auto verify = [](std::string_view expected_raw, std::expected<script::token, script::token_ec> result) {
    EXPECT_TRUE(result.has_value());
    EXPECT_TRUE(std::holds_alternative<script::token_integer_literal>(*result));
    EXPECT_EQ(expected_raw, std::get<script::token_integer_literal>(*result).raw);
  };
  script::lexer lex{"z'0 z'000 z'123\nz'0123\r\rz'10000 z'1 \n z'100"};
  verify("0", lex());
  verify("000", lex());
  verify("123", lex());
  verify("0123", lex());
  verify("10000", lex());
  verify("1", lex());
  verify("100", lex());
  EXPECT_TRUE(lex.drained());
}

TEST(lexer_tests, real_literal) {
  constexpr auto verify = [](std::string_view expected_raw, std::expected<script::token, script::token_ec> result) {
    EXPECT_TRUE(result.has_value());
    EXPECT_TRUE(std::holds_alternative<script::token_real_literal>(*result));
    EXPECT_EQ(expected_raw, std::get<script::token_real_literal>(*result).raw);
  };
  script::lexer lex{
      "r'0 r'000 r'123\nr'0123\r\rr'10000 r'1 \n r'100\n"
      "0 123 .123 123. 123.456"};
  verify("0", lex());
  verify("000", lex());
  verify("123", lex());
  verify("0123", lex());
  verify("10000", lex());
  verify("1", lex());
  verify("100", lex());
  verify("0", lex());
  verify("123", lex());
  verify(".123", lex());
  verify("123.", lex());
  verify("123.456", lex());
  EXPECT_TRUE(lex.drained());
}

}  // namespace epxut
