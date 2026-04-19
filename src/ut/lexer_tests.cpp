// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "lexer.hpp"

namespace script = epx::script;
using lex_result = decltype(std::declval<epx::script::lexer>().operator()());

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
    auto token = lex();
    EXPECT_FALSE(token.has_value());
    EXPECT_EQ(token.error(), script::token_ec::eof);
    EXPECT_EQ(token.error(), script::token_ec::eof);  // drained state can be queried for any times.
  }
}

TEST(lexer_tests, drained) {
  EXPECT_TRUE(script::lexer{""}.drained());
  EXPECT_TRUE(script::lexer{" "}.drained());
  EXPECT_TRUE(script::lexer{"\n\r\t"}.drained());
  EXPECT_TRUE(script::lexer{"\n \r\t   "}.drained());
  EXPECT_FALSE(script::lexer{"abc"}.drained());
  EXPECT_FALSE(script::lexer{"   ++"}.drained());
}

TEST(lexer_tests, operators) {
  script::lexer lex{"+-*/%().=,"};
  EXPECT_TRUE(std::holds_alternative<script::token_op_plus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_minus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_mul>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_div>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_percent>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_lparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_rparen>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_dot>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_eq>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_comma>(lex().value()));
  EXPECT_TRUE(lex.drained());
}

TEST(lexer_tests, integer_literal) {
  constexpr auto verify = [](std::string_view expected_raw, lex_result result) {
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<script::token_integer_literal>(*result));
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
  constexpr auto verify = [](std::string_view expected_raw, lex_result result) {
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<script::token_real_literal>(*result));
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

TEST(lexer_tests, id) {
  constexpr auto verify = [](std::string_view expected_raw, lex_result result) {
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<script::token_id>(*result));
    EXPECT_EQ(expected_raw, std::get<script::token_id>(*result).raw);
  };
  script::lexer lex{"a r z a123 foo\n\r\fbar baz pi gamma"};
  verify("a", lex());
  verify("r", lex());
  verify("z", lex());
  verify("a123", lex());
  verify("foo", lex());
  verify("bar", lex());
  verify("baz", lex());
  verify("pi", lex());
  verify("gamma", lex());
  EXPECT_TRUE(lex.drained());
}

TEST(lexer_tests, id_corners) {
  constexpr auto verify = [](std::string_view expected_raw, lex_result result) {
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(std::holds_alternative<script::token_id>(*result));
    EXPECT_EQ(expected_raw, std::get<script::token_id>(*result).raw);
  };
  {
    script::lexer lex{"a"};
    verify("a", lex());
    EXPECT_TRUE(lex.drained());
  }
  {
    script::lexer lex{"r"};
    verify("r", lex());
    EXPECT_TRUE(lex.drained());
  }
  {
    script::lexer lex{"z"};
    verify("z", lex());
    EXPECT_TRUE(lex.drained());
  }
}

TEST(lexer_tests, bad_inputs) {
  constexpr auto verify = [](std::string_view bad_input) {
    script::lexer lex{bad_input};
    auto res = lex();
    if (res.has_value()) return false;
    if (script::token_ec::bad_input != res.error()) return false;
    return true;
  };
  EXPECT_TRUE(verify("1a"));       // bad id
  EXPECT_TRUE(verify("r'"));       // bad literal prefix
  EXPECT_TRUE(verify("z'"));       // bad literal prefix
  EXPECT_TRUE(verify("r'1a"));     // bad real number
  EXPECT_TRUE(verify("r'1.1."));   // bad real number
  EXPECT_TRUE(verify("r'.1.1"));   // bad real number
  EXPECT_TRUE(verify("r'.1.1."));  // bad real number
  EXPECT_TRUE(verify("r'."));      // bad real number
  EXPECT_TRUE(verify("r'.."));     // bad real number
  EXPECT_TRUE(verify("1.1."));     // bad real number
  EXPECT_TRUE(verify(".1.1."));    // bad real number
  EXPECT_TRUE(verify(".1.1"));     // bad real number
  EXPECT_TRUE(verify("z'1b"));     // bad integer
  EXPECT_TRUE(verify("z'1.1"));    // bad integer
  EXPECT_TRUE(verify("z'.1.1"));   // bad integer
  EXPECT_TRUE(verify("z'.1.1."));  // bad integer
  EXPECT_TRUE(verify("z'1.1."));   // bad integer
  EXPECT_TRUE(verify("z'."));      // bad integer
  EXPECT_TRUE(verify("z'.."));     // bad integer
}

TEST(lexer_tests, compound) {
  script::lexer lex{"z'2/1.2++r'2=z'42342 foo(2, pi, z'432)"};
  EXPECT_EQ("2", std::get<script::token_integer_literal>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_op_div>(lex().value()));
  EXPECT_EQ("1.2", std::get<script::token_real_literal>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_op_plus>(lex().value()));
  EXPECT_TRUE(std::holds_alternative<script::token_op_plus>(lex().value()));
  EXPECT_EQ("2", std::get<script::token_real_literal>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_op_eq>(lex().value()));
  EXPECT_EQ("42342", std::get<script::token_integer_literal>(lex().value()).raw);
  EXPECT_EQ("foo", std::get<script::token_id>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_lparen>(lex().value()));
  EXPECT_EQ("2", std::get<script::token_real_literal>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_op_comma>(lex().value()));
  EXPECT_EQ("pi", std::get<script::token_id>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_op_comma>(lex().value()));
  EXPECT_EQ("432", std::get<script::token_integer_literal>(lex().value()).raw);
  EXPECT_TRUE(std::holds_alternative<script::token_rparen>(lex().value()));
  EXPECT_TRUE(lex.drained());
}

}  // namespace epxut
