// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "chars.hpp"
#include "lexer.hpp"
#include "parser.hpp"
#include "utils.hpp"

namespace script = epx::script;

namespace {

template <class T>
concept is_valid_parser_translate = requires(T&& t) {
  { std::move(t).translate() };
} && !requires(T&& t) {
  { t.translate() };
};

}  // namespace

namespace epxut {

TEST(parser_tests, empty_input) {
  EXPECT_FALSE(script::translate("").has_value());
  EXPECT_FALSE(script::translate("   ").has_value());
  EXPECT_FALSE(script::translate("\n\r\n\n\n  \n\r\n ").has_value());
}
TEST(parser_tests, binary_ops) {
  {
    auto res = script::translate("1+2");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=+]\n"
        "          left=[val|real=1]\n"
        "          right=[val|real=2]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("z'321-123");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=-]\n"
        "          left=[val|integer=321]\n"
        "          right=[val|real=123]\n"
        "",
        script::dump(*res));
  }
  {
    auto res = script::translate("456*z'789");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=*]\n"
        "          left=[val|real=456]\n"
        "          right=[val|integer=789]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("321/r'987");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=/]\n"
        "          left=[val|real=321]\n"
        "          right=[val|real=987]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+2+3-4-5+6");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=+]\n"
        "          left=[bin_op_expr|op=-]\n"
        "                 left=[bin_op_expr|op=-]\n"
        "                        left=[bin_op_expr|op=+]\n"
        "                               left=[bin_op_expr|op=+]\n"
        "                                      left=[val|real=1]\n"
        "                                      right=[val|real=2]\n"
        "                               right=[val|real=3]\n"
        "                        right=[val|real=4]\n"
        "                 right=[val|real=5]\n"
        "          right=[val|real=6]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2/3/4/5*6");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=*]\n"
        "          left=[bin_op_expr|op=/]\n"
        "                 left=[bin_op_expr|op=/]\n"
        "                        left=[bin_op_expr|op=/]\n"
        "                               left=[bin_op_expr|op=*]\n"
        "                                      left=[val|real=1]\n"
        "                                      right=[val|real=2]\n"
        "                               right=[val|real=3]\n"
        "                        right=[val|real=4]\n"
        "                 right=[val|real=5]\n"
        "          right=[val|real=6]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+2*3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=+]\n"
        "          left=[val|real=1]\n"
        "          right=[bin_op_expr|op=*]\n"
        "                  left=[val|real=2]\n"
        "                  right=[val|real=3]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2-3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=-]\n"
        "          left=[bin_op_expr|op=*]\n"
        "                 left=[val|real=1]\n"
        "                 right=[val|real=2]\n"
        "          right=[val|real=3]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2-3/z'4/5");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=-]\n"
        "          left=[bin_op_expr|op=*]\n"
        "                 left=[val|real=1]\n"
        "                 right=[val|real=2]\n"
        "          right=[bin_op_expr|op=/]\n"
        "                  left=[bin_op_expr|op=/]\n"
        "                         left=[val|real=3]\n"
        "                         right=[val|integer=4]\n"
        "                  right=[val|real=5]\n",
        script::dump(*res));
  }
}

TEST(parser_tests, parenthesized_exprs) {
  {
    auto res = script::translate("(1)");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[paren_expr]\n"
        "          inner=[val|real=1]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("(1+2)*3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=*]\n"
        "          left=[paren_expr]\n"
        "                 inner=[bin_op_expr|op=+]\n"
        "                         left=[val|real=1]\n"
        "                         right=[val|real=2]\n"
        "          right=[val|real=3]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+(2/3)*(z'4-5)");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=+]\n"
        "          left=[val|real=1]\n"
        "          right=[bin_op_expr|op=*]\n"
        "                  left=[paren_expr]\n"
        "                         inner=[bin_op_expr|op=/]\n"
        "                                 left=[val|real=2]\n"
        "                                 right=[val|real=3]\n"
        "                  right=[paren_expr]\n"
        "                          inner=[bin_op_expr|op=-]\n"
        "                                  left=[val|integer=4]\n"
        "                                  right=[val|real=5]\n",
        script::dump(*res));
  }
}

TEST(parser_tests, value_terms) {
  {
    auto res = script::translate("1");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ("stmt[0]=[val|real=1]", script::dump(*res));
  }
  {
    auto res = script::translate("z'2");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ("stmt[0]=[val|integer=2]", script::dump(*res));
  }
  {
    auto res = script::translate("a0");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ("stmt[0]=[val|id=a0]", script::dump(*res));
  }
  {
    auto res = script::translate("a+b-z'3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=-]\n"
        "          left=[bin_op_expr|op=+]\n"
        "                 left=[val|id=a]\n"
        "                 right=[val|id=b]\n"
        "          right=[val|integer=3]\n",
        script::dump(*res));
  }
}

TEST(parser_tests, func_calls) {
  {
    auto res = script::translate("f()");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ("stmt[0]=[func_call|id=f, param_count=0]\n", script::dump(*res));
  }
  {
    auto res = script::translate("g(1)");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[func_call|id=g, param_count=1]\n"
        "          param[0]=[val|real=1]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("h(1,z'2,pi)");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[func_call|id=h, param_count=3]\n"
        "          param[0]=[val|real=1]\n"
        "          param[1]=[val|integer=2]\n"
        "          param[2]=[val|id=pi]\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+f(2, g(3, h(4)))*z(5)");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "stmt[0]=[bin_op_expr|op=+]\n"
        "          left=[val|real=1]\n"
        "          right=[bin_op_expr|op=*]\n"
        "                  left=[func_call|id=f, param_count=2]\n"
        "                         param[0]=[val|real=2]\n"
        "                         param[1]=[func_call|id=g, param_count=2]\n"
        "                                    param[0]=[val|real=3]\n"
        "                                    param[1]=[func_call|id=h, param_count=1]\n"
        "                                               param[0]=[val|real=4]\n"
        "                  right=[func_call|id=z, param_count=1]\n"
        "                          param[0]=[val|real=5]\n",
        script::dump(*res));
  }
}

TEST(parser_tests, parse_real_literal_integer) {
  auto res = script::parse_real_literal("42");
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(epx::to_string(res->num), "42");
  EXPECT_EQ(epx::to_string(res->den), "1");
}

TEST(parser_tests, parse_real_literal_decimal) {
  auto res = script::parse_real_literal("3.14");
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(epx::to_string(res->num), "314");
  EXPECT_EQ(epx::to_string(res->den), "100");
}

TEST(parser_tests, parse_real_literal_leading_dot) {
  auto res = script::parse_real_literal(".5");
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(epx::to_string(res->num), "5");
  EXPECT_EQ(epx::to_string(res->den), "10");
}

TEST(parser_tests, parse_real_literal_trailing_dot) {
  auto res = script::parse_real_literal("3.");
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(epx::to_string(res->num), "3");
  EXPECT_EQ(epx::to_string(res->den), "1");
}

TEST(parser_tests, parse_real_literal_zero) {
  auto res = script::parse_real_literal("0.0");
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(epx::to_string(res->num), "0");
  EXPECT_EQ(epx::to_string(res->den), "10");
}

TEST(parser_tests, parse_real_literal_bad_input) {
  EXPECT_FALSE(script::parse_real_literal("abc").has_value());
}

}  // namespace epxut
