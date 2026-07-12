// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "chars.hpp"
#include "machine.hpp"
#include "parser.hpp"

namespace script = epx::script;

namespace epxut {

namespace {

std::string eval_to_string(std::string_view expr, unsigned int precision = 5) {
  auto script_res = script::translate(expr);
  if (!script_res.has_value()) return "";
  script::machine m;
  auto exec_res = m.execute(*script_res);
  if (!exec_res.has_value() || exec_res->empty()) return "";
  auto& result = (*exec_res)[0];
  if (std::holds_alternative<script::real>(result)) {
    return epx::to_string(std::move(std::get<script::real>(result)), precision);
  }
  if (std::holds_alternative<script::integer>(result)) {
    return epx::to_string(std::get<script::integer>(result));
  }
  return "";
}

bool eval_fails(std::string_view expr) {
  auto script_res = script::translate(expr);
  if (!script_res.has_value()) return true;
  script::machine m;
  auto exec_res = m.execute(*script_res);
  return !exec_res.has_value();
}

}  // namespace

TEST(machine_tests, integer_literals) {
  EXPECT_EQ(eval_to_string("z'0"), "0");
  EXPECT_EQ(eval_to_string("z'42"), "42");
  EXPECT_EQ(eval_to_string("z'100"), "100");
}

TEST(machine_tests, real_literals) {
  EXPECT_EQ(eval_to_string("3.14"), "3.14000");
  EXPECT_EQ(eval_to_string("1.0"), "1.00000");
  EXPECT_EQ(eval_to_string(".5"), "0.50000");
  EXPECT_EQ(eval_to_string("0.0"), "0.00000");
}

TEST(machine_tests, integer_arithmetic) {
  EXPECT_EQ(eval_to_string("z'1+z'2"), "3");
  EXPECT_EQ(eval_to_string("z'5-z'3"), "2");
  EXPECT_EQ(eval_to_string("z'2*z'3"), "6");
  EXPECT_EQ(eval_to_string("z'1/z'3"), "0.33333");
  EXPECT_EQ(eval_to_string("z'10-z'3"), "7");
}

TEST(machine_tests, real_arithmetic) {
  EXPECT_EQ(eval_to_string("1.5+2.5"), "4.00000");
  EXPECT_EQ(eval_to_string("3.0-1.0"), "2.00000");
  EXPECT_EQ(eval_to_string("2.0*3.0"), "6.00000");
  EXPECT_EQ(eval_to_string("1.0/2.0"), "0.50000");
}

TEST(machine_tests, mixed_arithmetic) {
  EXPECT_EQ(eval_to_string("z'1+2.5"), "3.50000");
  EXPECT_EQ(eval_to_string("3.0+z'4"), "7.00000");
  EXPECT_EQ(eval_to_string("z'2*1.5"), "3.00000");
  EXPECT_EQ(eval_to_string("z'6/2.0"), "3.00000");
}

TEST(machine_tests, function_calls) {
  EXPECT_EQ(eval_to_string("sin(0)"), "0.00000");
  EXPECT_EQ(eval_to_string("cos(0)"), "1.00000");
  EXPECT_EQ(eval_to_string("exp(0)"), "1.00000");
  EXPECT_EQ(eval_to_string("ln(1)"), "0.00000");
  EXPECT_EQ(eval_to_string("log(1)"), "0.00000");
  EXPECT_EQ(eval_to_string("sqrt(4)"), "2.00000");
  EXPECT_EQ(eval_to_string("sqrt(2.25)"), "1.50000");
  EXPECT_EQ(eval_to_string("arctan(0)"), "0.00000");
}

TEST(machine_tests, nested_expressions) {
  EXPECT_EQ(eval_to_string("(1.0+2.0)*3.0"), "9.00000");
  EXPECT_EQ(eval_to_string("1.0+(2.0*3.0)"), "7.00000");
  EXPECT_EQ(eval_to_string("sqrt(1.0+3.0)"), "2.00000");
}

TEST(machine_tests, precedence) {
  EXPECT_EQ(eval_to_string("1.0+2.0*3.0"), "7.00000");
  EXPECT_EQ(eval_to_string("1.0*2.0+3.0"), "5.00000");
  EXPECT_EQ(eval_to_string("z'1+z'2*z'3"), "7");
  EXPECT_EQ(eval_to_string("z'1*z'2+z'3"), "5");
}

TEST(machine_tests, predefined_constants) {
  auto pi_str = eval_to_string("pi", 2);
  EXPECT_TRUE(pi_str.starts_with("3.14")) << "pi = " << pi_str;
  auto e_str = eval_to_string("e", 2);
  EXPECT_TRUE(e_str.starts_with("2.71") || e_str.starts_with("2.72")) << "e = " << e_str;
}

TEST(machine_tests, error_undefined_variable) {
  auto script_res = script::translate("foo");
  ASSERT_TRUE(script_res.has_value());
  script::machine m;
  auto exec_res = m.execute(*script_res);
  EXPECT_FALSE(exec_res.has_value());
}

TEST(machine_tests, error_wrong_arity) {
  auto script_res = script::translate("sin(1,2)");
  ASSERT_TRUE(script_res.has_value());
  script::machine m;
  auto exec_res = m.execute(*script_res);
  EXPECT_FALSE(exec_res.has_value());
}

TEST(machine_tests, error_undefined_function) {
  auto script_res = script::translate("foobar(1)");
  ASSERT_TRUE(script_res.has_value());
  script::machine m;
  auto exec_res = m.execute(*script_res);
  EXPECT_FALSE(exec_res.has_value());
}

TEST(machine_tests, trig_functions) {
  EXPECT_EQ(eval_to_string("tan(0)"), "0.00000");
  EXPECT_EQ(eval_to_string("sinh(0)"), "0.00000");
  EXPECT_EQ(eval_to_string("cosh(0)"), "1.00000");
  EXPECT_EQ(eval_to_string("tanh(0)"), "0.00000");
}

TEST(machine_tests, inv_trig_functions) {
  EXPECT_EQ(eval_to_string("arcsin(0)"), "0.00000");
  EXPECT_EQ(eval_to_string("arccos(1)"), "0.00000");
  EXPECT_EQ(eval_to_string("arcsinh(0)"), "0.00000");
  EXPECT_EQ(eval_to_string("arccosh(1)"), "0.00000");
  EXPECT_EQ(eval_to_string("arctanh(0)"), "0.00000");
}

TEST(machine_tests, binary_funcs) {
  EXPECT_EQ(eval_to_string("pow(2.0, 3.0)"), "8.00000");
  EXPECT_EQ(eval_to_string("pow(2.0, 0.0)"), "1.00000");
  EXPECT_EQ(eval_to_string("pow(4.0, 0.5)"), "2.00000");
}

TEST(machine_tests, integer_edge_cases) {
  EXPECT_EQ(eval_to_string("z'0+z'0"), "0");
  EXPECT_EQ(eval_to_string("z'0*z'5"), "0");
  EXPECT_EQ(eval_to_string("z'3-z'7"), "-4");
  EXPECT_EQ(eval_to_string("z'0/z'5"), "0.00000");
  EXPECT_EQ(eval_to_string("z'2/z'4"), "0.50000");
}

TEST(machine_tests, real_negative_results) {
  EXPECT_EQ(eval_to_string("1.0-3.0"), "-2.00000");
  EXPECT_EQ(eval_to_string("0.0-1.0"), "-1.00000");
  EXPECT_EQ(eval_to_string("2.0-5.0"), "-3.00000");
}

TEST(machine_tests, chained_operations) {
  EXPECT_EQ(eval_to_string("1.0+2.0+3.0"), "6.00000");
  EXPECT_EQ(eval_to_string("1.0-2.0-3.0"), "-4.00000");
  EXPECT_EQ(eval_to_string("2.0*3.0*4.0"), "24.00000");
  EXPECT_EQ(eval_to_string("z'1+z'2+z'3"), "6");
  EXPECT_EQ(eval_to_string("z'2*z'3*z'4"), "24");
}

TEST(machine_tests, nested_func_calls) {
  EXPECT_EQ(eval_to_string("sin(cos(0))", 10), "0.8414709848");
  EXPECT_EQ(eval_to_string("exp(ln(2.0))"), "2.00000");
  EXPECT_EQ(eval_to_string("sqrt(sqrt(16.0))"), "2.00000");
}

TEST(machine_tests, func_with_integer_arg) {
  EXPECT_EQ(eval_to_string("sin(z'0)"), "0.00000");
  EXPECT_EQ(eval_to_string("sqrt(z'9)"), "3.00000");
  EXPECT_EQ(eval_to_string("exp(z'0)"), "1.00000");
}

TEST(machine_tests, mixed_func_and_arithmetic) {
  EXPECT_EQ(eval_to_string("2.0*sin(0)+1.0"), "1.00000");
  EXPECT_EQ(eval_to_string("cos(0)*3.0"), "3.00000");
  EXPECT_EQ(eval_to_string("pow(sqrt(4.0), 2.0)"), "4.00000");
}

TEST(machine_tests, binary_func_arity_errors) {
  EXPECT_TRUE(eval_fails("pow(1.0)"));
  EXPECT_TRUE(eval_fails("pow(1.0,2.0,3.0)"));
  EXPECT_TRUE(eval_fails("log_base(1.0)"));
}

TEST(machine_tests, zero_arg_func_errors) {
  EXPECT_TRUE(eval_fails("sin()"));
  EXPECT_TRUE(eval_fails("cos()"));
  EXPECT_TRUE(eval_fails("sqrt()"));
}

}  // namespace epxut
