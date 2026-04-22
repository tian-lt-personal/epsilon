// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
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

TEST(parser_tests, empty_input) { auto res = script::translate(""); }
TEST(parser_tests, binary_ops) {
  {
    auto res = script::translate("1+2");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr +]\n  left={\n    {\n      real: 1\n    }\n  }\n  right={\n    {\n      real: 2\n    }\n  "
        "}\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("z'321-123");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr -]\n  left={\n    {\n      integer: 321\n    }\n  }\n  right={\n    {\n      real: 123\n    "
        "}\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("456*z'789");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr *]\n  left={\n    {\n      real: 456\n    }\n  }\n  right={\n    {\n      integer: 789\n    "
        "}\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("321/r'987");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr /]\n  left={\n    {\n      real: 321\n    }\n  }\n  right={\n    {\n      real: 987\n    }\n "
        " }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+2+3-4-5+6");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr +]\n  left={\n    {\n      [binop_expr -]\n      left={\n        {\n          [binop_expr "
        "-]\n          left={\n            {\n              [binop_expr +]\n              left={\n                {\n  "
        "                [binop_expr +]\n                  left={\n                    {\n                      real: "
        "1\n                    }\n                  }\n                  right={\n                    {\n             "
        "         real: 2\n                    }\n                  }\n                }\n              }\n            "
        "  right={\n                {\n                  real: 3\n                }\n              }\n            }\n  "
        "        }\n          right={\n            {\n              real: 4\n            }\n          }\n        }\n   "
        "   }\n      right={\n        {\n          real: 5\n        }\n      }\n    }\n  }\n  right={\n    {\n      "
        "real: 6\n    }\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2/3/4/5*6");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr *]\n  left={\n    {\n      [binop_expr /]\n      left={\n        {\n          [binop_expr "
        "/]\n          left={\n            {\n              [binop_expr /]\n              left={\n                {\n  "
        "                [binop_expr *]\n                  left={\n                    {\n                      real: "
        "1\n                    }\n                  }\n                  right={\n                    {\n             "
        "         real: 2\n                    }\n                  }\n                }\n              }\n            "
        "  right={\n                {\n                  real: 3\n                }\n              }\n            }\n  "
        "        }\n          right={\n            {\n              real: 4\n            }\n          }\n        }\n   "
        "   }\n      right={\n        {\n          real: 5\n        }\n      }\n    }\n  }\n  right={\n    {\n      "
        "real: 6\n    }\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1+2*3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr +]\n  left={\n    {\n      real: 1\n    }\n  }\n  right={\n    {\n      [binop_expr *]\n     "
        " left={\n        {\n          real: 2\n        }\n      }\n      right={\n        {\n          real: 3\n      "
        "  }\n      }\n    }\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2-3");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr -]\n  left={\n    {\n      [binop_expr *]\n      left={\n        {\n          real: 1\n      "
        "  }\n      }\n      right={\n        {\n          real: 2\n        }\n      }\n    }\n  }\n  right={\n    {\n "
        "     real: 3\n    }\n  }\n}\n",
        script::dump(*res));
  }
  {
    auto res = script::translate("1*2-3/z'4/5");
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(
        "{\n  [binop_expr -]\n  left={\n    {\n      [binop_expr *]\n      left={\n        {\n          real: 1\n      "
        "  }\n      }\n      right={\n        {\n          real: 2\n        }\n      }\n    }\n  }\n  right={\n    {\n "
        "     [binop_expr /]\n      left={\n        {\n          [binop_expr /]\n          left={\n            {\n     "
        "         real: 3\n            }\n          }\n          right={\n            {\n              integer: 4\n    "
        "        }\n          }\n        }\n      }\n      right={\n        {\n          real: 5\n        }\n      }\n "
        "   }\n  }\n}\n",
        script::dump(*res));
  }
}

}  // namespace epxut
