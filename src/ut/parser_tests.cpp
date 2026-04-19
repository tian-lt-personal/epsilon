// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "lexer.hpp"
#include "parser.hpp"

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

}  // namespace epxut
