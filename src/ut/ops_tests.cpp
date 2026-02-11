// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include <ops.hpp>
#include <z.hpp>

// ut
#include "def.hpp"

namespace epxut {

TEST(ops_tests, add) {
  lz zero, one = {.digits = {1}};
  auto res = zero + one;
  EXPECT_EQ(one, res);
}

TEST(ops_tests, sub) {
  lz zero, one = {.digits = {1}}, minus_one = {.digits = {1}, .sgn = epx::sign::negative};
  {
    auto res = zero - one;
    EXPECT_EQ(minus_one, res);
  }
  {
    auto res = zero - minus_one;
    EXPECT_EQ(one, res);
  }
}

TEST(ops_tests, mul) {
  sz z3 = {.digits = {3}}, z9 = {.digits = {9}};
  auto res = z3 * z3;
  EXPECT_EQ(z9, res);
}

TEST(ops_tests, div) {
  sz z6 = {.digits = {6}}, z3 = {.digits = {3}}, z2 = {.digits = {2}};
  auto res = z6 / z3;
  EXPECT_EQ(z2, res);
}

TEST(ops_tests, mod) {
  sz z7 = {.digits = {7}}, z3 = {.digits = {3}}, z4 = {.digits = {4}};
  auto res = z7 % z4;
  EXPECT_EQ(z3, res);
}

}  // namespace epxut
