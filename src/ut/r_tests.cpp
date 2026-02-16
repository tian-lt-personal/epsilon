// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "chars.hpp"
#include "r.hpp"

// ut
#include "def.hpp"

namespace epxut {

TEST(r_tests, make_q) {
  auto one = stosz("1");
  {
    sz zero;
    auto q = epx::make_q(zero, one);  // 0/1
    auto eval1 = epx::coro::sync_get(q(1));
    EXPECT_TRUE(epx::is_zero(eval1));
    auto eval10 = epx::coro::sync_get(q(10));
    EXPECT_TRUE(epx::is_zero(eval10));
    auto eval100 = epx::coro::sync_get(q(100));
    EXPECT_TRUE(epx::is_zero(eval100));
  }
  {
    auto q = epx::make_q(one, one);  // 1/1
    auto eval1 = epx::coro::sync_get(q(1));
    EXPECT_EQ(stosz("4"), eval1);
    auto eval2 = epx::coro::sync_get(q(2));
    EXPECT_EQ(stosz("16"), eval2);
    auto eval10 = epx::coro::sync_get(q(10));
    EXPECT_EQ(stosz("1048576"), eval10);
  }
  {
    auto num = stosz("3");
    auto q = epx::make_q(one, num);  // 1/3
    auto eval1 = epx::coro::sync_get(q(1));
    EXPECT_EQ(one, eval1);
    auto eval5 = epx::coro::sync_get(q(5));
    EXPECT_EQ(stosz("341"), eval5);
  }

  {
    auto num = stosz("-3");
    auto q = epx::make_q(one, num);  // 1/-3
    auto eval1 = epx::coro::sync_get(q(1));
    EXPECT_EQ(stosz("-2"), eval1);
    auto eval5 = epx::coro::sync_get(q(5));
    EXPECT_EQ(stosz("-342"), eval5);
  }
}

TEST(r_tests, add) {
  {
    auto x = epx::make_q(stosz("0"), stosz("1"));
    auto y = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::add(x, y);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.00000", epx::to_string(expr, 5));
  }
  {
    auto x = epx::make_q(stosz("0"), stosz("-1"));
    auto y = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::add(x, y);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.00000", epx::to_string(expr, 5));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::add(x, y);
    EXPECT_EQ("2", epx::to_string(expr, 0));
    EXPECT_EQ("2.00000000", epx::to_string(expr, 8));
  }
  {
    auto x = epx::make_q(stosz("-1"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("-1"));
    auto expr = epx::add(x, y);
    EXPECT_EQ("-2", epx::to_string(expr, 0));
    EXPECT_EQ("-2.00000000", epx::to_string(expr, 8));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("1"));
    auto y = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::add(x, y);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.00000000", epx::to_string(expr, 8));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("100000000"));
    auto y = epx::make_q(stosz("99999997"), stosz("100000000"));
    auto expr = epx::add(epx::add(epx::add(x, y), x), x);
    EXPECT_EQ("1", epx::to_string(expr, 0));
    EXPECT_EQ("1.00000000", epx::to_string(expr, 8));
  }
}

}  // namespace epxut
