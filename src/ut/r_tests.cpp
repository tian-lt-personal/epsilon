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
    EXPECT_TRUE(epx::is_zero(q(1).get()));
    EXPECT_TRUE(epx::is_zero(q(10).get()));
    EXPECT_TRUE(epx::is_zero(q(100).get()));
  }
  {
    auto q = epx::make_q(one, one);  // 1/1
    EXPECT_EQ(stosz("4"), q(1).get());
    EXPECT_EQ(stosz("16"), q(2).get());
    EXPECT_EQ(stosz("1048576"), q(10).get());
  }
  {
    auto num = stosz("3");
    auto q = epx::make_q(one, num);  // 1/3
    EXPECT_EQ(one, q(1).get());
    EXPECT_EQ(stosz("341"), q(5).get());
  }

  {
    auto num = stosz("-3");
    auto q = epx::make_q(one, num);  // 1/-3
    EXPECT_EQ(stosz("-2"), q(1).get());
    EXPECT_EQ(stosz("-342"), q(5).get());
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

TEST(r_tests, opp) {
  {
    auto x = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::opp(x);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.00000", epx::to_string(expr, 5));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::opp(x);
    EXPECT_EQ("-1", epx::to_string(expr, 0));
    EXPECT_EQ("-1.00000", epx::to_string(expr, 5));
  }
  {
    auto x = epx::make_q(stosz("-1"), stosz("7"));
    auto expr = epx::opp(x);
    EXPECT_EQ("0.14285714285714285714286", epx::to_string(expr, 23));
  }
}

TEST(r_tests, msd) {
  {
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    EXPECT_EQ(10, epx::msd(zero, 10).get());
    EXPECT_EQ(1000, epx::msd(zero, 1000).get());
  }
  {
    auto r = epx::make_q(stosz("1"), stosz("1"));
    EXPECT_EQ(1, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("11"), stosz("10"));
    EXPECT_EQ(1, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("2"), stosz("1"));
    EXPECT_EQ(0, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("1"), stosz("2"));
    EXPECT_EQ(1, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("5"), stosz("1"));
    EXPECT_EQ(0, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("10"), stosz("1"));
    EXPECT_EQ(-1, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("127"), stosz("1"));
    EXPECT_EQ(-2, epx::msd(r, 10).get());
  }
  {
    auto r = epx::make_q(stosz("128"), stosz("1"));
    EXPECT_EQ(-3, epx::msd(r, 10).get());
  }
}

TEST(r_tests, mul) {
  {
    auto x = epx::make_q(stosz("0"), stosz("1"));
    auto y = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ(
        "0."
        "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000",
        epx::to_string(expr, 500));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("1.00000000000000000000", epx::to_string(expr, 20));
  }
  {
    auto x = epx::make_q(stosz("-2"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("-2.0000000000000000000000000000000000000000", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("3"));
    auto y = epx::make_q(stosz("3"), stosz("1"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("1.0000000000000000000000000000000000000000", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("7"), stosz("5"));
    auto y = epx::make_q(stosz("2"), stosz("9"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("0.3111111111111111111111111111111111111111", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("11"), stosz("7"));
    auto y = epx::make_q(stosz("1"), stosz("121"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("0.0129870129870129870129870129870129870130", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("10000000000"));
    auto y = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("0.00000000005000000000", epx::to_string(expr, 20));
  }
  {
    auto x = epx::make_q(stosz("827368917649287346"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("792873649187263413"));
    auto expr = epx::mul(x, y);
    EXPECT_EQ("1.04350664005214875176", epx::to_string(expr, 20));
  }
}

TEST(r_tests, inv) {
  {
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::inv(zero);
    EXPECT_THROW(expr(10).get(), epx::msd_overflow_error);
  }
  {
    auto r = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::inv(r);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  {
    auto r = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::inv(r);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
  {
    auto r = epx::make_q(stosz("1"), stosz("5"));
    auto expr = epx::inv(r);
    EXPECT_EQ("5.0000000000", epx::to_string(expr, 10));
  }
  {
    auto x = epx::make_q(stosz("341872364619278364918726351937561023895"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("782635487"));
    auto expr = epx::inv(epx::mul(x, y));
    EXPECT_EQ("0.0000000000000000000000000000022892622159", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("37"), stosz("1"));
    auto y = epx::make_q(stosz("1"), stosz("2141829365987369817236491872364918723641"));
    auto expr = epx::inv(epx::mul(x, y));
    EXPECT_EQ("57887280161820805871256537090943749287.59459459459459459459", epx::to_string(expr, 20));
  }
}

}  // namespace epxut
