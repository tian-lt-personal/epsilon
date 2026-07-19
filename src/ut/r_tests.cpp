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
    EXPECT_TRUE(epx::is_zero(q.approx(1).get()));
    EXPECT_TRUE(epx::is_zero(q.approx(10).get()));
    EXPECT_TRUE(epx::is_zero(q.approx(100).get()));
  }
  {
    auto q = epx::make_q(one, one);  // 1/1
    EXPECT_EQ(stosz("4"), q.approx(1).get());
    EXPECT_EQ(stosz("16"), q.approx(2).get());
    EXPECT_EQ(stosz("1048576"), q.approx(10).get());
  }
  {
    auto num = stosz("3");
    auto q = epx::make_q(one, num);  // 1/3
    EXPECT_EQ(one, q.approx(1).get());
    EXPECT_EQ(stosz("341"), q.approx(5).get());
  }

  {
    auto num = stosz("-3");
    auto q = epx::make_q(one, num);  // 1/-3
    EXPECT_EQ(stosz("-2"), q.approx(1).get());
    EXPECT_EQ(stosz("-342"), q.approx(5).get());
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
    EXPECT_THROW(expr.approx(10).get(), epx::msd_overflow_error);
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

TEST(r_tests, root) {
  {
    auto x = epx::make_q(stosz("2"), stosz("3"));
    EXPECT_THROW(epx::root(x, 0).approx(2).get(), epx::kthroot_too_small_error);
    EXPECT_THROW(epx::root(x, 1).approx(2).get(), epx::kthroot_too_small_error);
  }
  {
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::root(zero, 2);
    EXPECT_EQ("0.00000", epx::to_string(expr, 5));
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  {
    auto x = epx::make_q(stosz("1"), stosz("1"));
    EXPECT_EQ("1.00000", epx::to_string(epx::root(x, 2), 5));
    EXPECT_EQ("1.000", epx::to_string(epx::root(x, 11), 3));
  }
  {
    auto x = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::root(x, 2);
    EXPECT_EQ("1.41421", epx::to_string(expr, 5));
    EXPECT_EQ("1.4142135624", epx::to_string(expr, 10));
    EXPECT_EQ("1.4142135624", epx::to_string(expr, 10));
    EXPECT_EQ("1.4142135623730950488016887242096980785697", epx::to_string(expr, 40));
  }
  {
    auto x = epx::make_q(stosz("34012224"), stosz("1000000"));
    auto expr = epx::root(x, 6);
    EXPECT_EQ("1.8000000000", epx::to_string(expr, 10));
  }
}

// Generated by AI — exp tests
TEST(r_tests, exp) {
  {
    // exp(0) = 1
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::exp(zero);
    EXPECT_EQ("1", epx::to_string(expr, 0));
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
    EXPECT_EQ("1.00000000000000000000", epx::to_string(expr, 20));
  }
  {
    // exp(1) = e ~ 2.71828182845904523536...
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::exp(one);
    EXPECT_EQ("2.71828182845904523536", epx::to_string(expr, 20));
  }
  {
    // exp(-1) ~ 0.36787944117144232159...
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::exp(neg_one);
    EXPECT_EQ("0.36787944117144232160", epx::to_string(expr, 20));
  }
  {
    // exp(1/2) ~ 1.64872127070012814684...
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::exp(half);
    EXPECT_EQ("1.6487212707", epx::to_string(expr, 10));
  }
  {
    // exp(10) ~ 22026.46579480671651695...
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto expr = epx::exp(ten);
    EXPECT_EQ("22026.4657948067", epx::to_string(expr, 10));
  }
  {
    // exp(-10) ~ 0.00004539992976248485...
    auto neg_ten = epx::make_q(stosz("-10"), stosz("1"));
    auto expr = epx::exp(neg_ten);
    EXPECT_EQ("0.0000453999", epx::to_string(expr, 10));
  }
  {
    // exp(x) * exp(-x) = 1 identity
    auto x = epx::make_q(stosz("3"), stosz("7"));
    auto ex = epx::exp(x);
    auto emx = epx::exp(epx::opp(x));
    auto product = epx::mul(ex, emx);
    EXPECT_EQ("1.000000000000000000000000000000", epx::to_string(product, 30));
  }
  {
    // exp(-1/2) ~ 0.60653065971263342360...
    auto neg_half = epx::make_q(stosz("-1"), stosz("2"));
    auto expr = epx::exp(neg_half);
    EXPECT_EQ("0.60653065971263342360", epx::to_string(expr, 20));
  }
  {
    // exp(2) ~ 7.38905609893065022723...
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::exp(two);
    EXPECT_EQ("7.38905609893065022723", epx::to_string(expr, 20));
  }
  {
    // exp(-2) ~ 0.13533528323661270232...
    auto neg_two = epx::make_q(stosz("-2"), stosz("1"));
    auto expr = epx::exp(neg_two);
    EXPECT_EQ("0.13533528323661269189", epx::to_string(expr, 20));
  }
  {
    // exp(1/10) ~ 1.10517091807564762481...
    auto tenth = epx::make_q(stosz("1"), stosz("10"));
    auto expr = epx::exp(tenth);
    EXPECT_EQ("1.10517091807564762481", epx::to_string(expr, 20));
  }
  {
    // exp(-1/10) ~ 0.90483741803595957316...
    auto neg_tenth = epx::make_q(stosz("-1"), stosz("10"));
    auto expr = epx::exp(neg_tenth);
    EXPECT_EQ("0.90483741803595957316", epx::to_string(expr, 20));
  }
  {
    // exp(1/100) ~ 1.01005016708416805755...
    auto hundredth = epx::make_q(stosz("1"), stosz("100"));
    auto expr = epx::exp(hundredth);
    EXPECT_EQ("1.01005016708416805754", epx::to_string(expr, 20));
  }
  {
    // exp(5) ~ 148.41315910357660143...
    auto five = epx::make_q(stosz("5"), stosz("1"));
    auto expr = epx::exp(five);
    EXPECT_EQ("148.4131591026", epx::to_string(expr, 10));
  }
  {
    // exp(-5) ~ 0.00673794699908546709...
    auto neg_five = epx::make_q(stosz("-5"), stosz("1"));
    auto expr = epx::exp(neg_five);
    EXPECT_EQ("0.0067379470", epx::to_string(expr, 10));
  }
  {
    // exp(1) at higher precision (30 digits)
    // e = 2.718281828459045235360287471352...
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::exp(one);
    EXPECT_EQ("2.718281828459045235360287471353", epx::to_string(expr, 30));
  }
  {
    // exp(a+b) = exp(a) * exp(b) identity
    auto a = epx::make_q(stosz("1"), stosz("3"));
    auto b = epx::make_q(stosz("2"), stosz("5"));
    auto ea = epx::exp(a);
    auto eb = epx::exp(b);
    auto eab = epx::exp(epx::add(a, b));
    auto product = epx::mul(ea, eb);
    EXPECT_EQ(epx::to_string(eab, 10), epx::to_string(product, 10));
  }
  {
    // very small positive: exp(1/1000) ~ 1.00100050016670833...
    auto small = epx::make_q(stosz("1"), stosz("1000"));
    auto expr = epx::exp(small);
    EXPECT_EQ("1.0010005002", epx::to_string(expr, 10));
  }
  {
    // very small negative: exp(-1/1000) ~ 0.99900049983337499...
    auto small_neg = epx::make_q(stosz("-1"), stosz("1000"));
    auto expr = epx::exp(small_neg);
    EXPECT_EQ("0.9990004998", epx::to_string(expr, 10));
  }
  {
    // exp(20) ~ 485165195.40979027797...
    auto x = epx::make_q(stosz("20"), stosz("1"));
    auto expr = epx::exp(x);
    EXPECT_EQ("485165195.4097902780", epx::to_string(expr, 10));
  }
  {
    // exp(50) ~ 5184705528587072464087.45...
    auto x = epx::make_q(stosz("50"), stosz("1"));
    auto expr = epx::exp(x);
    EXPECT_EQ("5184705528587072464087.4533229335", epx::to_string(expr, 10));
  }
  {
    // exp(100) ~ 2.688e43
    auto x = epx::make_q(stosz("100"), stosz("1"));
    auto expr = epx::exp(x);
    EXPECT_EQ("26881171418161354484126255515800135873611118.7737419224", epx::to_string(expr, 10));
  }
  {
    // exp(-20) ~ 2.061e-9
    auto x = epx::make_q(stosz("-20"), stosz("1"));
    auto expr = epx::exp(x);
    EXPECT_EQ("0.0000000021", epx::to_string(expr, 10));
  }
  {
    // exp(-50) ~ 1.929e-22
    auto x = epx::make_q(stosz("-50"), stosz("1"));
    auto expr = epx::exp(x);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  {
    // Branch: very negative x -> exp(x) rounds to 0 at internal precision
    auto neg100 = epx::make_q(stosz("-100"), stosz("1"));
    auto expr = epx::exp(neg100);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  {
    // Branch: n <= 0, neg_n == 0 — direct approx(0) call
    // exp(0) at n=0: |exp(0) - result/4^0| < 1, so result = 1
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::exp(zero);
    auto val = expr.approx(0).get();
    EXPECT_EQ("1", epx::to_string(val));
  }
  {
    // Branch: n <= 0, neg_n > 0 — direct approx(-1) call
    // exp(0) at n=-1: result satisfies |1 - result*4| < 4
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::exp(zero);
    auto val = expr.approx(-1).get();
    EXPECT_EQ("0", epx::to_string(val));
  }
  {
    // Branch: n <= 0, x_le_bound true — negative x at low precision
    // exp(-2) at n=0: x0=-2 <= bound=-1, returns 4^0=1
    // |exp(-2) - 1| = |0.135 - 1| = 0.865 < 1. Valid.
    auto neg2 = epx::make_q(stosz("-2"), stosz("1"));
    auto expr = epx::exp(neg2);
    auto val = expr.approx(0).get();
    EXPECT_EQ("0", epx::to_string(val));
  }
  {
    // Branch: n <= 0, x_le_bound false — falls through to general case
    // exp(1) at n=0: x0=1 > bound=-1, general case gives result ~3
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::exp(one);
    auto val = expr.approx(0).get();
    // exp(1) ~ 2.718, result = 2
    EXPECT_EQ("2", epx::to_string(val));
  }
  {
    // Branch: n < 0, general case — exp(1) at negative precision
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::exp(one);
    auto val = expr.approx(-1).get();
    // exp(1) ~ 2.718, at n=-1: result*4 ~ 2.718, result = 0
    EXPECT_EQ("0", epx::to_string(val));
  }
}

// Generated by AI — log tests (Algorithm 8)
TEST(r_tests, log) {
  {
    // ln(1) = 0
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::log(one);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
    EXPECT_EQ("0.00000000000000000000", epx::to_string(expr, 20));
  }
  {
    // ln(e) = 1
    auto e_val = epx::exp(epx::make_q(stosz("1"), stosz("1")));
    auto expr = epx::log(e_val);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
    EXPECT_EQ("1.00000000000000000000", epx::to_string(expr, 20));
  }
  {
    // ln(2) ~ 0.69314718055994530941...
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::log(two);
    EXPECT_EQ("0.6931471806", epx::to_string(expr, 10));
    EXPECT_EQ("0.69314718055994530942", epx::to_string(expr, 20));
  }
  {
    // ln(1/2) ~ -0.69314718055994530941...
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::log(half);
    EXPECT_EQ("-0.6931471806", epx::to_string(expr, 10));
  }
  {
    // ln(10) ~ 2.30258509299404568401...
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto expr = epx::log(ten);
    EXPECT_EQ("2.3025850930", epx::to_string(expr, 10));
  }
  {
    // ln(1/10) ~ -2.30258509299404568401...
    auto tenth = epx::make_q(stosz("1"), stosz("10"));
    auto expr = epx::log(tenth);
    EXPECT_EQ("-2.3025850930", epx::to_string(expr, 10));
  }
  {
    // Identity: ln(exp(x)) = x for x = 3/7
    auto x = epx::make_q(stosz("3"), stosz("7"));
    auto expr = epx::log(epx::exp(x));
    EXPECT_EQ("0.4285714286", epx::to_string(expr, 10));
  }
  {
    // Identity: ln(a*b) = ln(a) + ln(b)
    auto a = epx::make_q(stosz("3"), stosz("1"));
    auto b = epx::make_q(stosz("7"), stosz("1"));
    auto ln_ab = epx::log(epx::mul(a, b));
    auto ln_a_plus_ln_b = epx::add(epx::log(a), epx::log(b));
    EXPECT_EQ(epx::to_string(ln_ab, 50), epx::to_string(ln_a_plus_ln_b, 50));
  }
  {
    // ln(2) at higher precision (30 digits)
    // ln(2) = 0.693147180559945309417232121458...
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::log(two);
    EXPECT_EQ("0.693147180559945309417232121458", epx::to_string(expr, 30));
  }
  {
    // Edge: very small positive — ln(1/1000) ~ -6.90775527898213705205...
    auto small = epx::make_q(stosz("1"), stosz("1000"));
    auto expr = epx::log(small);
    EXPECT_EQ("-6.9077552790", epx::to_string(expr, 10));
  }
  {
    // Edge: very large — ln(1000) ~ 6.90775527898213705205...
    auto large = epx::make_q(stosz("1000"), stosz("1"));
    auto expr = epx::log(large);
    EXPECT_EQ("6.9077552790", epx::to_string(expr, 10));
  }
  {
    // ln(4) = 2 * ln(2) ~ 1.38629436111989...
    auto four = epx::make_q(stosz("4"), stosz("1"));
    auto expr = epx::log(four);
    EXPECT_EQ("1.3862943611", epx::to_string(expr, 10));
  }
  {
    // ln(0) diverges — msd overflow
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::log(zero_val);
    EXPECT_THROW(expr.approx(10).get(), epx::msd_overflow_error);
  }
  {
    // ln(negative) — non-positive log error
    auto neg = epx::make_q(stosz("-2"), stosz("1"));
    auto expr = epx::log(neg);
    EXPECT_THROW(expr.approx(10).get(), epx::non_positive_log_error);
  }
}

// Generated by AI — arctan tests (Algorithm 9)
TEST(r_tests, arctan) {
  {
    // arctan(0) = 0
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::arctan(zero_val);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  {
    // arctan(1) = pi/4 ~ 0.7853981633974483096156608...
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arctan(one);
    EXPECT_EQ("0.7853981634", epx::to_string(expr, 10));
    EXPECT_EQ("0.78539816339744830962", epx::to_string(expr, 20));
  }
  {
    // arctan(-1) = -pi/4
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::arctan(neg_one);
    EXPECT_EQ("-0.7853981634", epx::to_string(expr, 10));
  }
  {
    // arctan(sqrt(3)) = pi/3 ~ 1.0471975511965977461542144...
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto expr = epx::arctan(sqrt3);
    EXPECT_EQ("1.0471975512", epx::to_string(expr, 10));
  }
  {
    // arctan(1/sqrt(3)) = pi/6 ~ 0.5235987755982988730771072...
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    auto expr = epx::arctan(inv_sqrt3);
    EXPECT_EQ("0.5235987756", epx::to_string(expr, 10));
  }
  {
    // arctan(1000) ~ pi/2 - 0.001 ~ 1.5697963271...
    auto large = epx::make_q(stosz("1000"), stosz("1"));
    auto expr = epx::arctan(large);
    EXPECT_EQ("1.5697963271", epx::to_string(expr, 10));
  }
  {
    // arctan(1/1000) ~ 0.000999999666...
    auto small = epx::make_q(stosz("1"), stosz("1000"));
    auto expr = epx::arctan(small);
    EXPECT_EQ("0.0009999997", epx::to_string(expr, 10));
  }
  {
    // Identity: arctan(x) + arctan(1/x) = pi/2 for x > 0
    auto x = epx::make_q(stosz("3"), stosz("7"));
    auto atan_x = epx::arctan(x);
    auto atan_inv_x = epx::arctan(epx::inv(x));
    auto sum = epx::add(atan_x, atan_inv_x);
    // pi/2 ~ 1.5707963267948966...
    EXPECT_EQ("1.5707963268", epx::to_string(sum, 10));
  }
  {
    // High precision (30 digits): arctan(1) = pi/4
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arctan(one);
    // pi/4 = 0.785398163397448309615660845819...
    EXPECT_EQ("0.785398163397448309615660845820", epx::to_string(expr, 30));
  }
  {
    // Gauss's formula: 12*arctan(1/18) + 8*arctan(1/57) - 5*arctan(1/239) = pi/4
    auto a18 = epx::arctan(epx::make_q(stosz("1"), stosz("18")));
    auto a57 = epx::arctan(epx::make_q(stosz("1"), stosz("57")));
    auto a239 = epx::arctan(epx::make_q(stosz("1"), stosz("239")));
    auto twelve = epx::make_q(stosz("12"), stosz("1"));
    auto eight = epx::make_q(stosz("8"), stosz("1"));
    auto five = epx::make_q(stosz("5"), stosz("1"));
    auto gauss = epx::add(epx::add(epx::mul(twelve, a18), epx::mul(eight, a57)), epx::opp(epx::mul(five, a239)));
    // Should equal pi/4
    EXPECT_EQ("0.7853981634", epx::to_string(gauss, 10));
  }
  {
    // arctan(-1000) ~ -pi/2 + 0.001 ~ -1.5697963271...
    auto neg_large = epx::make_q(stosz("-1000"), stosz("1"));
    auto expr = epx::arctan(neg_large);
    EXPECT_EQ("-1.5697963271", epx::to_string(expr, 10));
  }
  {
    // Identity: arctan(a) + arctan(b) = arctan((a+b)/(1-ab)) when ab < 1
    // Use a = 1/3, b = 1/5: ab = 1/15 < 1
    // (a+b)/(1-ab) = (8/15)/(14/15) = 8/14 = 4/7
    auto a = epx::make_q(stosz("1"), stosz("3"));
    auto b = epx::make_q(stosz("1"), stosz("5"));
    auto lhs = epx::add(epx::arctan(a), epx::arctan(b));
    auto rhs = epx::arctan(epx::make_q(stosz("4"), stosz("7")));
    EXPECT_EQ(epx::to_string(lhs, 20), epx::to_string(rhs, 20));
  }
}

// Generated by AI — arctan tests using mz (uint16_t container)
TEST(r_tests, arctan_mz) {
  auto one = epx::make_q(stomz("1"), stomz("1"));
  auto expr = epx::arctan(one);
  EXPECT_EQ("0.78539816339744830962", epx::to_string(expr, 20));

  auto neg_one = epx::make_q(stomz("-1"), stomz("1"));
  EXPECT_EQ("-0.7853981634", epx::to_string(epx::arctan(neg_one), 10));

  // arctan(sqrt(3)) = pi/3 — triggers |y| > 1 branch with denom/num > 1/2
  auto three = epx::make_q(stomz("3"), stomz("1"));
  auto sqrt3 = epx::root(three, 2);
  EXPECT_EQ("1.0471975512", epx::to_string(epx::arctan(sqrt3), 10));

  // arctan(1000) — triggers |y| > 1 branch with denom/num <= 1/2
  auto large = epx::make_q(stomz("1000"), stomz("1"));
  EXPECT_EQ("1.5697963271", epx::to_string(epx::arctan(large), 10));

  // arctan(1/sqrt(3)) = pi/6 — triggers |y| < 1, |y| > 1/2 branch
  auto inv_sqrt3 = epx::inv(sqrt3);
  EXPECT_EQ("0.5235987756", epx::to_string(epx::arctan(inv_sqrt3), 10));

  // arctan(1/1000) — triggers |y| <= 1/2 direct series
  auto small = epx::make_q(stomz("1"), stomz("1000"));
  EXPECT_EQ("0.0009999997", epx::to_string(epx::arctan(small), 10));
}

// Generated by AI — arctan tests using lz (uint32_t container)
TEST(r_tests, arctan_lz) {
  auto one = epx::make_q(stolz("1"), stolz("1"));
  auto expr = epx::arctan(one);
  EXPECT_EQ("0.78539816339744830962", epx::to_string(expr, 20));

  auto three = epx::make_q(stolz("3"), stolz("1"));
  auto sqrt3 = epx::root(three, 2);
  EXPECT_EQ("1.0471975512", epx::to_string(epx::arctan(sqrt3), 10));

  // arctan(1000) — triggers |y| > 1 branch with denom/num <= 1/2
  auto large = epx::make_q(stolz("1000"), stolz("1"));
  EXPECT_EQ("1.5697963271", epx::to_string(epx::arctan(large), 10));

  // arctan(1/sqrt(3)) = pi/6 — triggers |y| < 1, |y| > 1/2 branch
  auto inv_sqrt3 = epx::inv(sqrt3);
  EXPECT_EQ("0.5235987756", epx::to_string(epx::arctan(inv_sqrt3), 10));
}

// Generated by AI — pi precision tests via 4*arctan(1)
TEST(r_tests, pi) {
  auto one = epx::make_q(stosz("1"), stosz("1"));
  auto four = epx::make_q(stosz("4"), stosz("1"));
  auto pi = epx::mul(four, epx::arctan(one));

  // pi = 3.14159265358979323846264338327950288419716939937510...
  EXPECT_EQ("3", epx::to_string(pi, 0));
  EXPECT_EQ("3.1415926536", epx::to_string(pi, 10));
  EXPECT_EQ("3.14159265358979323846", epx::to_string(pi, 20));
  EXPECT_EQ("3.141592653589793238462643383280", epx::to_string(pi, 30));
  EXPECT_EQ("3.1415926535897932384626433832795028841972", epx::to_string(pi, 40));
  EXPECT_EQ("3.14159265358979323846264338327950288419716939937511", epx::to_string(pi, 50));
}

// Generated by AI — pi tests using mz (uint16_t container)
TEST(r_tests, pi_mz) {
  auto one = epx::make_q(stomz("1"), stomz("1"));
  auto four = epx::make_q(stomz("4"), stomz("1"));
  auto pi = epx::mul(four, epx::arctan(one));

  EXPECT_EQ("3.1415926536", epx::to_string(pi, 10));
  EXPECT_EQ("3.14159265358979323846", epx::to_string(pi, 20));
  EXPECT_EQ("3.141592653589793238462643383280", epx::to_string(pi, 30));
  EXPECT_EQ("3.1415926535897932384626433832795028841972", epx::to_string(pi, 40));
}

// Generated by AI — pi tests using lz (uint32_t container)
TEST(r_tests, pi_lz) {
  auto one = epx::make_q(stolz("1"), stolz("1"));
  auto four = epx::make_q(stolz("4"), stolz("1"));
  auto pi = epx::mul(four, epx::arctan(one));

  EXPECT_EQ("3.1415926536", epx::to_string(pi, 10));
  EXPECT_EQ("3.14159265358979323846", epx::to_string(pi, 20));
  EXPECT_EQ("3.141592653589793238462643383280", epx::to_string(pi, 30));
  EXPECT_EQ("3.1415926535897932384626433832795028841972", epx::to_string(pi, 40));
}

// Generated by AI — sine function tests
TEST(r_tests, sin) {
  // sin(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::sin(zero_val);
    EXPECT_EQ("0", epx::to_string(expr, 0));
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // sin(pi) = 0
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto four = epx::make_q(stosz("4"), stosz("1"));
    auto pi = epx::mul(four, epx::arctan(one));
    auto expr = epx::sin(pi);
    EXPECT_EQ("0", epx::to_string(expr, 0));
  }
  // sin(pi/2) = 1
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto half_pi = epx::mul(two, epx::arctan(one));
    auto expr = epx::sin(half_pi);
    EXPECT_EQ("1", epx::to_string(expr, 0));
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // sin(pi/6) = 0.5
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    // pi/6 = arctan(1/sqrt(3))
    auto pi_6 = epx::arctan(inv_sqrt3);
    auto expr = epx::sin(pi_6);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
  // sin(-pi/2) = -1
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto neg_two = epx::make_q(stosz("-2"), stosz("1"));
    auto neg_half_pi = epx::mul(neg_two, epx::arctan(one));
    auto expr = epx::sin(neg_half_pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
}

// Generated by AI — cosine function tests
TEST(r_tests, cos) {
  // cos(0) = 1
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::cos(zero_val);
    EXPECT_EQ("1", epx::to_string(expr, 0));
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // cos(pi) = -1
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto four = epx::make_q(stosz("4"), stosz("1"));
    auto pi = epx::mul(four, epx::arctan(one));
    auto expr = epx::cos(pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
  // cos(pi/3) = 0.5
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    // pi/3 = arctan(sqrt(3))
    auto pi_3 = epx::arctan(sqrt3);
    auto expr = epx::cos(pi_3);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
  // sin^2(x) + cos^2(x) = 1
  {
    auto x = epx::make_q(stosz("7"), stosz("10"));  // x = 0.7
    auto sx = epx::sin(x);
    auto cx = epx::cos(x);
    auto sum = epx::add(epx::mul(sx, sx), epx::mul(cx, cx));
    EXPECT_EQ("1.0000000000", epx::to_string(sum, 10));
  }
}

// Generated by AI — sin boundary tests near zone classification thresholds
TEST(r_tests, sin_boundary) {
  auto one = epx::make_q(stosz("1"), stosz("1"));

  // sin(pi/4) = sqrt(2)/2 — direct series path (zk < pi_k)
  {
    auto pi_4 = epx::arctan(one);
    auto expr = epx::sin(pi_4);
    EXPECT_EQ("0.7071067812", epx::to_string(expr, 10));
  }
  // sin(2*pi/3) = sqrt(3)/2 — series with pi-reduction (zk > pi_k)
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto two_pi_3 = epx::mul(two, epx::arctan(sqrt3));
    auto expr = epx::sin(two_pi_3);
    EXPECT_EQ("0.8660254038", epx::to_string(expr, 10));
  }
  // sin(5*pi/4) = -sqrt(2)/2 — p=1 odd, negated result
  {
    auto five = epx::make_q(stosz("5"), stosz("1"));
    auto five_pi_4 = epx::mul(five, epx::arctan(one));
    auto expr = epx::sin(five_pi_4);
    EXPECT_EQ("-0.7071067812", epx::to_string(expr, 10));
  }
  // sin(3*pi/2) = -1 — p=1, zk = pi/2 hits zone_pik, then negated
  {
    auto six = epx::make_q(stosz("6"), stosz("1"));
    auto three_pi_2 = epx::mul(six, epx::arctan(one));
    auto expr = epx::sin(three_pi_2);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
  // sin(2*pi) = 0 — p=2, zk near 0 hits zone_0
  {
    auto eight = epx::make_q(stosz("8"), stosz("1"));
    auto two_pi = epx::mul(eight, epx::arctan(one));
    auto expr = epx::sin(two_pi);
    EXPECT_EQ("0", epx::to_string(expr, 0));
  }
  // sin(1) — rational input, direct series
  {
    auto expr = epx::sin(one);
    EXPECT_EQ("0.8414709848", epx::to_string(expr, 10));
  }
  // sin(3) — near pi, tests series reduction since 3 > pi/2
  {
    auto three_val = epx::make_q(stosz("3"), stosz("1"));
    auto expr = epx::sin(three_val);
    EXPECT_EQ("0.1411200081", epx::to_string(expr, 10));
  }
  // sin(1/1000) — very small angle, near zone_0 boundary
  {
    auto small = epx::make_q(stosz("1"), stosz("1000"));
    auto expr = epx::sin(small);
    EXPECT_EQ("0.0009999998", epx::to_string(expr, 10));
  }
  // sin(100) — large argument, stress-tests argument reduction
  {
    auto big = epx::make_q(stosz("100"), stosz("1"));
    auto expr = epx::sin(big);
    EXPECT_EQ("-0.5063656411", epx::to_string(expr, 10));
  }
}

// Generated by AI — cos boundary tests
TEST(r_tests, cos_boundary) {
  auto one = epx::make_q(stosz("1"), stosz("1"));

  // cos(pi/4) = sqrt(2)/2
  {
    auto pi_4 = epx::arctan(one);
    auto expr = epx::cos(pi_4);
    EXPECT_EQ("0.7071067812", epx::to_string(expr, 10));
  }
  // cos(pi/6) = sqrt(3)/2
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    auto pi_6 = epx::arctan(inv_sqrt3);
    auto expr = epx::cos(pi_6);
    EXPECT_EQ("0.8660254038", epx::to_string(expr, 10));
  }
  // cos(2*pi/3) = -0.5
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto two_pi_3 = epx::mul(two, epx::arctan(sqrt3));
    auto expr = epx::cos(two_pi_3);
    EXPECT_EQ("-0.5000000000", epx::to_string(expr, 10));
  }
  // cos(2*pi) = 1
  {
    auto eight = epx::make_q(stosz("8"), stosz("1"));
    auto two_pi = epx::mul(eight, epx::arctan(one));
    auto expr = epx::cos(two_pi);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // cos(3*pi/2) = 0
  {
    auto six = epx::make_q(stosz("6"), stosz("1"));
    auto three_pi_2 = epx::mul(six, epx::arctan(one));
    auto expr = epx::cos(three_pi_2);
    EXPECT_EQ("0", epx::to_string(expr, 0));
  }
  // cos(1) — rational input
  {
    auto expr = epx::cos(one);
    EXPECT_EQ("0.5403023059", epx::to_string(expr, 10));
  }
  // sin^2(3) + cos^2(3) = 1 — identity at a near-pi boundary point
  {
    auto three_val = epx::make_q(stosz("3"), stosz("1"));
    auto sx = epx::sin(three_val);
    auto cx = epx::cos(three_val);
    auto sum = epx::add(epx::mul(sx, sx), epx::mul(cx, cx));
    EXPECT_EQ("1.0000000000", epx::to_string(sum, 10));
  }
}

// Generated by AI — sin tests using mz (uint16_t container)
TEST(r_tests, sin_mz) {
  // sin(0) = 0
  {
    auto zero_val = epx::make_q(stomz("0"), stomz("1"));
    auto expr = epx::sin(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // sin(pi/2) = 1
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    auto two = epx::make_q(stomz("2"), stomz("1"));
    auto half_pi = epx::mul(two, epx::arctan(one));
    auto expr = epx::sin(half_pi);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // sin(pi/6) = 0.5
  {
    auto three = epx::make_q(stomz("3"), stomz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    auto pi_6 = epx::arctan(inv_sqrt3);
    auto expr = epx::sin(pi_6);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
  // sin(-pi/2) = -1
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    auto neg_two = epx::make_q(stomz("-2"), stomz("1"));
    auto neg_half_pi = epx::mul(neg_two, epx::arctan(one));
    auto expr = epx::sin(neg_half_pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
}

// Generated by AI — sin tests using lz (uint32_t container)
TEST(r_tests, sin_lz) {
  // sin(pi/2) = 1
  {
    auto one = epx::make_q(stolz("1"), stolz("1"));
    auto two = epx::make_q(stolz("2"), stolz("1"));
    auto half_pi = epx::mul(two, epx::arctan(one));
    auto expr = epx::sin(half_pi);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // sin(pi/6) = 0.5
  {
    auto three = epx::make_q(stolz("3"), stolz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    auto pi_6 = epx::arctan(inv_sqrt3);
    auto expr = epx::sin(pi_6);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
  // sin(-pi/2) = -1
  {
    auto one = epx::make_q(stolz("1"), stolz("1"));
    auto neg_two = epx::make_q(stolz("-2"), stolz("1"));
    auto neg_half_pi = epx::mul(neg_two, epx::arctan(one));
    auto expr = epx::sin(neg_half_pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
}

// Generated by AI — cos tests using mz (uint16_t container)
TEST(r_tests, cos_mz) {
  // cos(0) = 1
  {
    auto zero_val = epx::make_q(stomz("0"), stomz("1"));
    auto expr = epx::cos(zero_val);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // cos(pi) = -1
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    auto four = epx::make_q(stomz("4"), stomz("1"));
    auto pi = epx::mul(four, epx::arctan(one));
    auto expr = epx::cos(pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
  // cos(pi/3) = 0.5
  {
    auto three = epx::make_q(stomz("3"), stomz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto pi_3 = epx::arctan(sqrt3);
    auto expr = epx::cos(pi_3);
    EXPECT_EQ("0.5000000000", epx::to_string(expr, 10));
  }
}

// Generated by AI — cos tests using lz (uint32_t container)
TEST(r_tests, cos_lz) {
  // cos(0) = 1
  {
    auto zero_val = epx::make_q(stolz("0"), stolz("1"));
    auto expr = epx::cos(zero_val);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // cos(pi) = -1
  {
    auto one = epx::make_q(stolz("1"), stolz("1"));
    auto four = epx::make_q(stolz("4"), stolz("1"));
    auto pi = epx::mul(four, epx::arctan(one));
    auto expr = epx::cos(pi);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
  // sin^2(x) + cos^2(x) = 1
  {
    auto x = epx::make_q(stolz("7"), stolz("10"));
    auto sx = epx::sin(x);
    auto cx = epx::cos(x);
    auto sum = epx::add(epx::mul(sx, sx), epx::mul(cx, cx));
    EXPECT_EQ("1.0000000000", epx::to_string(sum, 10));
  }
}

// Generated by AI — Section 4.3.7: Other elementary functions tests

TEST(r_tests, tan) {
  // tan(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::tan(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // tan(pi/4) = 1
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto pi_4 = epx::arctan(one);
    auto expr = epx::tan(pi_4);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // tan(pi/6) = 1/sqrt(3) ≈ 0.5773502692
  {
    auto three = epx::make_q(stosz("3"), stosz("1"));
    auto sqrt3 = epx::root(three, 2);
    auto inv_sqrt3 = epx::inv(sqrt3);
    auto pi_6 = epx::arctan(inv_sqrt3);
    auto expr = epx::tan(pi_6);
    EXPECT_EQ("0.5773502692", epx::to_string(expr, 10));
  }
  // tan(-pi/4) = -1
  {
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto neg_pi_4 = epx::arctan(neg_one);
    auto expr = epx::tan(neg_pi_4);
    EXPECT_EQ("-1.0000000000", epx::to_string(expr, 10));
  }
}

TEST(r_tests, arcsin) {
  // arcsin(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::arcsin(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // arcsin(1/2) = pi/6 ≈ 0.5235987756
  {
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::arcsin(half);
    EXPECT_EQ("0.5235987756", epx::to_string(expr, 10));
  }
  // arcsin(sqrt(2)/2) = pi/4 ≈ 0.7853981634
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto sqrt2 = epx::root(two, 2);
    auto val = epx::mul(sqrt2, epx::inv(epx::make_q(stosz("2"), stosz("1"))));
    auto expr = epx::arcsin(val);
    EXPECT_EQ("0.7853981634", epx::to_string(expr, 10));
  }
  // arcsin(-1/2) = -pi/6
  {
    auto neg_half = epx::make_q(stosz("-1"), stosz("2"));
    auto expr = epx::arcsin(neg_half);
    EXPECT_EQ("-0.5235987756", epx::to_string(expr, 10));
  }
  // arcsin(1) = pi/2 ≈ 1.5707963268
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arcsin(one);
    EXPECT_EQ("1.5707963268", epx::to_string(expr, 10));
  }
  // arcsin(-1) = -pi/2 ≈ -1.5707963268
  {
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::arcsin(neg_one);
    EXPECT_EQ("-1.5707963268", epx::to_string(expr, 10));
  }
}

TEST(r_tests, arccos) {
  // arccos(1) = 0
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arccos(one);
    EXPECT_EQ("0", epx::to_string(expr, 0));
  }
  // arccos(1/2) = pi/3 ≈ 1.0471975512
  {
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::arccos(half);
    EXPECT_EQ("1.0471975512", epx::to_string(expr, 10));
  }
  // arccos(sqrt(2)/2) = pi/4 ≈ 0.7853981634
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto sqrt2 = epx::root(two, 2);
    auto val = epx::mul(sqrt2, epx::inv(epx::make_q(stosz("2"), stosz("1"))));
    auto expr = epx::arccos(val);
    EXPECT_EQ("0.7853981634", epx::to_string(expr, 10));
  }
  // arccos(0) = pi/2 ≈ 1.5707963268
  {
    auto zero = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::arccos(zero);
    EXPECT_EQ("1.5707963268", epx::to_string(expr, 10));
  }
  // arcsin(x) + arccos(x) = pi/2
  {
    auto val = epx::make_q(stosz("3"), stosz("7"));
    auto sum = epx::add(epx::arcsin(val), epx::arccos(val));
    EXPECT_EQ("1.5707963268", epx::to_string(sum, 10));
  }
}

TEST(r_tests, sinh) {
  // sinh(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::sinh(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // sinh(1) ≈ 1.1752011936
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::sinh(one);
    EXPECT_EQ("1.1752011936", epx::to_string(expr, 10));
  }
  // sinh(-1) = -sinh(1) (odd function)
  {
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::sinh(neg_one);
    EXPECT_EQ("-1.1752011936", epx::to_string(expr, 10));
  }
}

TEST(r_tests, cosh) {
  // cosh(0) = 1
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::cosh(zero_val);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
  // cosh(1) ≈ 1.5430806348
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::cosh(one);
    EXPECT_EQ("1.5430806348", epx::to_string(expr, 10));
  }
  // cosh(-1) = cosh(1) (even function)
  {
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::cosh(neg_one);
    EXPECT_EQ("1.5430806348", epx::to_string(expr, 10));
  }
  // cosh^2(x) - sinh^2(x) = 1
  {
    auto x = epx::make_q(stosz("2"), stosz("1"));
    auto ch = epx::cosh(x);
    auto sh = epx::sinh(x);
    auto diff = epx::add(epx::mul(ch, ch), epx::opp(epx::mul(sh, sh)));
    EXPECT_EQ("1.0000000000", epx::to_string(diff, 10));
  }
}

TEST(r_tests, tanh) {
  // tanh(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::tanh(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // tanh(1) ≈ 0.7615941560
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::tanh(one);
    EXPECT_EQ("0.7615941560", epx::to_string(expr, 10));
  }
  // tanh(-1) = -tanh(1) (odd function)
  {
    auto neg_one = epx::make_q(stosz("-1"), stosz("1"));
    auto expr = epx::tanh(neg_one);
    EXPECT_EQ("-0.7615941560", epx::to_string(expr, 10));
  }
}

TEST(r_tests, arcsinh) {
  // arcsinh(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::arcsinh(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // arcsinh(1) ≈ 0.8813735870
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arcsinh(one);
    EXPECT_EQ("0.8813735870", epx::to_string(expr, 10));
  }
  // sinh(arcsinh(x)) = x (roundtrip)
  {
    auto val = epx::make_q(stosz("3"), stosz("5"));
    auto expr = epx::sinh(epx::arcsinh(val));
    EXPECT_EQ("0.6000000000", epx::to_string(expr, 10));
  }
}

TEST(r_tests, arccosh) {
  // arccosh(1) = 0
  {
    auto one = epx::make_q(stosz("1"), stosz("1"));
    auto expr = epx::arccosh(one);
    EXPECT_EQ("0", epx::to_string(expr, 0));
  }
  // arccosh(2) ≈ 1.3169578969
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::arccosh(two);
    EXPECT_EQ("1.3169578969", epx::to_string(expr, 10));
  }
  // cosh(arccosh(x)) = x (roundtrip)
  {
    auto val = epx::make_q(stosz("3"), stosz("1"));
    auto expr = epx::cosh(epx::arccosh(val));
    EXPECT_EQ("3.0000000000", epx::to_string(expr, 10));
  }
}

TEST(r_tests, arctanh) {
  // arctanh(0) = 0
  {
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::arctanh(zero_val);
    EXPECT_EQ("0.0000000000", epx::to_string(expr, 10));
  }
  // arctanh(1/2) ≈ 0.5493061443
  {
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::arctanh(half);
    EXPECT_EQ("0.5493061443", epx::to_string(expr, 10));
  }
  // tanh(arctanh(x)) = x (roundtrip)
  {
    auto val = epx::make_q(stosz("1"), stosz("3"));
    auto expr = epx::tanh(epx::arctanh(val));
    EXPECT_EQ("0.3333333333", epx::to_string(expr, 10));
  }
}

TEST(r_tests, pow) {
  // 2^10 = 1024
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto expr = epx::pow(two, ten);
    EXPECT_EQ("1024", epx::to_string(expr, 0));
  }
  // 2^0.5 = sqrt(2) ≈ 1.4142135624
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto half = epx::make_q(stosz("1"), stosz("2"));
    auto expr = epx::pow(two, half);
    EXPECT_EQ("1.4142135624", epx::to_string(expr, 10));
  }
  // 10^0 = 1
  {
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto zero_val = epx::make_q(stosz("0"), stosz("1"));
    auto expr = epx::pow(ten, zero_val);
    EXPECT_EQ("1.0000000000", epx::to_string(expr, 10));
  }
}

TEST(r_tests, log_base) {
  // log_2(8) = 3
  {
    auto eight = epx::make_q(stosz("8"), stosz("1"));
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto expr = epx::log_base(eight, two);
    EXPECT_EQ("3.0000000000", epx::to_string(expr, 10));
  }
  // log_10(100) = 2
  {
    auto hundred = epx::make_q(stosz("100"), stosz("1"));
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto expr = epx::log_base(hundred, ten);
    EXPECT_EQ("2.0000000000", epx::to_string(expr, 10));
  }
  // log_10(2) ≈ 0.3010299957
  {
    auto two = epx::make_q(stosz("2"), stosz("1"));
    auto ten = epx::make_q(stosz("10"), stosz("1"));
    auto expr = epx::log_base(two, ten);
    EXPECT_EQ("0.3010299957", epx::to_string(expr, 10));
  }
}

// Generated by AI — mz/lz tests for section 4.3.7 functions
TEST(r_tests, elementary_mz) {
  // tan(pi/4) = 1
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    auto pi_4 = epx::arctan(one);
    EXPECT_EQ("1.0000000000", epx::to_string(epx::tan(pi_4), 10));
  }
  // arcsin(1/2) = pi/6
  {
    auto half = epx::make_q(stomz("1"), stomz("2"));
    EXPECT_EQ("0.5235987756", epx::to_string(epx::arcsin(half), 10));
  }
  // arccos(1/2) = pi/3
  {
    auto half = epx::make_q(stomz("1"), stomz("2"));
    EXPECT_EQ("1.0471975512", epx::to_string(epx::arccos(half), 10));
  }
  // sinh(1) ≈ 1.1752011936
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    EXPECT_EQ("1.1752011936", epx::to_string(epx::sinh(one), 10));
  }
  // cosh(1) ≈ 1.5430806348
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    EXPECT_EQ("1.5430806348", epx::to_string(epx::cosh(one), 10));
  }
  // tanh(1) ≈ 0.7615941560
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    EXPECT_EQ("0.7615941560", epx::to_string(epx::tanh(one), 10));
  }
  // arcsinh(1) ≈ 0.8813735870
  {
    auto one = epx::make_q(stomz("1"), stomz("1"));
    EXPECT_EQ("0.8813735870", epx::to_string(epx::arcsinh(one), 10));
  }
  // arccosh(2) ≈ 1.3169578969
  {
    auto two = epx::make_q(stomz("2"), stomz("1"));
    EXPECT_EQ("1.3169578969", epx::to_string(epx::arccosh(two), 10));
  }
  // arctanh(1/2) ≈ 0.5493061443
  {
    auto half = epx::make_q(stomz("1"), stomz("2"));
    EXPECT_EQ("0.5493061443", epx::to_string(epx::arctanh(half), 10));
  }
  // pow(2, 0.5) = sqrt(2)
  {
    auto two = epx::make_q(stomz("2"), stomz("1"));
    auto half = epx::make_q(stomz("1"), stomz("2"));
    EXPECT_EQ("1.4142135624", epx::to_string(epx::pow(two, half), 10));
  }
  // log_2(8) = 3
  {
    auto eight = epx::make_q(stomz("8"), stomz("1"));
    auto two = epx::make_q(stomz("2"), stomz("1"));
    EXPECT_EQ("3.0000000000", epx::to_string(epx::log_base(eight, two), 10));
  }
}

TEST(r_tests, elementary_lz) {
  // tan(pi/4) = 1
  {
    auto one = epx::make_q(stolz("1"), stolz("1"));
    auto pi_4 = epx::arctan(one);
    EXPECT_EQ("1.0000000000", epx::to_string(epx::tan(pi_4), 10));
  }
  // arcsin(1/2) = pi/6
  {
    auto half = epx::make_q(stolz("1"), stolz("2"));
    EXPECT_EQ("0.5235987756", epx::to_string(epx::arcsin(half), 10));
  }
  // sinh(1) ≈ 1.1752011936
  {
    auto one = epx::make_q(stolz("1"), stolz("1"));
    EXPECT_EQ("1.1752011936", epx::to_string(epx::sinh(one), 10));
  }
  // cosh^2 - sinh^2 = 1
  {
    auto x = epx::make_q(stolz("2"), stolz("1"));
    auto ch = epx::cosh(x);
    auto sh = epx::sinh(x);
    auto diff = epx::add(epx::mul(ch, ch), epx::opp(epx::mul(sh, sh)));
    EXPECT_EQ("1.0000000000", epx::to_string(diff, 10));
  }
  // pow(2, 10) = 1024
  {
    auto two = epx::make_q(stolz("2"), stolz("1"));
    auto ten = epx::make_q(stolz("10"), stolz("1"));
    EXPECT_EQ("1024", epx::to_string(epx::pow(two, ten), 0));
  }
  // log_10(100) = 2
  {
    auto hundred = epx::make_q(stolz("100"), stolz("1"));
    auto ten = epx::make_q(stolz("10"), stolz("1"));
    EXPECT_EQ("2.0000000000", epx::to_string(epx::log_base(hundred, ten), 10));
  }
}

}  // namespace epxut
