// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// std
#include <limits>

// epx
#include "z.hpp"

// ut
#include "def.hpp"

namespace epxut {

TEST(z_tests, create) {
  {
    auto num = epx::create<sz::container_type>(0);
    EXPECT_TRUE(epx::is_zero(num));
  }
  {
    auto num = epx::create<sz::container_type>(0x030201);
    sz expected{.digits = {1, 2, 3}};
    EXPECT_EQ(expected, num);
  }
  {
    auto num = epx::create<mz::container_type>(-2);
    mz expected{.digits = {2}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected, num);
  }
}

TEST(z_tests, normalize) {
  {
    sz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));
    epx::normalize(zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));
  }
  {
    mz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));

    zero.digits = {0, 0};
    EXPECT_FALSE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));
    epx::normalize(zero);
    epx::add_n(zero, zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));
  }
  {
    lz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));

    zero.digits = {0, 0};
    zero.sgn = epx::sign::negative;
    EXPECT_FALSE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_negative(zero));
    epx::normalize(zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_negative(zero));
  }
  {
    sz num = {.digits = {0, 1, 2, 0, 0}, .sgn = epx::sign::negative};
    EXPECT_FALSE(epx::is_zero(num));
    EXPECT_TRUE(epx::is_negative(num));
    epx::normalize(num);
    auto expected_digits = sz::container_type{0, 1, 2};
    EXPECT_EQ(expected_digits, num.digits);
    EXPECT_TRUE(epx::is_negative(num));
  }
}

TEST(z_tests, negate) {
  {
    sz zero;
    epx::negate(zero);
    EXPECT_TRUE(epx::is_zero(zero));
  }
  {
    sz num = {.digits = {1}};
    sz expected = {.digits = {1}, .sgn = epx::sign::negative};
    epx::negate(num);
    EXPECT_EQ(expected, num);
  }
  {
    sz num = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected = {.digits = {1}};
    epx::negate(num);
    EXPECT_EQ(expected, num);
  }
}

TEST(z_tests, add) {
  {
    sz zero, one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative};
    EXPECT_TRUE(epx::is_zero(epx::add(zero, zero)));
    EXPECT_TRUE(epx::is_zero(epx::add(one, minus_one)));
    EXPECT_EQ(one, epx::add(zero, one));
    EXPECT_EQ(one, epx::add(one, zero));
  }
  {
    sz a = {.digits = {255}};
    sz b = {.digits = {1}};
    sz expected = {.digits = {0, 1}};
    EXPECT_EQ(expected, epx::add(a, b));
    EXPECT_EQ(expected, epx::add(b, a));
  }
  {
    sz a = {.digits = {255}, .sgn = epx::sign::negative};
    sz b = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected = {.digits = {0, 1}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected, epx::add(a, b));
    EXPECT_EQ(expected, epx::add(b, a));
  }
}

TEST(z_tests, sub) {
  sz zero, one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative}, two{.digits = {2}},
      minus_two{.digits = {2}, .sgn = epx::sign::negative};
  EXPECT_TRUE(epx::is_zero(epx::sub(zero, zero)));
  EXPECT_TRUE(epx::is_zero(epx::sub(one, one)));
  EXPECT_TRUE(epx::is_zero(epx::sub(minus_one, minus_one)));
  EXPECT_EQ(two, epx::sub(one, minus_one));
  EXPECT_EQ(minus_two, epx::sub(minus_one, one));
}

TEST(z_tests, mul) {
  {
    sz zero, one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative};
    auto neg_one = sz{.digits = {1}, .sgn = epx::sign::negative};

    EXPECT_TRUE(epx::is_zero(epx::mul(zero, zero)));
    EXPECT_TRUE(epx::is_zero(epx::mul(zero, one)));
    EXPECT_TRUE(epx::is_zero(epx::mul(one, zero)));
    EXPECT_FALSE(epx::is_negative(epx::mul(minus_one, zero)));

    EXPECT_EQ(one, epx::mul(one, one));
    EXPECT_EQ(one, epx::mul(minus_one, minus_one));
    EXPECT_EQ(neg_one, epx::mul(one, minus_one));
    EXPECT_EQ(neg_one, epx::mul(minus_one, one));
  }
  {
    sz a = {.digits = {255}};
    sz b = {.digits = {2}};
    sz expected = {.digits = {254, 1}};
    EXPECT_EQ(expected, epx::mul(a, b));

    a.sgn = epx::sign::negative;
    auto neg_expected = expected;
    neg_expected.sgn = epx::sign::negative;
    EXPECT_EQ(neg_expected, epx::mul(a, b));
    EXPECT_EQ(neg_expected, epx::mul(b, a));
    b.sgn = epx::sign::negative;
    EXPECT_EQ(expected, epx::mul(a, b));
  }
  {
    sz a = {.digits = {1, 1}};            // 257
    sz b = {.digits = {2, 1}};            // 258
    sz expected = {.digits = {2, 3, 1}};  // 66'306
    EXPECT_EQ(expected, epx::mul(a, b));
    EXPECT_EQ(expected, epx::mul(b, a));

    a.sgn = epx::sign::negative;
    auto neg_expected = expected;
    neg_expected.sgn = epx::sign::negative;
    EXPECT_EQ(neg_expected, epx::mul(a, b));
    EXPECT_EQ(neg_expected, epx::mul(b, a));

    b.sgn = epx::sign::negative;
    EXPECT_EQ(expected, epx::mul(a, b));
  }
  {
    using dz = epx::z<>;
    using D = typename dz::digit_type;
    dz a = {.digits = {std::numeric_limits<D>::max()}};
    dz b = {.digits = {static_cast<D>(2)}};
    dz expected = {.digits = {static_cast<D>(std::numeric_limits<D>::max() - 1), static_cast<D>(1)}};
    EXPECT_EQ(expected, epx::mul(a, b));

    a.sgn = epx::sign::negative;
    auto neg_expected = expected;
    neg_expected.sgn = epx::sign::negative;
    EXPECT_EQ(neg_expected, epx::mul(a, b));
    EXPECT_EQ(neg_expected, epx::mul(b, a));

    b.sgn = epx::sign::negative;
    EXPECT_EQ(expected, epx::mul(a, b));
  }
}

TEST(z_tests, div) {
  {
    sz u = {.digits = {11}}, v = {.digits = {7}};  // 11 / 7 = 1 ... 4
    auto [q, r] = epx::div(u, v);
    EXPECT_EQ(sz{.digits = {1}}, q);
    EXPECT_EQ(sz{.digits = {4}}, r);
  }
  {
    sz u = {.digits = {11}}, v = {.digits = {7}, .sgn = epx::sign::negative};  // 11/ -7 = -1 ... 4
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {4}};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative}, v = {.digits = {7}};  //-11 / 7 = -1 ... -4
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {4}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative},
       v = {.digits = {7}, .sgn = epx::sign::negative};  // -11/-7 = 1 ... -4
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}};
    sz expected_r = {.digits = {4}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
}

TEST(z_tests, floor_div) {
  {
    sz u = {.digits = {11}}, v = {.digits = {7}};  // 11 / 7 = 1 ... 4
    auto [q, r] = epx::floor_div(u, v);
    EXPECT_EQ(sz{.digits = {1}}, q);
    EXPECT_EQ(sz{.digits = {4}}, r);
  }
  {
    sz u = {.digits = {11}}, v = {.digits = {7}, .sgn = epx::sign::negative};  // 11/ -7 = -2 ... -3
    auto [q, r] = epx::floor_div(u, v);
    sz expected_q = {.digits = {2}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {3}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative}, v = {.digits = {7}};  //-11 / 7 = -2 ... 3
    auto [q, r] = epx::floor_div(u, v);
    sz expected_q = {.digits = {2}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {3}};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative},
       v = {.digits = {7}, .sgn = epx::sign::negative};  // -11/-7 = 1 ... -4
    auto [q, r] = epx::floor_div(u, v);
    sz expected_q = {.digits = {1}};
    sz expected_r = {.digits = {4}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
}

TEST(z_tests, ceil_div) {
  {
    sz u = {.digits = {11}}, v = {.digits = {7}};  // 11 / 7 = 2 ... -3
    auto [q, r] = epx::ceil_div(u, v);
    sz expected_q = {.digits = {2}};
    sz expected_r = {.digits = {3}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}}, v = {.digits = {7}, .sgn = epx::sign::negative};  // 11/ -7 = -1 ... 4
    auto [q, r] = epx::ceil_div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {4}};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative}, v = {.digits = {7}};  //-11 / 7 = -1 ... -4
    auto [q, r] = epx::ceil_div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {4}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative},
       v = {.digits = {7}, .sgn = epx::sign::negative};  // -11/-7 = 2 ... 3
    auto [q, r] = epx::ceil_div(u, v);
    sz expected_q = {.digits = {2}};
    sz expected_r = {.digits = {3}};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
}

TEST(z_tests, mul_4exp) {
  {
    sz num{.digits = {1, 2, 3}};
    sz expected = num;
    epx::mul_4exp(num, 0);
    EXPECT_EQ(expected, num);
  }
  {
    sz num{.digits = {1}};       // 1
    sz expected{.digits = {4}};  // 1 << 2 = 4
    epx::mul_4exp(num, 1);
    EXPECT_EQ(expected, num);
  }
  {
    sz num{.digits = {1}};          // 1
    sz expected{.digits = {0, 1}};  // 1 << 8 = 256
    epx::mul_4exp(num, 4);
    EXPECT_EQ(expected, num);
  }
  {
    sz num{.digits = {4, 1}};       // 260
    sz expected{.digits = {0x41}};  // 260 >> 2 = 65
    epx::mul_4exp(num, -1);
    EXPECT_EQ(expected, num);
  }
  {
    sz num{.digits = {1, 2, 3}, .sgn = epx::sign::negative};
    epx::mul_4exp(num, -32);  // Large negative exponent
    sz expected{};
    EXPECT_EQ(expected, num);
    EXPECT_FALSE(epx::is_negative(num));  // Should reset sign to positive
  }
  {
    sz num{};
    epx::mul_4exp(num, 5);
    EXPECT_TRUE(epx::is_zero(num));
    epx::mul_4exp(num, -5);
    EXPECT_TRUE(epx::is_zero(num));
  }
}

TEST(z_tests, pow) {
  {
    sz zero;
    sz one = {.digits = {1}};
    auto pwr = epx::pow(zero, 0);  // 0^0 = 1
    EXPECT_EQ(one, pwr);
  }
  {
    sz zero;
    auto pwr = epx::pow(zero, 1);  // 0^1 = 0
    EXPECT_EQ(zero, pwr);
  }
  {
    auto num = create_sz(2);
    auto expected = create_sz(64);
    auto pwr = epx::pow(num, 6);  //  2^6 = 64
    EXPECT_EQ(expected, pwr);
  }
  {
    sz one = {.digits = {1}};
    sz minus_one{.digits = {1}, .sgn = epx::sign::negative};
    EXPECT_EQ(minus_one, epx::pow(minus_one, 9));  // (-1)^9 = -1
    EXPECT_EQ(one, epx::pow(minus_one, 10));       // (-1)^10 = 1
  }
}

TEST(z_tests, root) {
  {
    sz zero;
    EXPECT_TRUE(epx::is_zero(epx::root(zero, 1)));
    EXPECT_TRUE(epx::is_zero(epx::root(zero, 2)));
    EXPECT_TRUE(epx::is_zero(epx::root(zero, 7)));
  }
  {
    sz one = {.digits = {1}};
    EXPECT_EQ(one, epx::root(sz{}, 0));
    EXPECT_EQ(one, epx::root(one, 0));
    EXPECT_EQ(one, epx::root(create_sz(2), 0));
    EXPECT_EQ(one, epx::root(create_sz(9), 0));
  }
  {
    EXPECT_EQ(create_sz(2), epx::root(create_sz(4), 2));
    EXPECT_EQ(create_sz(3), epx::root(create_sz(9), 2));
    EXPECT_EQ(create_sz(19), epx::root(create_sz(361), 2));
    EXPECT_EQ(create_sz(1000), epx::root(create_sz(1000000), 2));
  }
  {
    EXPECT_EQ(create_sz(5), epx::root(create_sz(125), 3));
    EXPECT_EQ(create_sz(37398), epx::root(create_sz(1956111062177043216), 4));
  }
  {
    EXPECT_EQ(create_sz(111), epx::root(create_sz(12345), 2));
    EXPECT_EQ(create_sz(23), epx::root(create_sz(12345), 3));
  }
}

// Generated by AI (Claude Opus 4.6) — Additional tests for improving code coverage

TEST(z_tests, root_negative_radicand) {
  // z::root should throw on negative input
  sz neg = {.digits = {5}, .sgn = epx::sign::negative};
  EXPECT_THROW(epx::root(neg, 2), epx::negative_radicand_error);
  EXPECT_THROW(epx::root(neg, 3), epx::negative_radicand_error);
}

TEST(z_tests, mz_add) {
  // Test add/add_n with mz type (uint16_t digits) for template branch coverage
  {
    mz zero, one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative};
    EXPECT_TRUE(epx::is_zero(epx::add(zero, zero)));
    EXPECT_TRUE(epx::is_zero(epx::add(one, minus_one)));
    EXPECT_EQ(one, epx::add(zero, one));
    EXPECT_EQ(one, epx::add(one, zero));
  }
  {
    // Carry propagation: 65535 + 1 = 65536 = {0, 1}
    mz a = {.digits = {65535}};
    mz b = {.digits = {1}};
    mz expected = {.digits = {0, 1}};
    EXPECT_EQ(expected, epx::add_n(a, b));
    EXPECT_EQ(expected, epx::add_n(b, a));
  }
  {
    // Different signs, |lhs| < |rhs| — covers swap branch in signed add()
    mz small = {.digits = {3}};
    mz big = {.digits = {10}, .sgn = epx::sign::negative};
    mz expected = {.digits = {7}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected, epx::add(small, big));
  }
}

TEST(z_tests, mz_sub) {
  // Test sub/sub_n with mz type for template branch coverage
  {
    mz a = {.digits = {0, 1}};  // 65536
    mz b = {.digits = {1}};
    mz expected = {.digits = {65535}};  // 65535
    EXPECT_EQ(expected, epx::sub_n(a, b));
  }
  {
    mz one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative};
    mz two{.digits = {2}}, minus_two{.digits = {2}, .sgn = epx::sign::negative};
    EXPECT_EQ(two, epx::sub(one, minus_one));
    EXPECT_EQ(minus_two, epx::sub(minus_one, one));
  }
}

TEST(z_tests, mz_mul) {
  // Test mul with mz type for template branch coverage
  {
    mz a = {.digits = {1000}}, b = {.digits = {2000}};
    auto result = epx::mul(a, b);
    // 1000 * 2000 = 2,000,000 = 0x1E8480. In 16-bit words: low=0x8480, high=0x1E
    mz expected = {.digits = {0x8480, 0x001E}};
    EXPECT_EQ(expected, result);
  }
  {
    mz zero, one{.digits = {1}}, minus_one{.digits = {1}, .sgn = epx::sign::negative};
    EXPECT_TRUE(epx::is_zero(epx::mul(zero, one)));
    EXPECT_EQ(one, epx::mul(one, one));
    EXPECT_EQ(one, epx::mul(minus_one, minus_one));
  }
}

TEST(z_tests, mz_div) {
  // Test div/floor_div/ceil_div with mz type for template branch coverage
  {
    mz u = {.digits = {100}}, v = {.digits = {7}};  // 100 / 7 = 14 r 2
    auto [q, r] = epx::div(u, v);
    EXPECT_EQ(mz{.digits = {14}}, q);
    EXPECT_EQ(mz{.digits = {2}}, r);
  }
  {
    mz u = {.digits = {100}, .sgn = epx::sign::negative}, v = {.digits = {7}};  // -100 / 7
    auto [q, r] = epx::floor_div(u, v);
    mz expected_q = {.digits = {15}, .sgn = epx::sign::negative};  // -15
    mz expected_r = {.digits = {5}};                               // 5 (-100 = -15*7 + 5)
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    mz u = {.digits = {100}}, v = {.digits = {7}};  // ceil(100/7) = 15
    auto [q, r] = epx::ceil_div(u, v);
    mz expected_q = {.digits = {15}};
    mz expected_r = {.digits = {5}, .sgn = epx::sign::negative};  // -5 (100 = 15*7 - 5)
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
}

TEST(z_tests, mz_mul_4exp) {
  // Test mul_4exp with mz type for template branch coverage
  {
    mz num{.digits = {1}};
    mz expected{.digits = {16}};  // 1 * 4^2 = 16
    epx::mul_4exp(num, 2);
    EXPECT_EQ(expected, num);
  }
  {
    mz num{.digits = {64}};
    mz expected{.digits = {4}};  // 64 / 4^2 = 4
    epx::mul_4exp(num, -2);
    EXPECT_EQ(expected, num);
  }
}

TEST(z_tests, mz_pow) {
  // Test pow with mz type for template branch coverage
  {
    mz base = {.digits = {3}};
    auto result = epx::pow(base, 4);  // 3^4 = 81
    mz expected = {.digits = {81}};
    EXPECT_EQ(expected, result);
  }
  {
    mz base = {.digits = {5}};
    auto result = epx::pow(base, 0);  // 5^0 = 1
    mz expected = {.digits = {1}};
    EXPECT_EQ(expected, result);
  }
}

TEST(z_tests, mz_root) {
  // Generated by AI (Claude Opus 4.6) — Tests z::root with mz type including
  // Newton iteration. This also covers the bit_length() fix (sizeof(D)*CHAR_BIT).
  // All expected values cross-checked with Node.js BigInt arithmetic.
  auto stomz = [](std::string_view s) { return epx::try_from_chars<mz::container_type>(s).value(); };

  // --- Trivial cases: k=0, k=1, zero ---
  {
    mz one = {.digits = {1}};
    EXPECT_EQ(one, epx::root(mz{}, 0));
    EXPECT_EQ(one, epx::root(one, 0));
  }
  {
    mz num = {.digits = {42}};
    EXPECT_EQ(num, epx::root(num, 1));
  }
  {
    mz zero;
    EXPECT_TRUE(epx::is_zero(epx::root(zero, 2)));
    EXPECT_TRUE(epx::is_zero(epx::root(zero, 5)));
  }

  // --- Small perfect squares ---
  {
    auto num = epx::create<mz::container_type>(144);
    auto expected = epx::create<mz::container_type>(12);
    EXPECT_EQ(expected, epx::root(num, 2));
  }
  {
    auto num = epx::create<mz::container_type>(10000);
    auto expected = epx::create<mz::container_type>(100);
    EXPECT_EQ(expected, epx::root(num, 2));
  }
  // --- Small perfect cube ---
  {
    auto num = epx::create<mz::container_type>(125);
    auto expected = epx::create<mz::container_type>(5);
    EXPECT_EQ(expected, epx::root(num, 3));
  }

  // --- Large perfect squares (multi-digit mz) ---
  EXPECT_EQ(stomz("12345"), epx::root(stomz("152399025"), 2));       // sqrt(12345^2)
  EXPECT_EQ(stomz("1000000"), epx::root(stomz("1000000000000"), 2)); // sqrt(10^12)
  EXPECT_EQ(stomz("9999"), epx::root(stomz("99980001"), 2));         // sqrt(9999^2)

  // --- Large perfect cubes ---
  EXPECT_EQ(stomz("1000"), epx::root(stomz("1000000000"), 3));       // cbrt(10^9)
  EXPECT_EQ(stomz("20000"), epx::root(stomz("8000000000000"), 3));   // cbrt(20000^3)
  EXPECT_EQ(stomz("123"), epx::root(stomz("1860867"), 3));           // cbrt(123^3)
  EXPECT_EQ(stomz("12345"), epx::root(stomz("1881365963625"), 3));   // cbrt(12345^3)

  // --- 4th roots ---
  EXPECT_EQ(stomz("100"), epx::root(stomz("100000000"), 4));         // 4th-root(10^8)
  EXPECT_EQ(stomz("6"), epx::root(stomz("1296"), 4));                // 4th-root(6^4)
  EXPECT_EQ(stomz("999"), epx::root(stomz("996005996001"), 4));      // 4th-root(999^4)
  EXPECT_EQ(stomz("1024"), epx::root(stomz("1099511627776"), 4));    // 4th-root(2^40)

  // --- 5th roots ---
  EXPECT_EQ(stomz("20"), epx::root(stomz("3200000"), 5));            // 5th-root(20^5)
  EXPECT_EQ(stomz("64"), epx::root(stomz("1073741824"), 5));         // 5th-root(2^30)
  EXPECT_EQ(stomz("1024"), epx::root(stomz("1125899906842624"), 5)); // 5th-root(2^50)

  // --- 7th root ---
  EXPECT_EQ(stomz("50"), epx::root(stomz("781250000000"), 7));       // 7th-root(50^7)

  // --- Non-perfect roots (floor) ---
  EXPECT_EQ(stomz("11111"), epx::root(stomz("123456789"), 2));       // floor(sqrt(123456789))
  EXPECT_EQ(stomz("995"), epx::root(stomz("987654321"), 3));         // floor(cbrt(987654321))
  EXPECT_EQ(stomz("999999"), epx::root(stomz("999999999999"), 2));   // floor(sqrt(~10^12))
  EXPECT_EQ(stomz("9999"), epx::root(stomz("999999999999"), 3));     // floor(cbrt(~10^12))
  EXPECT_EQ(stomz("999"), epx::root(stomz("999999999999"), 4));      // floor(4th-root(~10^12))

  // --- Very large numbers (>64-bit) ---
  EXPECT_EQ(stomz("4294967296"), epx::root(stomz("18446744073709551616"), 2)); // sqrt(2^64) = 2^32
  EXPECT_EQ(stomz("9999999999"), epx::root(stomz("99999999999999999999"), 2)); // floor(sqrt(~10^20))
  EXPECT_EQ(stomz("4641588"), epx::root(stomz("99999999999999999999"), 3));    // floor(cbrt(~10^20))
}

}  // namespace epxut
