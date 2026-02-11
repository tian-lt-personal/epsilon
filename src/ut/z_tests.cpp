// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include <z.hpp>

// std
#include <limits>

// ut
#include "def.hpp"

namespace epxut {

TEST(z_tests, normalize) {
  {
    sz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));
    epx::normalize(zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));
  }
  {
    mz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));

    zero.digits = {0, 0};
    EXPECT_FALSE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));
    epx::normalize(zero);
    epx::add_n(zero, zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));
  }
  {
    lz zero;
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));

    zero.digits = {0, 0};
    zero.sgn = epx::sign::negative;
    EXPECT_FALSE(epx::is_zero(zero));
    EXPECT_FALSE(epx::is_positive(zero));
    epx::normalize(zero);
    EXPECT_TRUE(epx::is_zero(zero));
    EXPECT_TRUE(epx::is_positive(zero));
  }
  {
    sz num = {.digits = {0, 1, 2, 0, 0}, .sgn = epx::sign::negative};
    EXPECT_FALSE(epx::is_zero(num));
    EXPECT_FALSE(epx::is_positive(num));
    epx::normalize(num);
    auto expected_digits = sz::container_type{0, 1, 2};
    EXPECT_EQ(expected_digits, num.digits);
    EXPECT_FALSE(epx::is_positive(num));
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
    EXPECT_TRUE(epx::is_positive(epx::mul(minus_one, zero)));

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
    sz u = {.digits = {11}}, v = {.digits = {7}};
    auto [q, r] = epx::div(u, v);
    EXPECT_EQ(sz{.digits = {1}}, q);
    EXPECT_EQ(sz{.digits = {4}}, r);
  }
  {
    sz u = {.digits = {11}}, v = {.digits = {7}, .sgn = epx::sign::negative};
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {3}, .sgn = epx::sign::negative};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative}, v = {.digits = {7}};
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}, .sgn = epx::sign::negative};
    sz expected_r = {.digits = {3}};
    EXPECT_EQ(expected_q, q);
    EXPECT_EQ(expected_r, r);
  }
  {
    sz u = {.digits = {11}, .sgn = epx::sign::negative}, v = {.digits = {7}, .sgn = epx::sign::negative};
    auto [q, r] = epx::div(u, v);
    sz expected_q = {.digits = {1}};
    sz expected_r = {.digits = {4}, .sgn = epx::sign::negative};
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
    EXPECT_TRUE(epx::is_positive(num));  // Should reset sign to positive
  }
  {
    sz num{};
    epx::mul_4exp(num, 5);
    EXPECT_TRUE(epx::is_zero(num));
    epx::mul_4exp(num, -5);
    EXPECT_TRUE(epx::is_zero(num));
  }
}

}  // namespace epxut
