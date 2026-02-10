// gtest
#include <gtest/gtest.h>

// epx
#include <z.hpp>

// ut
#include "def.hpp"

namespace epxut {

TEST(n_tests, cmp_n) {
  {
    sz zero, one{.digits = {1}};
    EXPECT_GT(epx::cmp_n(one, zero), 0);
    EXPECT_LT(epx::cmp_n(zero, one), 0);
    EXPECT_EQ(0, epx::cmp_n(zero, zero));
    EXPECT_EQ(0, epx::cmp_n(one, one));
  }
  {
    sz a = {.digits = {1, 2}};
    sz b = {.digits = {2, 1}};
    EXPECT_GT(epx::cmp_n(a, b), 0);
    EXPECT_LT(epx::cmp_n(b, a), 0);
  }
  {
    sz a = {.digits = {0, 0, 1}};
    sz b = {.digits = {0, 1}};
    EXPECT_GT(epx::cmp_n(a, b), 0);
    EXPECT_LT(epx::cmp_n(b, a), 0);
  }
}

TEST(n_tests, add_n) {
  {
    sz zero;
    auto r = epx::add_n(zero, zero);
    EXPECT_TRUE(epx::is_zero(zero));
  }
  {
    sz zero, one{.digits = {1}};
    auto r = epx::add_n(zero, one);
    EXPECT_EQ(one, r);
  }
  {
    sz zero, one{.digits = {1}};
    auto r = epx::add_n(one, zero);
    EXPECT_EQ(one, r);
  }
  {
    sz a = {.digits = {255}};
    sz b = {.digits = {1}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {0, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {1}};
    sz b = {.digits = {255}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {0, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {255}};
    sz b = {.digits = {255}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {254, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {0, 1}};
    sz b = {.digits = {255}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {255, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {255, 255, 255}};
    sz b = {.digits = {1}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {0, 0, 0, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {1}};
    sz b = {.digits = {255, 255, 255}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {0, 0, 0, 1}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {12}};
    sz b = {.digits = {34}};
    auto r = epx::add_n(a, b);
    sz expected = {.digits = {46}};
    EXPECT_EQ(expected, r);
  }
}

TEST(n_tests, sub_n) {
  {
    sz zero;
    EXPECT_TRUE(epx::is_zero(epx::sub_n(zero, zero)));
  }
  {
    sz zero, one{.digits = {1}};
    auto r = epx::sub_n(one, zero);
    EXPECT_EQ(one, r);
  }
  {
    sz a = {.digits = {0, 1}};
    sz b = {.digits = {1}};
    auto r = epx::sub_n(a, b);
    sz expected = {.digits = {255}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {0, 0, 1}};
    sz b = {.digits = {1}};
    auto r = epx::sub_n(a, b);
    sz expected = {.digits = {255, 255}};
    EXPECT_EQ(expected, r);
  }
  {
    sz a = {.digits = {46}};
    sz b = {.digits = {12}};
    auto r = epx::sub_n(a, b);
    sz expected = {.digits = {34}};
    EXPECT_EQ(expected, r);
  }
}

TEST(n_tests, mul_n) {
  {
    sz zero;
    auto r = epx::mul_n(zero, zero);
    EXPECT_TRUE(epx::is_zero(r));
  }
  {
    sz zero, one{.digits = {1}};
    auto r1 = epx::mul_n(zero, one);
    auto r2 = epx::mul_n(one, zero);
    EXPECT_TRUE(epx::is_zero(r1));
    EXPECT_TRUE(epx::is_zero(r2));
  }
  {
    sz a = {.digits = {7}};
    sz b = {.digits = {3}};
    sz expected = {.digits = {21}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
    EXPECT_EQ(expected, epx::mul_n(b, a));
  }
  {
    sz a = {.digits = {255}};
    sz b = {.digits = {2}};
    sz expected = {.digits = {254, 1}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
  }
  {
    sz a = {.digits = {0, 1}};
    sz b = {.digits = {0, 1}};
    sz expected = {.digits = {0, 0, 1}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
  }
  {
    sz a = {.digits = {1, 1}};
    sz b = {.digits = {1, 0, 1}};
    sz expected = {.digits = {1, 1, 1, 1}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
    EXPECT_EQ(expected, epx::mul_n(b, a));
  }
  {
    sz a = {.digits = {255, 255}};
    sz b = {.digits = {255, 255}};
    sz expected = {.digits = {1, 0, 254, 255}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
  }
  {
    sz a = {.digits = {0, 0}};
    sz b = {.digits = {5}};
    EXPECT_TRUE(epx::is_zero(epx::mul_n(a, b)));
  }
  {
    mz a = {.digits = {65535}};
    mz b = {.digits = {65535}};
    mz expected = {.digits = {1, 65534}};
    EXPECT_EQ(expected, epx::mul_n(a, b));
  }
}

TEST(n_tests, div_n) {
  {
    {
      sz zero, one = {.digits = {1}};
      auto res = epx::div_n(zero, one);
      EXPECT_TRUE(epx::is_zero(res.q));
      EXPECT_TRUE(epx::is_zero(res.r));
      EXPECT_THROW(epx::div_n(one, zero), epx::divide_by_zero_error);
    }
    {
      sz zero;
      sz a = {.digits = {0, 1}};
      sz b = {.digits = {2}};
      auto res = epx::div_n(a, b);  // 256 / 2 = 128 ... 0
      EXPECT_EQ(sz{.digits = {0x80}}, res.q);
      EXPECT_EQ(zero, res.r);
    }
    {
      sz zero;
      sz a = {.digits = {0, 0, 1}};
      sz b = {.digits = {0, 2}};
      auto res = epx::div_n(a, b);  // 65536 / 512 = 128 ... 0
      EXPECT_EQ(sz{.digits = {0x80}}, res.q);
      EXPECT_EQ(zero, res.r);
    }
    {
      sz a = {.digits = {5, 0, 0x10}};
      sz b = {.digits = {0x3f, 3}};
      auto res = epx::div_n(a, b);  // 1048581 / 831 = 1261 ... 690
      sz expected_q{.digits = {0xed, 0x4}};
      sz expected_r{.digits = {0xb2, 0x2}};
      EXPECT_EQ(expected_q, res.q);
      EXPECT_EQ(expected_r, res.r);
    }
    {
      sz a = {.digits = {5, 0, 0, 1}};
      sz b = {.digits = {0x3f, 3}};
      auto res = epx::div_n(a, b);  // 16777221 / 831 = 20189 ... 162
      sz expected_q{.digits = {0xdd, 0x4e}};
      EXPECT_EQ(expected_q, res.q);
      EXPECT_EQ(sz{.digits = {0xa2}}, res.r);
    }
    {
      sz one = {.digits = {1}};
      sz a = {.digits = {0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xff}};
      sz b = {.digits = {0x11, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xff}};
      auto res = epx::div_n(a, b);  // 18437381296131290130 / 18437381296131290129 = 1 ... 1
      EXPECT_EQ(one, res.q);
      EXPECT_EQ(one, res.r);
    }
    {
      sz zero, one = {.digits = {1}};
      sz a = {.digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12}};
      sz b = {.digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12}};
      auto res = epx::div_n(a, b);  // 1311768467463790335 / 1311768467463790335= 1 ... 0
      EXPECT_EQ(one, res.q);
      EXPECT_EQ(zero, res.r);
    }
    {
      sz one = {.digits = {1}};
      sz a = {.digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12}};
      sz b = {.digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x11}};
      auto res = epx::div_n(a, b);  // 1311768467463790335 / 1239710873425862399 = 1 ... 72057594037927936
      sz expected_r{.digits = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 1}};
      EXPECT_EQ(one, res.q);
      EXPECT_EQ(expected_r, res.r);
    }
    {
      sz a = {
          .digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12, 0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12}};
      sz b = {.digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x11}};
      auto res = epx::div_n(a, b);  // 24197857203266735141494478776148221695 / 1239710873425862399 =
                                    // 19518952137926717496 ... 1080293969106388791
      sz expected_q{.digits = {0x38, 0x30, 0xf9, 0xb8, 0x8d, 0x3f, 0xe1, 0xe, 1}};
      sz expected_r{.digits = {0x37, 0x47, 0x95, 0x29, 0xa4, 0xf9, 0xfd, 0xe}};
      EXPECT_EQ(expected_q, res.q);
      EXPECT_EQ(expected_r, res.r);
    }
    {
      sz a = {
          .digits = {0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12, 0xff, 0xde, 0xbc, 0x9a, 0x78, 0x56, 0x34, 0x12}};
      auto res = epx::div_n(a, sz{.digits = {0x25}});  // 24197857203266735141494478776148221695 / 37 =
                                                       // 653996140628830679499850777733735721 ... 18
      sz expected_q{
          .digits = {0x29, 0xa5, 0x81, 0xe8, 0xcb, 0x85, 0x4d, 0xbb, 0x14, 0xac, 0xd4, 0xcc, 0x7f, 0xf4, 0x7d}};
      sz expected_r{.digits = {0x12}};
      EXPECT_EQ(expected_q, res.q);
      EXPECT_EQ(expected_r, res.r);
    }
  }
}

}  // namespace epxut
