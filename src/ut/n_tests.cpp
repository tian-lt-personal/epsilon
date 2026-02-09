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
      sz zero, one{.digits = {1}};
      epx::div_n(zero, one);
    }
  }
}

}  // namespace epxut
