// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// gtest
#include <gtest/gtest.h>

// epx
#include "chars.hpp"
#include "z.hpp"

// ut
#include "def.hpp"

namespace epxut {

TEST(chars_tests, try_from_decimal_chars) {
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("0");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("1");
    sz expected = {.digits = {1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("10");
    sz expected = {.digits = {0xa}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("255");
    sz expected = {.digits = {0xff}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("256");
    sz expected = {.digits = {0, 1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("65536");
    sz expected = {.digits = {0, 0, 1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("16777216");
    sz expected = {.digits = {0, 0, 0, 1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<mz::container_type>("4294967296");
    mz expected = {.digits = {0, 0, 1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
}

TEST(chars_tests, try_from_chars_bad_chars) {
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("abc");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("12*34");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("12 34");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("+1");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-1");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("12.34");
    EXPECT_FALSE(opt_num.has_value());
  }
}

}  // namespace epxut
