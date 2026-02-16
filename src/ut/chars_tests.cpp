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
    EXPECT_TRUE(epx::is_positive(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("0");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_TRUE(epx::is_positive(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-0");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_TRUE(epx::is_positive(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("000000");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_TRUE(epx::is_positive(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("1");
    sz expected = {.digits = {1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-1");
    sz expected = {.digits = {1}, .sgn = epx::sign::negative};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("000001");
    sz expected = {.digits = {1}};
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_EQ(expected, *opt_num);
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-000001");
    sz expected = {.digits = {1}, .sgn = epx::sign::negative};
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
  {
    auto opt_num = epx::try_from_chars<mz::container_type>("-4294967296");
    mz expected = {.digits = {0, 0, 1}, .sgn = epx::sign::negative};
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
    auto opt_num = epx::try_from_chars<sz::container_type>("12.34");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("+");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("++1");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("+a");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("--2");
    EXPECT_FALSE(opt_num.has_value());
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-b");
    EXPECT_FALSE(opt_num.has_value());
  }
}

TEST(chars_tests, to_decimal_string) {
  {
    auto num = epx::try_from_chars<sz::container_type>("").value();
    EXPECT_TRUE(epx::is_zero(num));
    auto str = epx::to_string(num);
    EXPECT_EQ("0", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("0").value();
    EXPECT_TRUE(epx::is_zero(num));
    auto str = epx::to_string(num);
    EXPECT_EQ("0", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("-0").value();
    EXPECT_TRUE(epx::is_zero(num));
    auto str = epx::to_string(num);
    EXPECT_EQ("0", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("0000").value();
    EXPECT_TRUE(epx::is_zero(num));
    auto str = epx::to_string(num);
    EXPECT_EQ("0", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("1").value();
    auto str = epx::to_string(num);
    EXPECT_EQ("1", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("-1").value();
    auto str = epx::to_string(num);
    EXPECT_EQ("-1", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("-0000100").value();
    auto str = epx::to_string(num);
    EXPECT_EQ("-100", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("12345678910").value();
    auto str = epx::to_string(num);
    EXPECT_EQ("12345678910", str);
  }
  {
    auto num = epx::try_from_chars<sz::container_type>("-12345678910").value();
    auto str = epx::to_string(num);
    EXPECT_EQ("-12345678910", str);
  }
  {
    auto num = epx::try_from_chars<mz::container_type>(
                   "000009837417236498716239487169238746192783461928736491287364192873461928736416341997812934781692387"
                   "46192374")
                   .value();
    auto str = epx::to_string(num);
    EXPECT_EQ("983741723649871623948716923874619278346192873649128736419287346192873641634199781293478169238746192374",
              str);
  }
}

}  // namespace epxut
