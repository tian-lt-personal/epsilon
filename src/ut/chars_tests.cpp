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

TEST(chars_tests, z_to_decimal_string) {
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

TEST(chars_tests, r_to_decimal_string) {
  {
    auto q = epx::make_q(stosz("0"), stosz("1"));
    EXPECT_EQ("0", epx::to_string(q, 0));
    EXPECT_EQ("0.0", epx::to_string(q, 1));
  }
  {
    auto q = epx::make_q(stosz("1"), stosz("1"));
    EXPECT_EQ("1", epx::to_string(q, 0));
    EXPECT_EQ("1.0", epx::to_string(q, 1));
    EXPECT_EQ("1.00", epx::to_string(q, 2));
    EXPECT_EQ("1.000000000000000", epx::to_string(q, 15));
  }
  {
    auto q = epx::make_q(stosz("2"), stosz("1"));
    EXPECT_EQ("2", epx::to_string(q, 0));
    EXPECT_EQ("2.0", epx::to_string(q, 1));
    EXPECT_EQ("2.00000", epx::to_string(q, 5));
  }
  {
    auto q = epx::make_q(stosz("1"), stosz("2"));
    auto s = epx::to_string(q, 0);
    EXPECT_EQ("1", epx::to_string(q, 0));
    EXPECT_EQ("0.5", epx::to_string(q, 1));
    EXPECT_EQ("0.5000", epx::to_string(q, 4));
  }
  {
    auto q = epx::make_q(stosz("1"), stosz("3"));
    EXPECT_EQ("0", epx::to_string(q, 0));
    EXPECT_EQ("0.3", epx::to_string(q, 1));
    EXPECT_EQ("0.333333", epx::to_string(q, 6));
  }
  {
    auto q = epx::make_q(stosz("4"), stosz("5"));
    EXPECT_EQ("1", epx::to_string(q, 0));
    EXPECT_EQ("0.8", epx::to_string(q, 1));
    EXPECT_EQ("0.8000", epx::to_string(q, 4));
  }
  {
    auto q = epx::make_q(stosz("-1"), stosz("1"));
    EXPECT_EQ("-1", epx::to_string(q, 0));
    EXPECT_EQ("-1.0", epx::to_string(q, 1));
  }
  {
    auto q = epx::make_q(stosz("4"), stosz("-5"));
    EXPECT_EQ("-1", epx::to_string(q, 0));
    EXPECT_EQ("-0.8", epx::to_string(q, 1));
    EXPECT_EQ("-0.8000", epx::to_string(q, 4));
  }
  {
    auto q = epx::make_q(stosz("832798712356987419287340"), stosz("819725918374510348"));
    EXPECT_EQ("1015948", epx::to_string(q, 0));
    EXPECT_EQ("1015947.8", epx::to_string(q, 1));
    EXPECT_EQ("1015947.76", epx::to_string(q, 2));
    EXPECT_EQ("1015947.76", epx::to_string(q, 2));
    EXPECT_EQ("1015947.76215982044", epx::to_string(q, 11));
    EXPECT_EQ("1015947.7621598204356532794380011", epx::to_string(q, 25));
  }
}

}  // namespace epxut
