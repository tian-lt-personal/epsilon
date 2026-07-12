// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <fstream>
#include <string>

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
    EXPECT_FALSE(epx::is_negative(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("0");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_FALSE(epx::is_negative(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("-0");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_FALSE(epx::is_negative(*opt_num));
  }
  {
    auto opt_num = epx::try_from_chars<sz::container_type>("000000");
    EXPECT_TRUE(opt_num.has_value());
    EXPECT_TRUE(epx::is_zero(*opt_num));
    EXPECT_FALSE(epx::is_negative(*opt_num));
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

// Ad hoc reproduction tests for double-rounding bug in to_string(r<C>, k).
// These tests exercise rounding boundary behavior that the old formula got wrong.

TEST(chars_tests, r_to_string_repro_double_rounding_2_3) {
  auto q = epx::make_q(stosz("2"), stosz("3"));
  // 2/3 = 0.666... All digits are 6, so the (k+1)th digit is always 6 >= 5.
  // The k-th digit always rounds UP from 6 to 7.
  // Pattern: k-1 sixes, then a 7.
  EXPECT_EQ("0.7", epx::to_string(q, 1));
  EXPECT_EQ("0.67", epx::to_string(q, 2));
  EXPECT_EQ("0.6666666667", epx::to_string(q, 10));
  EXPECT_EQ("0.66666666666666666667", epx::to_string(q, 20));
  EXPECT_EQ("0.6666666666666666666666666666666666666667", epx::to_string(q, 40));
  EXPECT_EQ("0.66666666666666666666666666666666666666666666666667", epx::to_string(q, 50));
}

TEST(chars_tests, r_to_string_repro_double_rounding_1_3) {
  auto q = epx::make_q(stosz("1"), stosz("3"));
  // 1/3 = 0.333... All digits are 3 < 5, so we NEVER round up.
  // Pattern: k threes.
  EXPECT_EQ("0.3", epx::to_string(q, 1));
  EXPECT_EQ("0.33", epx::to_string(q, 2));
  EXPECT_EQ("0.3333333333", epx::to_string(q, 10));
  EXPECT_EQ("0.33333333333333333333", epx::to_string(q, 20));
  EXPECT_EQ("0.3333333333333333333333333333333333333333", epx::to_string(q, 40));
  EXPECT_EQ("0.33333333333333333333333333333333333333333333333333", epx::to_string(q, 50));
}

TEST(chars_tests, r_to_string_repro_halfway) {
  // Verify correct rounding at exact halfway points.
  // 1/2 = 0.5 -> at k=0: round(0.5) = 1
  EXPECT_EQ("1", epx::to_string(epx::make_q(stosz("1"), stosz("2")), 0));
  EXPECT_EQ("0.5", epx::to_string(epx::make_q(stosz("1"), stosz("2")), 1));
  EXPECT_EQ("0.50", epx::to_string(epx::make_q(stosz("1"), stosz("2")), 2));

  // 1/4 = 0.25 -> at k=1: round(2.5) = 3 -> "0.3"
  EXPECT_EQ("0.3", epx::to_string(epx::make_q(stosz("1"), stosz("4")), 1));

  // 3/4 = 0.75 -> at k=1: round(7.5) = 8 -> "0.8"
  EXPECT_EQ("0.8", epx::to_string(epx::make_q(stosz("3"), stosz("4")), 1));

  // 1/8 = 0.125 -> at k=2: round(12.5) = 13 -> "0.13"
  EXPECT_EQ("0.13", epx::to_string(epx::make_q(stosz("1"), stosz("8")), 2));

  // 5/8 = 0.625 -> at k=2: round(62.5) = 63 -> "0.63"
  EXPECT_EQ("0.63", epx::to_string(epx::make_q(stosz("5"), stosz("8")), 2));
}

TEST(chars_tests, r_to_string_pi_exact_match) {
  // Read the hundred_k_pi.txt reference file.
  std::string path = std::string(EPSILON_UT_SRC_DIR) + "/hundred_k_pi.txt";
  std::ifstream file(path);
  ASSERT_TRUE(file.is_open()) << "Cannot open reference file: " << path;

  std::string content;
  std::getline(file, content);
  ASSERT_GE(content.size(), 3u) << "Reference file too short";
  ASSERT_EQ('3', content[0]);
  ASSERT_EQ('.', content[1]);

  // Fractional digits only (everything after "3.").
  std::string raw_digits = content.substr(2);

  // Full 100k-digit test only in release builds — the computation is too
  // expensive for debug.  In debug, a quick 500-digit exact match catches
  // regressions without slowing down the edit/compile/test cycle.
#ifdef NDEBUG
  int max_k = static_cast<int>(raw_digits.size()) - 1;  // need one extra digit for rounding
  int k = std::min(100000, max_k);
#else
  int k = 500;
#endif

  // Compute pi via Gauss's formula: pi = 4 * arctan(1).
  auto one = epx::make_q(stosz("1"), stosz("1"));
  auto four = epx::make_q(stosz("4"), stosz("1"));
  auto pi = epx::mul(four, epx::arctan(one));

  // Build the correctly-rounded expected string.
  std::string frac = raw_digits.substr(0, k);
  char next_digit = raw_digits[k];  // the (k+1)th fractional digit
  if (next_digit >= '5') {
    int i = k - 1;
    while (i >= 0 && frac[i] == '9') {
      frac[i] = '0';
      --i;
    }
    if (i >= 0) {
      frac[i]++;
    } else {
      frac.insert(frac.begin(), '1');
    }
  }
  std::string expected = "3." + frac;

  EXPECT_EQ(expected, epx::to_string(pi, k));
}

}  // namespace epxut
