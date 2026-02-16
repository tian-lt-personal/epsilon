// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_INC_CHARS_HPP
#define EPSILON_INC_CHARS_HPP

// std
#include <cctype>
#include <optional>
#include <string_view>

// epx
#include "ops.hpp"
#include "z.hpp"

namespace epx {

template <container C, int B = 10>
std::optional<z<C>> try_from_chars(std::string_view chars) {
  using D = typename z<C>::digit_type;

  if constexpr (B == 10) {
    const z<C> ten = {.digits = {10}};

    sign sgn;
    if (chars.length() > 1) {
      if (chars.front() == '+' && isdigit(static_cast<unsigned int>(chars[1]))) {
        sgn = sign::positive;
        chars = chars.substr(1);
      } else if (chars.front() == '-' && isdigit(static_cast<unsigned int>(chars[1]))) {
        sgn = sign::negative;
        chars = chars.substr(1);
      } else if (isdigit(static_cast<unsigned int>(chars.front()))) {
        sgn = sign::positive;
      } else {
        return std::nullopt;
      }
    } else {
      sgn = sign::positive;
    }

    z<C> res;
    for (auto ch : chars) {
      if (!isdigit(static_cast<unsigned int>(ch))) {
        return std::nullopt;
      }
      int digit = ch - '0';
      z<C> tmp;
      if (digit > 0) {
        tmp = {.digits = {static_cast<D>(digit)}};
      }
      res = res * ten + tmp;
    }
    if (!is_zero(res)) {
      res.sgn = sgn;
    }
    return res;
  } else {
    static_assert("not implemented.");
  }
}

template <container C, int B = 10>
std::string to_string(z<C> num) {
  if constexpr (B == 10) {
    if (is_zero(num)) {
      return "0";
    }

    std::string res;
    sign sgn = num.sgn;
    const z<C> ten = {.digits = {10}};
    while (!is_zero(num)) {
      auto [q, r] = div_n(num, ten);
      if (is_zero(r)) {
        res.push_back('0');
      } else {
        res.push_back(static_cast<char>(r.digits[0]) + '0');
      }
      num = std::move(q);
    }
    if (sgn == sign::negative) {
      res.push_back('-');
    }
    std::ranges::reverse(res);
    return res;
  } else {
    static_assert("not implemented.");
  }
}

}  // namespace epx

#endif  // EPSILON_INC_CHARS_HPP
