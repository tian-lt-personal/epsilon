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
    z<C> ten = {.digits = {10}};
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
    return res;
  } else {
    static_assert("not implemented.");
  }
}

}  // namespace epx

#endif  // EPSILON_INC_CHARS_HPP
