// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_INC_R_HPP
#define EPSILON_INC_R_HPP

// std
#include <functional>

// epx
#include "coro.hpp"
#include "z.hpp"

namespace epx {

template <container C>
using r = std::function<coro::lazy<z<C>>(unsigned int)>;  // q (n) = floor(p * 4^n / q)

template <container C>
constexpr r<C> make_q(z<C> p, z<C> q) {
  return [p = std::move(p), q = std::move(q)](unsigned int n) -> coro::lazy<z<C>> {
    auto rr = mul_4exp(p, n);
    auto [quo, _] = floor_div(rr, q);
    co_return quo;
  };
}

template <container C>
constexpr r<C> add(r<C> x, r<C> y) {
  return [x = std::move(x), y = std::move(y)](unsigned int n) -> coro::lazy<z<C>> {
    // w = 4 when B = 4.
    auto xn = co_await x(n + 1);
    auto yn = co_await y(n + 1);
    co_return mul_4exp(xn + yn, -1);
  };
}

}  // namespace epx

#endif  // EPSILON_INC_R_HPP
