// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <functional>

// epx
#include "coro.hpp"
#include "z.hpp"

namespace epx {

template <container C>
using r = std::function<coro::lazy<z<C>>(unsigned int)>;  // q (n) = floor(p * 4^n / q)

template <container C>
r<C> make_q(z<C> p, z<C> q) {
  return [p = std::move(p), q = std::move(q)](unsigned int n) -> coro::lazy<z<C>> {
    auto rr = mul_4exp(p, n);
    auto [quo, _] = floor_div(rr, q);
    co_return quo;
  };
}

}  // namespace epx
