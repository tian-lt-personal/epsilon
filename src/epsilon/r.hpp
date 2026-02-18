// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_INC_R_HPP
#define EPSILON_INC_R_HPP

// std
#include <algorithm>
#include <functional>
#include <limits>

// epx
#include "coro.hpp"
#include "z.hpp"

namespace epx {

template <container C>
using r = std::function<coro::lazy<z<C>>(int)>;  // q (n) = floor(p * 4^n / q)

template <container C>
constexpr r<C> make_q(z<C> p, z<C> q) {
  return [p = std::move(p), q = std::move(q)](int n) -> coro::lazy<z<C>> {
    auto rr = mul_4exp(p, n);
    auto [quo, _] = floor_div(rr, q);
    co_return quo;
  };
}

template <container C>
constexpr r<C> add(r<C> x, r<C> y) {
  return [x = std::move(x), y = std::move(y)](int n) -> coro::lazy<z<C>> {
    if (n == std::numeric_limits<int>::max()) [[unlikely]] {
      throw precision_overflow_error{};
    }
    // w = 4 when B = 4.
    auto xn = co_await x(n + 1);
    auto yn = co_await y(n + 1);
    co_return mul_4exp(xn + yn, -1);
  };
}

template <container C>
constexpr r<C> opp(r<C> x) {
  return [x = std::move(x)](int n) -> coro::lazy<z<C>> {
    auto xn = co_await x(n);
    negate(xn);
    co_return xn;
  };
}

template <container C>
coro::lazy<int> msd(r<C> x, int max) {
  const auto zero = details::zero<C>();
  const auto one = details::one<C>();
  const auto base = details::base<C>();

  auto x0 = co_await x(0);
  if (cmp_n(x0, base) > 0) {  // 4 < x0
    int i = -1;
    for (;;) {
      auto xi = co_await x(i);
      if (cmp_n(xi, one) <= 0) {
        co_return i + 1;
      }
      --i;  // TODO: use faster search algo.
    }
  } else if (cmp_n(x0, one) > 0) {  // 1 < x0 <= 4
    co_return 0;
  } else if (cmp_n(x0, zero) > 0) {  // 0 < x0 <= 1
    co_return 1;
  } else {  // x0 <= 0
    int i = 0;
    auto xi = x0;
    for (;;) {
      auto c = cmp_n(xi, one);
      if (c > 0 || i >= max) {
        co_return i;
      }
      ++i;  // TODO: use faster search algo.
      if (i > max_msd<global_config_tag>) [[unlikely]] {
        throw msd_overflow_error{};
      }
      xi = co_await x(i);
    }
  }
}

template <container C>
coro::lazy<int> msd(r<C> x) {
  return msd(std::move(x), std::numeric_limits<int>::max());
}

template <container C>
constexpr r<C> mul(r<C> x, r<C> y) {
  struct {
    r<C> x;
    r<C> y;

    coro::lazy<z<C>> operator()(int n) {
      if (n > std::numeric_limits<int>::max() - 3) [[unlikely]] {
        throw precision_overflow_error{};
      }
      const auto one = details::one<C>();
      int max_msd = n + 3 - (n + 2) / 2;
      int px = n - (co_await msd(y, max_msd)) + 3;
      int py = n - (co_await msd(x, max_msd)) + 3;
      auto xpx = co_await x(px);
      auto ypy = co_await y(py);
      auto xy = mul_4exp(add_n(mul_n(xpx, ypy), one), n - px - py);
      xy.sgn = xpx.sgn == ypy.sgn ? sign::positive : sign::negative;
      co_return xy;
    }
  } expr{.x = std::move(x), .y = std::move(y)};
  return expr;
}

template <container C>
constexpr r<C> inv(r<C> x) {
  struct {
    r<C> x;

    coro::lazy<z<C>> operator()(int n) {
      int msdx = co_await msd(x);
      if (n <= -msdx) {
        co_return z<C>{};  // zero
      }

      assert(std::numeric_limits<int>::max() > 2 * max_msd<global_config_tag> + 1);
      if (n > std::numeric_limits<int>::max() - 2 * max_msd<global_config_tag> - 1) [[unlikely]] {
        throw precision_overflow_error{};
      }

      int k = n + 2 * msdx + 1;
      auto p = mul_4exp(details::one<C>(), k + n);
      auto q = add(co_await x(k), details::one<C>());
      auto [quo, _] = ceil_div(p, q);
      co_return add(quo, details::one<C>());
    }
  } expr{.x = std::move(x)};
  return expr;
}

}  // namespace epx

#endif  // EPSILON_INC_R_HPP
