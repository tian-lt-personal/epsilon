#ifndef EPSILON_INC_Z_HPP
#define EPSILON_INC_Z_HPP

// std
#include <algorithm>
#include <cassert>
#include <concepts>
#include <limits>
#include <ranges>

// epx
#include <t.hpp>

namespace epx {

namespace details {

template <class MaxD, std::unsigned_integral T>
constexpr T umul(T lhs, T rhs, T& o) {
  if constexpr (sizeof(T) >= sizeof(MaxD)) {
    constexpr T s = sizeof(T) * 4u;
    constexpr T mask = std::numeric_limits<T>::max() >> s;
    T a0 = lhs & mask, b0 = rhs & mask;
    T a1 = lhs >> s, b1 = rhs >> s;
    T ll = a0 * b0, lh = a0 * b1, hl = a1 * b0, hh = a1 * b1;
    T w = lh + (ll >> s);
    w += hl;
    if (w < hl) hh += mask + 1u;
    o = hh + (w >> s);
    return (w << s) + (ll & mask);
  } else {
    constexpr unsigned s = sizeof(T) * 8u;
    MaxD l = lhs, r = rhs;
    auto prod = l * r;
    o = static_cast<T>(prod >> s);
    return static_cast<T>(prod);
  }
}

}  // namespace details

template <container C>
constexpr bool is_zero(const z<C>& num) noexcept {
  return std::ranges::size(num.digits) == 0;
}

template <container C>
constexpr bool is_positive(const z<C>& num) noexcept {
  return num.sgn == sign::positive;
}

template <container C>
constexpr z<C>& normalize(z<C>& num) noexcept {
  using D = typename z<C>::digit_type;
  auto pos = std::find_if(num.digits.rbegin(), num.digits.rend(), [](D x) { return x != 0; });
  num.digits.resize(std::ranges::size(num.digits) - std::ranges::distance(num.digits.rbegin(), pos));
  num.sgn = num.digits.empty() == false ? num.sgn : sign::positive;
  return num;
}

template <container C>
constexpr z<C>& negate(z<C>& num) noexcept {
  num.sgn = is_positive(num) ? sign::negative : sign::positive;
  return num;
}

template <container C>
constexpr int cmp_n(const z<C>& lhs, const z<C>& rhs) {
  if (std::ranges::size(lhs.digits) != std::ranges::size(rhs.digits)) {
    return std::ranges::size(lhs.digits) < std::ranges::size(rhs.digits) ? -1 : 1;
  }

  auto res = std::lexicographical_compare_three_way(std::ranges::crbegin(lhs.digits), std::ranges::crend(lhs.digits),
                                                    std::ranges::crbegin(rhs.digits), std::ranges::crend(rhs.digits));
  if (res < 0) return -1;
  if (res > 0) return 1;
  return 0;
}

template <container C>
constexpr z<C> add_n(const z<C>& lhs, const z<C>& rhs) {
  using D = typename z<C>::digit_type;
  z<C> r;
  r.digits.reserve(std::max(std::ranges::size(lhs.digits), std::ranges::size(rhs.digits)) + 1);

  const auto& [a, b] = std::ranges::size(lhs.digits) <= std::ranges::size(rhs.digits)
                           ? std::tie(lhs.digits, rhs.digits)
                           : std::tie(rhs.digits, lhs.digits);
  size_t i = 0;
  D cy = 0;
  for (; i < std::ranges::size(a); ++i) {
    D s1 = a[i] + b[i];
    D cy1 = s1 < a[i];
    D s2 = s1 + cy;
    D cy2 = s2 < s1;
    cy = cy1 | cy2;
    r.digits.push_back(s2);
  }
  for (; i < std::ranges::size(b); ++i) {
    D s = b[i] + cy;
    cy = s < b[i];
    r.digits.push_back(s);
  }
  if (cy > 0) {
    r.digits.push_back(1u);
  }
  return r;
}

template <container C>
constexpr z<C> sub_n(const z<C>& lhs, const z<C>& rhs) {
  assert(std::ranges::size(lhs.digits) >= std::ranges::size(rhs.digits));
  using D = typename z<C>::digit_type;
  z<C> r;
  const auto& a = lhs.digits;
  const auto& b = rhs.digits;

  r.digits.reserve(a.size());
  D borrow = 0;
  size_t i = 0;

  for (; i < b.size(); ++i) {
    D d1 = a[i];
    D d2 = b[i];
    D diff = d1 - d2 - borrow;
    borrow = borrow ? (d1 <= d2) : (d1 < d2);
    r.digits.push_back(diff);
  }
  for (; i < a.size(); ++i) {
    D d1 = a[i];
    D diff = d1 - borrow;
    borrow = (d1 < borrow);
    r.digits.push_back(diff);
  }
  assert(borrow == 0);
  normalize(r);
  return r;
}

template <container C>
constexpr z<C> mul_n(const z<C>& lhs, const z<C>& rhs) {
  using D = typename z<C>::digit_type;
  if (is_zero(lhs) || is_zero(rhs)) {
    return z<C>{};
  }

  z<C> r;
  r.digits.resize(std::ranges::size(lhs.digits) + std::ranges::size(rhs.digits));

  const auto& a = lhs.digits;
  const auto& b = rhs.digits;
  for (size_t j = 0; j < std::ranges::size(b); ++j) {
    D cy = 0;
    for (size_t i = 0; i < std::ranges::size(a); ++i) {
      D o;
      auto prod = details::umul<default_digit_t>(a[i], b[j], o);
      prod += cy;
      cy = (prod < cy ? 1u : 0u) + o;
      r.digits[i + j] += prod;
      if (r.digits[i + j] < prod) ++cy;
    }
    r.digits[j + std::ranges::size(a)] = cy;
  }
  normalize(r);
  return r;
}

template <container C>
constexpr z<C> add(const z<C>& lhs, const z<C>& rhs) {
  z<C> r;
  if (lhs.sgn == rhs.sgn) {
    r = add_n(lhs, rhs);
    r.sgn = lhs.sgn;
  } else {
    const auto* minuend = &lhs;
    const auto* substrahend = &rhs;
    if (cmp_n(lhs, rhs) < 0) {
      minuend = &rhs;
      substrahend = &lhs;
    }
    r = sub_n(*minuend, *substrahend);
    r.sgn = std::ranges::size(r.digits) > 0 ? minuend->sgn : sign::positive;
  }
  return r;
}

template <container C>
constexpr z<C> sub(const z<C>& lhs, z<C> rhs) {
  return add(lhs, negate(rhs));
}

template <container C>
constexpr z<C> mul(const z<C>& lhs, const z<C>& rhs) {
  z<C> r = mul_n(lhs, rhs);
  if (is_zero(r)) {
    return r;
  }
  r.sgn = lhs.sgn == rhs.sgn ? sign::positive : sign::negative;
  return r;
}

}  // namespace epx

#endif
