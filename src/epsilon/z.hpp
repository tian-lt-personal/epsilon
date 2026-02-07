#ifndef EPSILON_INC_Z_HPP
#define EPSILON_INC_Z_HPP

// std
#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <limits>
#include <ranges>

// epx
#include <t.hpp>

namespace epx {

namespace details {

template <class MaxD, std::unsigned_integral T>
constexpr auto umul(T lhs, T rhs) {
  struct result_t {
    T p0;
    T p1;
  };

  if constexpr (sizeof(T) >= sizeof(MaxD)) {
    constexpr T s = sizeof(T) * 4u;
    constexpr T mask = std::numeric_limits<T>::max() >> s;
    T a0 = lhs & mask, b0 = rhs & mask;
    T a1 = lhs >> s, b1 = rhs >> s;
    T ll = a0 * b0, lh = a0 * b1, hl = a1 * b0, hh = a1 * b1;
    T w = lh + (ll >> s);
    w += hl;
    if (w < hl) hh += mask + 1u;
    return result_t{.p0 = (w << s) + (ll & mask), .p1 = hh + (w >> s)};
  } else {
    constexpr unsigned s = sizeof(T) * 8u;
    MaxD l = lhs, r = rhs;
    auto prod = l * r;
    return result_t{.p0 = static_cast<T>(prod), .p1 = static_cast<T>(prod >> s)};
  }
}

// u0 - LSD, u1 - MSD
constexpr auto div_2ul(uint32_t u0, uint32_t u1, uint32_t v) {
  struct result_t {
    uint64_t q;
    uint32_t r;
  };
  assert(v != 0);
  uint64_t u = (static_cast<uint64_t>(u1) << 32) | u0;
  return result_t{.q = u / v, .r = static_cast<uint32_t>(u % v)};
}

// u0 - LSD, u1 - MSD
constexpr auto div_2ull(uint64_t u0, uint64_t u1, uint64_t v) {
  struct result_t {
    uint64_t q;
    uint64_t r;
  };
  if (v == 0) [[unlikely]] {
    throw divide_by_zero{};
  }
  if (u1 >= v) [[unlikely]] {
    throw overflow_error{};
  }
  if (u1 == 0) {
    return result_t{.q = u0 / v, .r = u0 % v};
  }

  // D1. [Normalize] (and break them down int 32-bit parts)
  int s = std::countl_zero(v);
  assert(s < 64);
  v <<= s;

  uint32_t vn1 = static_cast<uint32_t>(v >> 32);
  uint32_t vn0 = static_cast<uint32_t>(v & 0xFFFFFFFFu);

  uint32_t un4, un3, un2, un1, un0;
  if (s > 0) {
    un4 = static_cast<uint32_t>(u1 >> (32 - s));
    un3 = static_cast<uint32_t>((u1 << s) | (u0 >> (32 - s)));
    un2 = static_cast<uint32_t>((u0 << s) & 0xFFFFFFFFu);
    un1 = static_cast<uint32_t>(u0 >> (32 - s));
    un0 = static_cast<uint32_t>((u0 << s) & 0xFFFFFFFFu);
  }

  // D2. [Initialize j] (omitted by unfolded loops)
  // D3. [Calculate qhat]
  constexpr uint64_t b = 0x100000000ull;
  auto [qhat64, rhat32] = div_2ul(un4, un3, vn1);
  while (qhat64 == b || qhat64 * vn0 > (static_cast<uint64_t>(rhat32) << 32 | un2)) {
    --qhat64;
    rhat32 += vn1;
    if (rhat32 < vn1) break;  // continue if rhat < b.
  }
  assert(qhat64 < b);

  // D4. [Multiply and subtract]
  auto [p0, p1] = umul<uint64_t>(qhat64, v);

  int64_t borrow = 0;
  int64_t diff = static_cast<int64_t>(un2) - (p0 & 0xFFFFFFFF);
  un2 = static_cast<uint32_t>(diff);
  borrow = (p0 >> 32) - (diff >> 63);

  diff = static_cast<int64_t>(un3) - (p1 & 0xFFFFFFFF) - borrow;
  un3 = static_cast<uint32_t>(diff);
  borrow = (p1 >> 32) - (diff >> 63);

  diff = static_cast<int64_t>(un4) - borrow;
  un4 = static_cast<uint32_t>(diff);

  if (diff < 0) {
    --qhat64;
    uint64_t carry = static_cast<uint64_t>(un2) + vn0;
    un2 = static_cast<uint32_t>(carry);
    carry = static_cast<uint64_t>(un3) + vn1 + (carry >> 32);
    un3 = static_cast<uint32_t>(carry);
    un4 += (carry >> 32);
  }
  assert(qhat64 < b);
  auto q1 = qhat64;

  //TODO

  return result_t{.q = q1 << 32, .r = u0 >> s};
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
      auto [p0, p1] = details::umul<default_digit_t>(a[i], b[j]);
      p0 += cy;
      cy = (p0 < cy ? 1u : 0u) + p1;
      r.digits[i + j] += p0;
      if (r.digits[i + j] < p0) ++cy;
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
