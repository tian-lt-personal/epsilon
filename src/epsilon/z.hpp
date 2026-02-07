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
    T v;
    T o;
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
    return result_t{.v = (w << s) + (ll & mask), .o = hh + (w >> s)};
  } else {
    constexpr unsigned s = sizeof(T) * 8u;
    MaxD l = lhs, r = rhs;
    auto prod = l * r;
    return result_t{.v = static_cast<T>(prod), .o = static_cast<T>(prod >> s)};
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

  if (s > 0) {
    u1 <<= s;
    u1 |= (u0 >> (64 - s));
    u0 <<= s;
  }
  uint32_t u1hi = static_cast<uint32_t>(u1 >> 32);
  uint32_t u1lo = static_cast<uint32_t>(u1 & 0xffffffffu);
  uint32_t u0hi = static_cast<uint32_t>(u0 >> 32);

  // D2. [Initialize j] (omitted by unfolded loops)
  // D3. [Calculate qhat]
  constexpr uint64_t b = 0x100000000ull;
  auto [qhat64, rhat32] = div_2ul(u1hi, u1lo, vn1);
  while (qhat64 == b || qhat64 * vn0 > (static_cast<uint64_t>(rhat32) << 32 | u0hi)) {
    --qhat64;
    rhat32 += vn1;
    if (rhat32 < vn1) break;  // continue if rhat < b.
  }

  // D4. [Multiply and subtract]
  // Compute p = qhat64 * v (p is 128-bit product of two 64-bit numbers)
  auto [p_low, p_high] = umul<uint64_t>(qhat64, v);
  
  // Subtract p from u1:u0
  // First, subtract low part from u0, tracking borrow to u1
  uint64_t borrow_to_u1 = (u0 < p_low) ? 1 : 0;
  u0 -= p_low;
  
  // Save original u1 to detect underflow
  uint64_t original_u1 = u1;
  uint64_t to_subtract_from_u1 = p_high + borrow_to_u1;
  
  // Subtract high part and borrow from u1
  u1 -= to_subtract_from_u1;
  
  // D6: [Add back] if underflow occurred
  // In unsigned arithmetic: a - b underflows if a < b
  // After the subtraction, we need to check if we borrowed beyond u1
  if (original_u1 < to_subtract_from_u1) {
    // Result is negative, so add back v and decrement quotient
    qhat64--;
    
    // Add v back to (u1, u0)
    uint64_t temp_u0 = u0 + v;
    uint64_t carry = (temp_u0 < v) ? 1 : 0;  // Detect overflow in u0 addition
    u0 = temp_u0;
    u1 += carry;
  }
  
  return result_t{.q = qhat64, .r = u0 >> s};
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
      auto [prod, o] = details::umul<default_digit_t>(a[i], b[j]);
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
