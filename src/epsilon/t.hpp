#ifndef EPSILON_INC_T_HPP
#define EPSILON_INC_T_HPP

// std
#include <concepts>
#include <cstdint>
#include <ranges>
#include <type_traits>

namespace epx {

using default_digit_t = uint32_t;
using max_digit_t = uint32_t;
using default_container_t = std::vector<default_digit_t>;

struct bad_digit_t;
template <class D>
using wide_digit_t =
    std::conditional_t<sizeof(D) == sizeof(uint8_t), uint16_t,
                       std::conditional_t<sizeof(D) == sizeof(uint16_t), uint32_t,
                                          std::conditional_t<sizeof(D) == sizeof(uint32_t), uint64_t, bad_digit_t>>>;

template <class T>
concept container = std::ranges::random_access_range<T> &&  //
                    std::ranges::sized_range<T> &&          //
                    requires {
                      typename T::value_type;
                      std::is_unsigned_v<typename T::value_type>;
                    } && requires(T c) {
                      c.push_back(typename T::value_type{});
                      c.reserve(size_t{});
                      sizeof(typename T::value_type) < sizeof(max_digit_t);
                    };

enum class sign : uint8_t { positive, negative };

template <container C = default_container_t>
struct z {
  using container_type = C;
  using digit_type = typename C::value_type;
  std::strong_ordering operator<=>(const z&) const = default;

  C digits;  // // least significant digits (LSD)
  sign sgn = sign::positive;
};

struct divide_by_zero : public std::runtime_error {
  divide_by_zero() : std::runtime_error("epx: divide by zero") {}
};

struct overflow_error : public std::runtime_error {
  overflow_error() : std::runtime_error("epx: overflow error") {}
};

}  // namespace epx

#endif
