// std headers
#include <exception>

// gtest headers
#include <gtest/gtest.h>

// epsilon headers
#include <coro.hpp>

namespace coro = epx::coro;

namespace {

struct test_error : std::exception {};

}  // namespace

namespace epxut {

TEST(coro_tests, lazy) {
  int side_effect = 0;
  [&]() -> coro::forget {
    EXPECT_EQ(side_effect, 0);
    auto two = co_await [&]() -> coro::lazy<int> {
      EXPECT_EQ(side_effect, 0);
      auto one = co_await [&]() -> coro::lazy<int> {
        EXPECT_EQ(side_effect, 0);
        side_effect = 1;
        co_return 1;
      }();
      EXPECT_EQ(side_effect, 1);
      EXPECT_EQ(one, 1);
      side_effect = 2;
      co_return 2;
    }();
    EXPECT_EQ(side_effect, 2);
    EXPECT_EQ(two, 2);
  }();
  EXPECT_EQ(side_effect, 2);
}

TEST(coro_tests, lazy_with_exception) {
  []() -> coro::forget {
    try {
      co_await [](bool x) -> coro::lazy<int> {
        if (x) {
          throw test_error{};
        } else {
          co_return 0;
        }
      }(true);
    } catch (const test_error&) {
      co_return;
    }
    std::terminate();
  }();
}

}  // namespace epxut
