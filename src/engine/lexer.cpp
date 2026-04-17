// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <cassert>
#include <cctype>

// epx
#include "lexer.hpp"

namespace epx::script {

namespace {
enum struct lex_r_ec {
  others,
  double_dot,
};
}

bool isdelim(char ch) {
  return std::isspace(static_cast<unsigned char>(ch)) || ch == '.' || ch == '+' || ch == '-' || ch == '*' ||
         ch == '/' || ch == '%' || ch == '(' || ch == ')' || ch == '=' || ch == ',';
}

std::expected<token, token_ec> lexer::operator()() noexcept {
  const char* end = input_.data() + input_.length();
  auto getch = [&]() noexcept -> std::expected<char, token_ec> {
    while (cursor_ < end && std::isspace(static_cast<unsigned char>(*cursor_))) {
      ++cursor_;
    }
    if (cursor_ == end) {
      return std::unexpected{token_ec::eof};
    }
    return *cursor_;
  };
  return getch().and_then([&](char ch) noexcept -> std::expected<token, token_ec> {
    auto try_lex_z = [&](const char* cursor) noexcept -> std::expected<token_integer_literal, int> {
      assert(cursor < end);
      if (*cursor == 'z') {
        if (cursor + 2 < end && *(cursor + 1) == '\'') {
          cursor += 2;
        } else {
          return std::unexpected{1};
        }
      }
      size_t len = 0;
      auto data = cursor;
      while (cursor < end && std::isdigit(static_cast<unsigned char>(*cursor))) {
        ++cursor;
        ++len;
      }
      if (len == 0) {
        return std::unexpected{2};
      } else if (cursor < end && (*cursor == '.' || !isdelim(*cursor))) {
        return std::unexpected{3};
      } else {
        cursor_ = cursor;
        return token_integer_literal{.raw = {data, len}};
      }
    };
    auto try_lex_r = [&](const char* cursor) noexcept -> std::expected<token_real_literal, lex_r_ec> {
      assert(cursor < end);
      if (*cursor == 'r') {
        if (cursor + 2 < end && *(cursor + 1) == '\'') {
          cursor += 2;
        } else {
          return std::unexpected{lex_r_ec::others};
        }
      }
      bool has_dot = false;
      size_t len = 0;
      auto data = cursor;
      while (cursor < end && (std::isdigit(static_cast<unsigned char>(*cursor)) || *cursor == '.')) {
        if (*cursor == '.') {
          if (!has_dot) {
            has_dot = true;
          } else {
            return std::unexpected{lex_r_ec::double_dot};
          }
        }
        ++cursor;
        ++len;
      }
      if (len == 0) {
        return std::unexpected{lex_r_ec::others};
      } else if (cursor < end && !isdelim(*cursor)) {
        return std::unexpected{lex_r_ec::others};
      } else if (has_dot && (*cursor == '.' || len == 1)) {
        return std::unexpected{lex_r_ec::others};
      } else {
        cursor_ = cursor;
        return token_real_literal{.raw = {data, len}};
      }
    };
    auto try_lex_id = [&](const char* cursor) noexcept -> std::expected<token_id, token_ec> {
      assert(cursor < end);
      assert(std::isalpha(static_cast<unsigned char>(*cursor)));
      auto data = cursor++;
      size_t len = 1;
      while (cursor < end &&
             (std::isalpha(static_cast<unsigned char>(*cursor)) || std::isdigit(static_cast<unsigned char>(*cursor)))) {
        ++cursor;
        ++len;
      }
      if (cursor == end || isdelim(*cursor)) {
        cursor_ = cursor;
        return token_id{.raw = {data, len}};
      } else {
        return std::unexpected{token_ec::bad_input};
      }
    };
    switch (ch) {
      case '=':
        ++cursor_;
        return token_op_eq{};
      case ',':
        ++cursor_;
        return token_op_comma{};
      case '+':
        ++cursor_;
        return token_op_plus{};
      case '-':
        ++cursor_;
        return token_op_minus{};
      case '*':
        ++cursor_;
        return token_op_mul{};
      case '/':
        ++cursor_;
        return token_op_div{};
      case '(':
        ++cursor_;
        return token_lparen{};
      case ')':
        ++cursor_;
        return token_rparen{};
      case '%':
        ++cursor_;
        return token_op_percent{};
      case 'z': {
        const char* next = cursor_ + 1;
        if (next == end) {
          auto data = cursor_++;
          return token_id{.raw = {data, 1uz}};
        } else if (*next == '\'') {
          return try_lex_z(cursor_)
              .transform([](token_integer_literal literal) { return token{literal}; })
              .transform_error([](int) { return token_ec::bad_input; });
        } else {
          return try_lex_id(cursor_).transform([](token_id id) { return token{id}; });
        }
      }
      case 'r': {
        const char* next = cursor_ + 1;
        if (next == end) {
          auto data = cursor_++;
          return token_id{.raw = {data, 1uz}};
        } else if (*next == '\'') {
          return try_lex_r(cursor_)
              .transform([](token_real_literal literal) { return token{literal}; })
              .transform_error([](lex_r_ec) { return token_ec::bad_input; });
        } else {
          return try_lex_id(cursor_).transform([](token_id id) { return token{id}; });
        }
      }
    }

    if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '.') {
      return try_lex_r(cursor_)
          .transform([](token_real_literal literal) { return token{literal}; })
          .or_else([&](lex_r_ec err) -> std::expected<token, lex_r_ec> {
            if (ch == '.' && err != lex_r_ec::double_dot) {
              ++cursor_;
              return token_op_dot{};
            }
            return std::unexpected{err};
          })
          .transform_error([](lex_r_ec) { return token_ec::bad_input; });
    } else if (std::isalpha(static_cast<unsigned char>(ch))) {
      return try_lex_id(cursor_).transform([](token_id id) { return token{id}; });
    }
    return std::unexpected{token_ec::bad_input};
  });
}

}  // namespace epx::script
