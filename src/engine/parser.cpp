// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <exception>
#include <optional>
#include <utility>

// epx
#include "parser.hpp"

namespace epx::script {

namespace details {
namespace {
#if defined(_MSC_VER) && !defined(__clang__)
#define _EPX_EMPTY_BASE __declspec(empty_bases)
#else
#define _EPX_EMPTY_BASE
#endif

template <class... Ts>
constexpr bool is(const token& tk) {
  return (std::holds_alternative<Ts>(tk) || ...);
};
constexpr bool is_bin_op(const token& tk) { return is<token_op_plus, token_op_minus, token_op_mul, token_op_div>(tk); }

template <class... Ts>
struct _EPX_EMPTY_BASE overloads : Ts... {
  using Ts::operator()...;
};

int get_precedence(const token& tk) noexcept {
  return std::visit(overloads{
                        [](token_op_plus) noexcept { return 10; },
                        [](token_op_minus) noexcept { return 10; },
                        [](token_op_mul) noexcept { return 20; },
                        [](token_op_div) noexcept { return 20; },
                        [](auto) noexcept -> int { std::terminate(); },
                    },
                    tk);
}

node_kind get_op_kind(const token& tk) noexcept {
  return std::visit(overloads{
                        [](token_op_plus) noexcept { return node_kind::add_expr; },
                        [](token_op_minus) noexcept { return node_kind::sub_expr; },
                        [](token_op_mul) noexcept { return node_kind::mul_expr; },
                        [](token_op_div) noexcept { return node_kind::div_expr; },
                        [](auto) noexcept -> node_kind { std::terminate(); },
                    },
                    tk);
}

using lex_result = std::invoke_result_t<epx::script::lexer>;
using parse_expr_result = std::expected<expr*, translate_ec>;

struct tu {
  lexer lex;
  token curtk;
  token nextk;
  details::ast_context ctx;

  explicit tu(lexer l) noexcept : lex(l) {}
  std::optional<translate_ec> initialize() {
    auto res = lex();
    if (!res.has_value()) {
      return translate_ec::unknown;
    }
    curtk = *res;
    res = lex();
    if (res.has_value()) {
      nextk = *res;
    }
    ctx.pool = std::make_unique<std::pmr::monotonic_buffer_resource>();
    return std::nullopt;
  }
  std::expected<mathscript, translate_ec> parse() && noexcept {
    // todo: support multiple statements
    return parse_stmt().transform([&](stmt* stmt) { return mathscript{.statements = {stmt}, .ctx_ = std::move(ctx)}; });
  }

 private:
  std::expected<stmt*, translate_ec> parse_stmt() noexcept {
    return parse_expr().transform([](expr* expr) { return static_cast<stmt*>(expr); });
  }
  std::expected<expr*, translate_ec> parse_expr() noexcept { return parse_expr_with_precedence(0); }
  std::expected<expr*, translate_ec> parse_expr_with_precedence(int min_precedence) noexcept {
    return parse_term().and_then([&](expr* left) -> parse_expr_result {
      while (is_bin_op(curtk) && get_precedence(curtk) > min_precedence) {
        auto op = curtk;
        consume_token();
        return parse_expr_with_precedence(get_precedence(op)).and_then([&](expr* right) -> parse_expr_result {
          return make<binop_expr>(get_op_kind(op), left, right);
        });
      }
      return left;
    });
  }
  std::expected<expr*, translate_ec> parse_term() noexcept {
    return std::visit(overloads{
                          [this](token_integer_literal integer) noexcept -> parse_expr_result {
                            consume_token();
                            return make<val_term>(integer);
                          },
                          [this](token_real_literal real) -> parse_expr_result {
                            consume_token();
                            return make<val_term>(real);
                          },
                          [this](token_id id) -> parse_expr_result {
                            consume_token();
                            return parse_func_call(id);
                          },
                          [this](token_lparen) -> parse_expr_result {
                            consume_token();
                            return parse_term().and_then([&](expr* expr) -> parse_expr_result {
                              consume_token();
                              return expr;
                            });
                          },
                          [](auto) -> parse_expr_result { return std::unexpected{translate_ec::unknown}; },
                      },
                      curtk);
  }
  std::expected<expr*, translate_ec> parse_func_call(token_id) noexcept { abort(); }
  void consume_token() noexcept {
    curtk = nextk;
    auto res = lex();
    nextk = res.has_value() ? *res : std::monostate{};
  };
  template <class T, class... Us>
  T* make(Us&&... params) noexcept {
    std::pmr::polymorphic_allocator<T> alloc{ctx.pool.get()};
    return alloc.new_object<T>(std::forward<Us>(params)...);  // let the exception terminates the program.
  }
};

}  // namespace
}  // namespace details

std::expected<mathscript, translate_ec> translate(std::string_view input) noexcept {
  details::tu tu{lexer{input}};
  if (auto err = tu.initialize(); err.has_value()) {
    return std::unexpected{*err};
  }
  return std::move(tu).parse();
}

}  // namespace epx::script
