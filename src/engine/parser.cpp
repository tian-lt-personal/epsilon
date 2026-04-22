// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <exception>
#include <optional>
#include <utility>

// epx
#include "parser.hpp"
#include "tmp.hpp"

namespace epx::script {

namespace details {
namespace {

template <class... Ts>
constexpr bool is(const token& tk) {
  return (std::holds_alternative<Ts>(tk) || ...);
};
constexpr bool is_bin_op(const token& tk) { return is<token_op_plus, token_op_minus, token_op_mul, token_op_div>(tk); }

int get_precedence(const token& tk) noexcept {
  return std::visit(tmp::overloads{
                        [](token_op_plus) noexcept { return 10; },
                        [](token_op_minus) noexcept { return 10; },
                        [](token_op_mul) noexcept { return 20; },
                        [](token_op_div) noexcept { return 20; },
                        [](auto) noexcept -> int { std::terminate(); },
                    },
                    tk);
}

node_kind get_op_kind(const token& tk) noexcept {
  return std::visit(tmp::overloads{
                        [](token_op_plus) noexcept { return node_kind::add_expr; },
                        [](token_op_minus) noexcept { return node_kind::sub_expr; },
                        [](token_op_mul) noexcept { return node_kind::mul_expr; },
                        [](token_op_div) noexcept { return node_kind::div_expr; },
                        [](auto) noexcept -> node_kind { std::terminate(); },
                    },
                    tk);
}

using expr_result_type = std::expected<expr*, translate_ec>;

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
    using result_type = std::expected<mathscript, translate_ec>;
    // todo: support multiple statements
    return parse_expr().and_then([&](expr* e) -> result_type {
      if (drained()) {
        return mathscript{.statements = {stmt{e}}, .ctx_ = std::move(ctx)};
      } else {
        return std::unexpected{translate_ec::unknown};
      }
    });
  }

 private:
  bool drained() const noexcept { return lex.drained() && is<std::monostate>(curtk); }
  std::expected<expr*, translate_ec> parse_expr() noexcept { return parse_expr_with_precedence(0); }
  std::expected<expr*, translate_ec> parse_expr_with_precedence(int min_precedence) noexcept {
    return parse_term().and_then([&](expr* left) -> expr_result_type {
      while (is_bin_op(curtk) && get_precedence(curtk) > min_precedence) {
        auto op = curtk;
        consume_token();
        if (auto right_result = parse_expr_with_precedence(get_precedence(op)); right_result.has_value()) {
          left = make<binop_expr>(get_op_kind(op), left, *right_result);
        } else {
          return std::unexpected{right_result.error()};
        }
      }
      return left;
    });
  }
  std::expected<expr*, translate_ec> parse_term() noexcept {
    return std::visit(tmp::overloads{
                          [this](token_integer_literal integer) noexcept -> expr_result_type {
                            consume_token();
                            return make<val_term>(integer);
                          },
                          [this](token_real_literal real) -> expr_result_type {
                            consume_token();
                            return make<val_term>(real);
                          },
                          [this](token_id id) -> expr_result_type {
                            consume_token();
                            return parse_func_call(id);
                          },
                          [this](token_lparen) -> expr_result_type {
                            consume_token();
                            return parse_term().and_then([&](expr* expr) -> expr_result_type {
                              consume_token();
                              return expr;
                            });
                          },
                          [](auto) -> expr_result_type { return std::unexpected{translate_ec::unknown}; },
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
    return alloc.template new_object<T>(std::forward<Us>(params)...);  // let the exception terminate the program.
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
