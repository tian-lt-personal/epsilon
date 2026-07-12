// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <exception>
#include <optional>
#include <utility>

// epx
#include "chars.hpp"
#include "tmp.hpp"

// engine
#include "parser.hpp"

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
                            if (is<token_lparen>(curtk)) {
                              return parse_func_call(id);
                            } else {
                              return make<val_term>(id);
                            }
                          },
                          [this](token_lparen) -> expr_result_type {
                            consume_token();
                            return parse_expr().and_then([&](expr* inner) -> expr_result_type {
                              if (!is<token_rparen>(curtk)) {
                                return std::unexpected{translate_ec::unknown};
                              }
                              consume_token();
                              return make<paren_expr>(inner);
                            });
                          },
                          [](auto) -> expr_result_type { return std::unexpected{translate_ec::unknown}; },
                      },
                      curtk);
  }
  std::expected<expr*, translate_ec> parse_func_call(token_id id) noexcept {
    if (!is<token_lparen>(curtk)) {
      std::terminate();
    }
    std::vector<expr*> param_seq;
    do {
      consume_token();
      if (is<token_rparen>(curtk)) {
        break;
      }
      auto param_res = parse_expr();
      if (!param_res.has_value()) {
        return std::unexpected{param_res.error()};
      }
      param_seq.push_back(*param_res);
    } while (is<token_op_comma>(curtk));
    if (!is<token_rparen>(curtk)) {
      return std::unexpected{translate_ec::unknown};
    }
    consume_token();
    return make<func_call>(id, std::move(param_seq));
  }
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

std::expected<parsed_real, translate_ec> parse_real_literal(std::string_view raw) noexcept {
  auto dot_pos = raw.find('.');
  if (dot_pos == std::string_view::npos) {
    auto num = epx::try_from_chars<default_container_type>(raw);
    if (!num.has_value()) return std::unexpected{translate_ec::unknown};
    return parsed_real{.num = std::move(*num), .den = epx::details::one<default_container_type>()};
  }
  auto int_str = raw.substr(0, dot_pos);
  auto frac_str = raw.substr(dot_pos + 1);

  auto int_part = int_str.empty() ? epx::details::zero<default_container_type>()
                                  : epx::try_from_chars<default_container_type>(int_str);
  auto frac_part = frac_str.empty() ? epx::details::zero<default_container_type>()
                                    : epx::try_from_chars<default_container_type>(frac_str);
  if (!int_part.has_value() || !frac_part.has_value()) {
    return std::unexpected{translate_ec::unknown};
  }

  auto k = static_cast<unsigned>(frac_str.length());
  auto ten_pow = epx::details::pow10<default_container_type>(k);
  auto num = epx::add(epx::mul(*int_part, ten_pow), *frac_part);
  return parsed_real{.num = std::move(num), .den = std::move(ten_pow)};
}

}  // namespace epx::script
