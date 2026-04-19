// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <optional>

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
struct _EPX_EMPTY_BASE overloads : Ts... {
  using Ts::operator()...;
};

using lex_result = decltype(std::declval<epx::script::lexer>().operator()());

struct tu {
  lexer lex;
  token curtk;
  details::ast_context ctx;

  explicit tu(lexer l) noexcept : lex(l) {}
  std::optional<translate_ec> initialize() { return consume_token(); }
  std::expected<mathscript, translate_ec> parse() && noexcept {
    // todo: support multiple statements
    return parse_stmt().transform([&](stmt* stmt) { return mathscript{.statements = {stmt}, .ctx_ = std::move(ctx)}; });
  }
  std::expected<stmt*, translate_ec> parse_stmt() {
    return parse_expr().transform([](expr* expr) { return static_cast<stmt*>(expr); });
  }
  std::expected<expr*, translate_ec> parse_expr() {
    using result_type = std::expected<expr*, translate_ec>;
    return std::visit(overloads{
                          [](token_integer_literal) -> result_type { abort(); },
                          [](auto) -> result_type { return std::unexpected{translate_ec::unknown}; },
                      },
                      curtk);
  }

  std::optional<translate_ec> consume_token() {
    auto res = lex();
    if (!res.has_value()) [[unlikely]] {
      return translate_ec::unknown;
    }
    curtk = *res;
    return std::nullopt;
  };
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
