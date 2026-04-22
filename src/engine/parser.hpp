// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_ENGINE_INC_PARSER_HPP
#define EPSILON_ENGINE_INC_PARSER_HPP

// std
#include <expected>
#include <memory>
#include <memory_resource>
#include <string_view>
#include <variant>
#include <vector>

// epx
#include "lexer.hpp"

namespace epx::script {

enum struct node_kind {
  val,
  paren_expr,
  add_expr,
  sub_expr,
  mul_expr,
  div_expr,
  func_call,
};

namespace details {

namespace {
struct tu;
}  // namespace

struct ast_context {
  std::unique_ptr<std::pmr::monotonic_buffer_resource> pool;
};

struct expr {
 public:
  explicit expr(node_kind kind) noexcept : kind_(kind) {}
  node_kind kind() const noexcept { return kind_; }

 private:
  node_kind kind_;
};

struct paren_expr : expr {
  paren_expr() noexcept : expr(node_kind::paren_expr) {}
  explicit paren_expr(expr* inner_content) noexcept : expr(node_kind::paren_expr), inner(inner_content) {}
  expr* inner = nullptr;
};
struct val_term : expr {
  val_term(std::variant<token_integer_literal, token_real_literal> value) : expr(node_kind::val), val(value) {}
  std::variant<token_integer_literal, token_real_literal> val;
};
struct func_call : expr {
  explicit func_call(token_id funcname, std::vector<expr*> parameters) noexcept
      : expr(node_kind::func_call), id(funcname), params(std::move(parameters)) {}
  token_id id;
  std::vector<expr*> params;
};
struct binop_expr : expr {
  explicit binop_expr(node_kind op, expr* left_expr, expr* right_expr) noexcept
      : expr(op), left(left_expr), right(right_expr) {}
  expr* left = nullptr;
  expr* right = nullptr;
};

}  // namespace details

enum struct translate_ec {
  unknown,
};

using stmt = std::variant<details::expr*>;

struct mathscript {
  friend struct details::tu;

  std::vector<stmt> statements;
  details::ast_context ctx_;
};

std::expected<mathscript, translate_ec> translate(std::string_view input) noexcept;

}  // namespace epx::script

#endif  // EPSILON_ENGINE_INC_PARSER_HPP
