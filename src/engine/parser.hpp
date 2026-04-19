
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

// epx
#include "lexer.hpp"

namespace epx::script {

enum struct node_kind {
  expr,
  val,
  add_op,
  sub_op,
  mul_op,
  div_op,
  func_call,
};

class stmt {
 public:
  explicit stmt(node_kind kind) noexcept : kind_(kind) {}
  node_kind kind() noexcept { return kind_; }

 private:
  node_kind kind_;
};

namespace details {

namespace {
struct tu;
}  // namespace

struct ast_context {
  std::unique_ptr<std::pmr::monotonic_buffer_resource> pool;
};

struct value : stmt {
  value() : stmt(node_kind::val) {}
  std::variant<token_integer_literal, token_real_literal> token;
};
struct expr : stmt {
  expr() : stmt(node_kind::expr) {}
  expr* content = nullptr;
};

}  // namespace details

enum struct translate_ec {
  unknown,
};

struct mathscript {
  friend struct details::tu;

  std::vector<stmt*> statements;
  details::ast_context ctx_;
};

std::expected<mathscript, translate_ec> translate(std::string_view input) noexcept;

}  // namespace epx::script

#endif  // EPSILON_ENGINE_INC_PARSER_HPP
