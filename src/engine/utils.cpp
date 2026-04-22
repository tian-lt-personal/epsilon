// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <format>
#include <optional>
#include <string_view>

// epx
#include "tmp.hpp"
#include "utils.hpp"

namespace epx::script {

namespace details {
namespace {
void dump_expr(std::string& out, const expr* e, unsigned depth);
void fillmargin(std::string& out, unsigned depth) { out.append(depth * 2uz, ' '); }
void scope_beg(std::string& out, unsigned depth, std::optional<std::string_view> scope_name = std::nullopt) {
  fillmargin(out, depth);
  if (scope_name.has_value()) {
    out += std::format("{}={{\n", *scope_name);
  } else {
    out += "{\n";
  }
}
void scope_end(std::string& out, unsigned depth) {
  fillmargin(out, depth);
  out += "}\n";
}

void dump_val(std::string& out, const val_term* val, unsigned depth) {
  std::visit(tmp::overloads{
                 [&](token_integer_literal literal) {
                   fillmargin(out, depth);
                   out += std::format("integer: {}\n", literal.raw);
                 },
                 [&](token_real_literal literal) {
                   fillmargin(out, depth);
                   out += std::format("real: {}\n", literal.raw);
                 },
                 [&](auto) {},
             },
             val->val);
}
void dump_binop_expr(std::string& out, const binop_expr* expr, unsigned depth) {
  char opch;
  switch (expr->kind()) {
    case node_kind::add_expr:
      opch = '+';
      break;
    case node_kind::sub_expr:
      opch = '-';
      break;
    case node_kind::mul_expr:
      opch = '*';
      break;
    case node_kind::div_expr:
      opch = '/';
      break;
    default:
      opch = '?';
  }
  fillmargin(out, depth);
  out += std::format("[binop_expr {}]\n", opch);
  scope_beg(out, depth, "left");
  dump_expr(out, expr->left, depth + 1);
  scope_end(out, depth);
  scope_beg(out, depth, "right");
  dump_expr(out, expr->right, depth + 1);
  scope_end(out, depth);
}
void dump_expr(std::string& out, const expr* e, unsigned depth) {
  scope_beg(out, depth);
  switch (e->kind()) {
    case node_kind::add_expr:
    case node_kind::sub_expr:
    case node_kind::mul_expr:
    case node_kind::div_expr:
      dump_binop_expr(out, static_cast<const binop_expr*>(e), depth + 1);
      break;
    case node_kind::val:
      dump_val(out, static_cast<const val_term*>(e), depth + 1);
      break;
    default:
      fillmargin(out, depth + 1);
      out += "[unknown expr]\n";
  }
  scope_end(out, depth);
}
void dump_stmt(std::string& out, const stmt stmt, unsigned depth) { dump_expr(out, std::get<expr*>(stmt), depth); }
}  // namespace
}  // namespace details

std::string dump(const mathscript& script) noexcept {
  std::string result;
  for (auto& stmt : script.statements) {
    details::dump_stmt(result, stmt, 0);
  }
  return result;
}

}  // namespace epx::script
