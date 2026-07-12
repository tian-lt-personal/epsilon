// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <expected>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>

// epx
#include "chars.hpp"
#include "tmp.hpp"

// engine
#include "machine.hpp"

namespace epx::script {

namespace {

using unary_func_t = real (*)(real);
using binary_func_t = real (*)(real, real);

real make_pi() {
  auto one = epx::details::one<default_container_type>();
  auto four = epx::make_q(epx::create<default_container_type>(4), one);
  auto one_r = epx::make_q(one, epx::details::one<default_container_type>());
  return epx::mul(std::move(four), epx::arctan(std::move(one_r)));
}

real make_e() {
  auto one = epx::details::one<default_container_type>();
  auto one_r = epx::make_q(one, epx::details::one<default_container_type>());
  return epx::exp(std::move(one_r));
}

const std::map<std::string, unary_func_t, std::less<>> unary_funcs = {
    {"sin", [](real x) -> real { return epx::sin(std::move(x)); }},
    {"cos", [](real x) -> real { return epx::cos(std::move(x)); }},
    {"tan", [](real x) -> real { return epx::tan(std::move(x)); }},
    {"exp", [](real x) -> real { return epx::exp(std::move(x)); }},
    {"log", [](real x) -> real { return epx::log(std::move(x)); }},
    {"ln", [](real x) -> real { return epx::log(std::move(x)); }},
    {"arcsin", [](real x) -> real { return epx::arcsin(std::move(x)); }},
    {"arccos", [](real x) -> real { return epx::arccos(std::move(x)); }},
    {"arctan", [](real x) -> real { return epx::arctan(std::move(x)); }},
    {"sinh", [](real x) -> real { return epx::sinh(std::move(x)); }},
    {"cosh", [](real x) -> real { return epx::cosh(std::move(x)); }},
    {"tanh", [](real x) -> real { return epx::tanh(std::move(x)); }},
    {"arcsinh", [](real x) -> real { return epx::arcsinh(std::move(x)); }},
    {"arccosh", [](real x) -> real { return epx::arccosh(std::move(x)); }},
    {"arctanh", [](real x) -> real { return epx::arctanh(std::move(x)); }},
    {"sqrt", [](real x) -> real { return epx::root(std::move(x), 2); }},
};

const std::map<std::string, binary_func_t, std::less<>> binary_funcs = {
    {"pow", [](real x, real y) -> real { return epx::pow(std::move(x), std::move(y)); }},
    {"log_base", [](real x, real y) -> real { return epx::log_base(std::move(x), std::move(y)); }},
};

const std::map<std::string, std::function<real()>, std::less<>> builtin_constants = {
    {"pi", make_pi},
    {"e", make_e},
};

class eval_session {
 public:
  explicit eval_session(eval_context* ctx) : ctx_(ctx) {}

  std::expected<eval_result, machine_ec> eval(const expr* e) const;

 private:
  std::expected<eval_result, machine_ec> eval_val(const details::val_term* v) const;
  std::expected<eval_result, machine_ec> eval_paren(const details::paren_expr* p) const;
  std::expected<eval_result, machine_ec> eval_binop(const details::binop_expr* b) const;
  std::expected<eval_result, machine_ec> eval_funcall(const details::func_call* f) const;

  eval_context* ctx_;
};

std::expected<eval_result, machine_ec> eval_session::eval(const expr* e) const {
  switch (e->kind()) {
    case node_kind::val:
      return eval_val(static_cast<const details::val_term*>(e));
    case node_kind::paren_expr:
      return eval_paren(static_cast<const details::paren_expr*>(e));
    case node_kind::add_expr:
    case node_kind::sub_expr:
    case node_kind::mul_expr:
    case node_kind::div_expr:
      return eval_binop(static_cast<const details::binop_expr*>(e));
    case node_kind::func_call:
      return eval_funcall(static_cast<const details::func_call*>(e));
  }
  return std::unexpected{machine_ec::unknown};
}

std::expected<eval_result, machine_ec> eval_session::eval_val(const details::val_term* v) const {
  return std::visit(
      tmp::overloads{
          [](const token_integer_literal& lit) -> std::expected<eval_result, machine_ec> {
            auto num = epx::try_from_chars<default_container_type>(lit.raw);
            if (!num.has_value()) return std::unexpected{machine_ec::unknown};
            return eval_result{integer{std::move(*num)}};
          },
          [](const token_real_literal& lit) -> std::expected<eval_result, machine_ec> {
            auto parsed = parse_real_literal(lit.raw);
            if (!parsed.has_value()) return std::unexpected{machine_ec::unknown};
            return eval_result{epx::make_q(std::move(parsed->num), std::move(parsed->den))};
          },
          [this](const token_id& id) -> std::expected<eval_result, machine_ec> {
            if (auto res = ctx_->lookup(std::string{id.raw}); res.has_value()) {
              return eval(*res);
            }
            auto name = std::string{id.raw};
            if (auto it = builtin_constants.find(name); it != builtin_constants.end()) {
              return eval_result{it->second()};
            }
            return std::unexpected{machine_ec::unknown};
          },
          [](auto) -> std::expected<eval_result, machine_ec> { return std::unexpected{machine_ec::unknown}; },
      },
      v->val);
}

std::expected<eval_result, machine_ec> eval_session::eval_paren(const details::paren_expr* p) const {
  return eval(p->inner);
}

std::expected<eval_result, machine_ec> eval_session::eval_binop(const details::binop_expr* b) const {
  auto lhs = eval(b->left);
  if (!lhs.has_value()) return std::unexpected{lhs.error()};
  auto rhs = eval(b->right);
  if (!rhs.has_value()) return std::unexpected{rhs.error()};

  auto lhs_is_real = std::holds_alternative<real>(*lhs);
  auto rhs_is_real = std::holds_alternative<real>(*rhs);
  auto kind = b->kind();

  // Integer division: always produce real.
  if (kind == node_kind::div_expr && !lhs_is_real && !rhs_is_real) {
    auto& li = std::get<integer>(*lhs);
    auto& ri = std::get<integer>(*rhs);
    return eval_result{epx::make_q(li, ri)};
  }

  // One or both operands are real: promote both to real.
  if (lhs_is_real || rhs_is_real) {
    auto lr = lhs_is_real ? std::move(std::get<real>(*lhs))
                          : epx::make_q(std::get<integer>(*lhs), epx::details::one<default_container_type>());
    auto rr = rhs_is_real ? std::move(std::get<real>(*rhs))
                          : epx::make_q(std::get<integer>(*rhs), epx::details::one<default_container_type>());
    switch (kind) {
      case node_kind::add_expr:
        return eval_result{epx::add(std::move(lr), std::move(rr))};
      case node_kind::sub_expr:
        return eval_result{epx::add(std::move(lr), epx::opp(std::move(rr)))};
      case node_kind::mul_expr:
        return eval_result{epx::mul(std::move(lr), std::move(rr))};
      case node_kind::div_expr:
        return eval_result{epx::mul(std::move(lr), epx::inv(std::move(rr)))};
      default:
        return std::unexpected{machine_ec::unknown};
    }
  }

  // Both integers.
  auto& li = std::get<integer>(*lhs);
  auto& ri = std::get<integer>(*rhs);
  switch (kind) {
    case node_kind::add_expr:
      return eval_result{epx::add(li, ri)};
    case node_kind::sub_expr:
      return eval_result{epx::sub(li, ri)};
    case node_kind::mul_expr:
      return eval_result{epx::mul(li, ri)};
    default:
      return std::unexpected{machine_ec::unknown};
  }
}

std::expected<eval_result, machine_ec> eval_session::eval_funcall(const details::func_call* f) const {
  std::vector<eval_result> args;
  args.reserve(f->params.size());
  for (auto* p : f->params) {
    auto arg = eval(p);
    if (!arg.has_value()) return std::unexpected{arg.error()};
    args.push_back(std::move(*arg));
  }

  auto name = std::string{f->id.raw};

  if (auto it = unary_funcs.find(name); it != unary_funcs.end()) {
    if (args.size() != 1) return std::unexpected{machine_ec::unknown};
    auto r = std::holds_alternative<real>(args[0])
                 ? std::move(std::get<real>(args[0]))
                 : epx::make_q(std::get<integer>(args[0]), epx::details::one<default_container_type>());
    return eval_result{it->second(std::move(r))};
  }

  if (auto it = binary_funcs.find(name); it != binary_funcs.end()) {
    if (args.size() != 2) return std::unexpected{machine_ec::unknown};
    auto r0 = std::holds_alternative<real>(args[0])
                  ? std::move(std::get<real>(args[0]))
                  : epx::make_q(std::get<integer>(args[0]), epx::details::one<default_container_type>());
    auto r1 = std::holds_alternative<real>(args[1])
                  ? std::move(std::get<real>(args[1]))
                  : epx::make_q(std::get<integer>(args[1]), epx::details::one<default_container_type>());
    return eval_result{it->second(std::move(r0), std::move(r1))};
  }

  return std::unexpected{machine_ec::unknown};
}

}  // namespace

std::optional<machine_ec> eval_context::define(std::string name, expr* e) noexcept {
  if (auto res = lookup(name); res.has_value()) {
    return machine_ec::unknown;
  } else {
    table_.emplace(std::move(name), e);
    return std::nullopt;
  }
}

std::expected<expr*, machine_ec> eval_context::lookup(const std::string& name) const noexcept {
  if (auto iter = table_.find(name); iter != table_.end()) {
    return iter->second;
  } else if (parent_) {
    return parent_->lookup(name);
  }
  return std::unexpected{machine_ec::unknown};
}

std::expected<std::vector<eval_result>, machine_ec> machine::execute(const mathscript& script) noexcept {
  eval_context ctx;
  eval_session session{&ctx};
  std::vector<eval_result> results;

  for (auto& s : script.statements) {
    auto* e = std::get<expr*>(s);
    auto result = session.eval(e);
    if (!result.has_value()) return std::unexpected{result.error()};
    results.push_back(std::move(*result));
  }

  return results;
}

}  // namespace epx::script
