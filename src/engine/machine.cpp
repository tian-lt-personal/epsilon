// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// epx
#include "tmp.hpp"

//  engine
#include "machine.hpp"

namespace epx::script {

namespace {

class eval_session {
 public:
  std::expected<eval_result, machine_ec> eval(const expr*) && {}

 private:
  eval_result result_;
};

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
  std::vector<eval_result> results;
  for (auto& s : script.statements) {
    (void)s;
    std::visit(
        tmp::overloads{
            //[](const expr*) { std::terminate(); },
            [](auto) { std::terminate(); },
        },
        s);
  }
  abort();
}

}  // namespace epx::script
