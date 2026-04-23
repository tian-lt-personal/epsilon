// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#include "machine.hpp"

namespace epx::script {

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

}  // namespace epx::script
