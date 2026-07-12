// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_ENGINE_INC_MACHINE_HPP
#define EPSILON_ENGINE_INC_MACHINE_HPP

// std
#include <expected>
#include <map>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

// epx
#include "parser.hpp"
#include "r.hpp"
#include "z.hpp"

namespace epx::script {

enum struct machine_ec {
  unknown,  // todo: improve diagnostics
};

using real = r<default_container_type>;
using integer = z<default_container_type>;

class eval_context {
 public:
  eval_context() = default;
  explicit eval_context(eval_context* parent) : parent_(parent) {}
  std::optional<machine_ec> define(std::string name, expr* e) noexcept;
  std::expected<expr*, machine_ec> lookup(const std::string& name) const noexcept;

 private:
  std::map<std::string, expr*> table_;
  const eval_context* parent_ = nullptr;
};

struct function {
  std::vector<token_id> var_seq;
  expr* body;
  std::shared_ptr<eval_context> context;
};

using eval_result = std::variant<std::monostate, real, integer, function>;

class machine {
 public:
  std::expected<std::vector<eval_result>, machine_ec> execute(const mathscript& script) noexcept;
};

}  // namespace epx::script

#endif  // EPSILON_ENGINE_INC_MACHINE_HPP
