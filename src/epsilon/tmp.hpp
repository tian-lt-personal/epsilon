// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_INC_TMP_HPP
#define EPSILON_INC_TMP_HPP

namespace epx::tmp {

template <class... Fs>
struct overloads: Fs... {
  using Fs::operator()...;
};

} // namespace epx::tmp

#endif  // EPSILON_INC_TMP_HPP
