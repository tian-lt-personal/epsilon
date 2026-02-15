// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_INC_TMP_HPP
#define EPSILON_INC_TMP_HPP

#ifdef _MSC_VER
#define EPSILON_DECL_EMPTY_BASES __declspec(empty_bases)
#else
#define EPSILON_DECL_EMPTY_BASES
#endif  // _MSC_VER

namespace epx::tmp {

template <class... Fs>
struct EPSILON_DECL_EMPTY_BASES overloads : Fs... {
  using Fs::operator()...;
};

}  // namespace epx::tmp

#endif  // EPSILON_INC_TMP_HPP
