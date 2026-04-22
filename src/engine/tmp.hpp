// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#ifndef EPSILON_ENGINE_INC_TMP_HPP
#define EPSILON_ENGINE_INC_TMP_HPP

#if defined(_MSC_VER) && !defined(__clang__)
#define _EPX_EMPTY_BASE __declspec(empty_bases)
#else
#define _EPX_EMPTY_BASE
#endif

namespace epx::tmp {

template <class... Ts>
struct _EPX_EMPTY_BASE overloads : Ts... {
  using Ts::operator()...;
};

}  // namespace epx::tmp

#endif  // EPSILON_ENGINE_INC_TMP_HPP
