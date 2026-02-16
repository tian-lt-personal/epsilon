// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#pragma once

// std
#include <cstdint>
#include <string_view>
#include <vector>

// epx
#include "chars.hpp"

namespace epxut {

using sz = epx::z<std::vector<uint8_t>>;
using mz = epx::z<std::vector<uint16_t>>;
using lz = epx::z<std::vector<uint32_t>>;

constexpr sz stosz(std::string_view chars) { return epx::try_from_chars<sz::container_type>(chars).value(); }

}  // namespace epxut
