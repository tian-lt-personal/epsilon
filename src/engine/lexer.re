// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#include "lexer.hpp"

namespace epx {

std::expected<token, token_error> lexer::operator()() {
  const char*& cursor = cursor_;
  const char* const limit = input_.data() + input_.length();
  const char* marker = nullptr;

lex_again:
  marker = cursor;

/*!re2c
   re2c:api:style          = free-form;
   re2c:indent:string      = "  ";
   re2c:define:YYCTYPE     = "char";
   re2c:define:YYCURSOR    = "cursor";
   re2c:define:YYLIMIT     = "limit";
   re2c:define:YYFILL      = "return std::unexpected{token_error{.code = token_ec::eof, .line = line_, .column = column_}};";

   [ \t\r\n]+ {
     for (auto p = marker; p != cursor; ++p) {
        if (*p == '\n') { ++line_; column_ = 1; }
        else ++column_;
     }
     goto lex_again;
   }

   "+" { ++column_; return token_op_add{}; }

   * { return std::unexpected{token_error{.code = token_ec::bad_input, .line = line_, .column = column_}}; }
*/
}

} // namespace epx
