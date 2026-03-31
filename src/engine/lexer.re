// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

#include "lexer.hpp"

namespace epx {

std::expected<token, token_error> lexer::operator()() noexcept {
  const char*& cursor = cursor_;
  const char* const limit = input_.data() + input_.length();
  const char* marker = nullptr;

lex_again:
  const char* start = cursor;

  if (cursor >= limit) {
    return std::unexpected{token_error{.code = token_ec::eof, .line = line_, .column = column_}};
  }

/*!re2c
   re2c:api:style          = free-form;
   re2c:yyfill:enable      = 0;
   re2c:indent:string      = "  ";
   re2c:define:YYPEEK      = "cursor < limit ? *cursor : 0";
   re2c:define:YYCTYPE     = "char";
   re2c:define:YYCURSOR    = "cursor";
   re2c:define:YYMARKER    = "marker";
   re2c:define:YYLIMIT     = "limit";

   // decimal value
   [+-]? ( ([0-9]+ ("." [0-9]*)?) | ("." [0-9]+) ) {
     column_ += (size_t)(cursor - start);
     return token_val_decimal {
       .raw = std::string_view(marker, (size_t)(cursor - start)),
     };
   }

   // whitespaces
   [ \t\r\n]+ {
     for (auto p = marker; p != cursor; ++p) {
        if (*p == '\n') { ++line_; column_ = 1; }
        else ++column_;
     }
     goto lex_again;
   }

   // operators and parenthesises
   "*" { ++column_; return token_op_mul{}; }
   "/" { ++column_; return token_op_div{}; }
   "%" { ++column_; return token_op_percent{}; }
   "(" { ++column_; return token_lparen{}; }
   ")" { ++column_; return token_rparen{}; }
   "+" { ++column_; return token_op_add{}; }
   "-" { ++column_; return token_op_sub{}; }

   // wildcard
   * { return std::unexpected{token_error{.code = token_ec::bad_input, .line = line_, .column = column_}}; }
*/
}

} // namespace epx
