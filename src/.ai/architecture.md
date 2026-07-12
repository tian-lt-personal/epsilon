# Epsilon — Architecture Summary

## Overview

Epsilon is a **C++23 arbitrary-precision arithmetic library** combined with a **math scripting engine**. It provides exact big-integer arithmetic (`z`), lazy real-number arithmetic (`r`) with arbitrary-precision elementary functions, and a parser/lexer for math expressions. The project is authored by Tian Liao (MIT license).

---

## Project Structure

```
src/
├── CMakeLists.txt              # Top-level CMake (project root)
├── .clang-format               # Google-based style, 120 column limit
├── epsilon/                    # Core math library (header-only, INTERFACE)
│   ├── CMakeLists.txt
│   ├── t.hpp                   # Core types, concepts, digit aliases
│   ├── z.hpp                   # Big integer arithmetic
│   ├── r.hpp                   # Lazy real numbers & elementary functions
│   ├── coro.hpp                # Lazy coroutine for R approximations
│   ├── ops.hpp                 # C++ operator overloads for z<C>
│   ├── chars.hpp               # String conversion (parse/format)
│   └── tmp.hpp                 # TMP utilities (overloads pattern)
├── engine/                     # Math expression parsing (static lib)
│   ├── CMakeLists.txt
│   ├── lexer.hpp / lexer.cpp   # Lexer: text → tokens
│   ├── parser.hpp / parser.cpp # Recursive-descent parser: tokens → AST
│   ├── utils.hpp / utils.cpp   # AST dump / pretty-print
│   └── tmp.hpp                 # TMP utilities (duplicate)
├── ut/                         # Unit tests (GoogleTest)
│   ├── CMakeLists.txt          # FetchContent for GTest 1.17.0
│   ├── def.hpp                 # Shared test type aliases
│   ├── z_tests.cpp             # Big integer tests
│   ├── r_tests.cpp             # Real number tests
│   ├── coro_tests.cpp          # Coroutine tests
│   ├── chars_tests.cpp         # String conversion tests
│   ├── ops_tests.cpp           # Operator overload tests
│   ├── n_tests.cpp             # (TBD)
│   ├── lexer_tests.cpp         # Lexer tests
│   └── parser_tests.cpp        # Parser/AST tests
└── b/                          # CMake build output (gitignored)
```

---

## Tech Stack

| Layer        | Technology                         |
|--------------|------------------------------------|
| Language     | C++23 (`cxx_std_23`)               |
| Build system | CMake ≥ 3.20                       |
| Test framework | GoogleTest 1.17.0 (FetchContent) |
| Code style   | clang-format, Google-based, 120 cols |
| Compilers    | MSVC (with `/W4 /WX`) and Clang   |
| CI           | GitHub Actions (per git history)   |

---

## Module 1: `epsilon` — Core Math Library

**Type**: Header-only INTERFACE library (`add_library(epsilon INTERFACE)`). No compiled sources. No external dependencies beyond the C++ standard library.

### 1a. Core Types (`t.hpp`)

- **`z<C>`**: Big integer template parameterized by container `C`. Fields: `digits` (LSD-first unsigned digit vector), `sgn` (sign enum: positive/negative). Supports `<=>` comparison.
- **`sign`**: `enum class` with `positive` and `negative`.
- **`container` concept**: Constrains `C` to be a random-access, sized range of unsigned integers supporting `push_back`, `reserve`, and digit size < max digit size.
- **Digit type aliases**:
  - `default_digit_type` = `uint32_t`
  - `max_digit_type` = `uint32_t`
  - `default_container_type` = `std::vector<uint32_t>`
  - `wide_digit_type<D>` — maps `uint8_t→uint16_t`, `uint16_t→uint32_t`, `uint32_t→uint64_t`
- **Error types**: `divide_by_zero_error`, `msd_overflow_error`, `precision_overflow_error`, `negative_radicand_error`, `negative_zpow_error`, `kthroot_too_small_error`, `non_positive_log_error` — all inherit `std::runtime_error`.

### 1b. Big Integer Arithmetic (`z.hpp`)

All functions are `constexpr` and templated on container `C`. Namespace: `epx`.

| Function      | Description |
|---------------|-------------|
| `create(v)`   | Convert integral `T` to `z<C>` |
| `is_zero(n)`  | Check if magnitude is zero |
| `normalize(n)`| Strip leading-zero digits |
| `negate(n)`   | Flip sign in place |
| `cmp_n(a,b)`  | Unsigned magnitude comparison |
| `add_n(a,b)`  | Unsigned addition with carry chain |
| `sub_n(a,b)`  | Unsigned subtraction (requires |a| ≥ |b|) |
| `mul_n(a,b)`  | Unsigned schoolbook multiplication |
| `div_n(u,v)`  | Unsigned division (Knuth Algorithm D for multi-digit; single-digit fast path) |
| `add/sub/mul` | Signed wrappers calling the unsigned `_n` variants |
| `div(a,b)`    | Signed division → `{q, r}` with truncation toward zero |
| `floor_div(a,b)` | Floor division (remainder has same sign as divisor) |
| `ceil_div(a,b)`  | Ceil division (remainder has opposite sign of divisor) |
| `mul_2exp/4exp`  | Left/right bit-shift by power of 2 or 4 |
| `pow(n, exp)` | Integer exponentiation |
| `root(n, k)`  | Integer k-th root (Newton's method) |
| `bit_length`  | Number of significant bits |

Internal helpers in `details` namespace: `zero()`, `one()`, `base()`, `ten()`, `umul()`, `div_2d()`, `bit_shift()`, `pow10()`.

### 1c. Lazy Real Numbers (`r.hpp`)

**`r<C>`**: Represents an exact real number as a **lazy function** `int → coro::lazy<z<C>>`. Given precision `n`, it returns `floor(x * 4^n)` — the value scaled by `4^n` and floored. This follows the Ménissier-Morain exact-real-arithmetic model with **base B=4**.

Key characteristics:
- **Memoization**: Once computed at precision `n`, results are cached; higher precision reuses lower via `mul_4exp`.
- **Constructors**: `make_q(p, q)` creates a rational `p/q` as a real.
- **Arithmetic**: `add`, `opp`, `mul`, `inv` — each returns a new `r<C>` that composes lazily.
- **Elementary functions** (all return `r<C>`):
  - `exp(x)` — via exp_rational decomposition + Taylor series
  - `log(x)` — via log_rational + arctanh-based decomposition
  - `sin(x)`, `cos(x)` — via argument reduction to `[0, π/2]` + Taylor series
  - `arctan(x)` — via argument reduction + direct series
  - `tan(x) = sin/cos`, `arcsin`, `arccos`
  - `sinh/cosh/tanh` — via `exp`
  - `arcsinh/arccosh/arctanh` — via `log` and `sqrt`
  - `pow(x, y) = exp(y*log(x))`, `log_base(x, b)`
  - `root(x, k)` — k-th root

Internal implementation details (all in `details` namespace, AI-generated following the paper's algorithms):
- `exp_taylor`, `compute_e`, `fp_mul`, `fp_pow`, `exp_rational`, `compute_d`, `compute_k` — for exponential
- `log_near_one_plus/minus`, `compute_ln2`, `log_rational` — for logarithm
- `atan_reciprocal`, `compute_pi`, `atan_series`, `atan_rational` — for arctangent
- `sin_series` — for sine
- All guarded with extra precision (`exp_guard=12`, `log_guard=12`, `atan_guard=12`, `sin_guard=12`)

### 1d. Lazy Coroutine (`coro.hpp`)

**`coro::lazy<T>`**: A custom C++20 coroutine type that suspends at `initial_suspend` and resumes when `co_await`ed. On first await, runs to completion, caches the result (or exception) in `lazy_promise::state_` (a `variant<monostate, T, exception_ptr>`), then resumes the awaiting coroutine via `final_suspend`.

Also provides `coro::forget` — a fire-and-forget coroutine that starts immediately (`suspend_never`) and terminates on exception.

### 1e. Operators (`ops.hpp`)

C++ operator overloads for `z<C>`: `+`, `-`, `*`, `/`, `%`. All delegate to the corresponding `epx::add/sub/mul/div` functions.

### 1f. String Conversion (`chars.hpp`)

- **`try_from_chars<C, B=10>(string_view) → optional<z<C>>`**: Base-10 string parsing with optional `+`/`-` sign. Returns `nullopt` on invalid input.
- **`to_string<C, B=10>(z<C>) → string`**: Repeated division by 10, reversed.
- **`to_string<C, B=10>(r<C>, k) → string`**: Real number to string with `k` decimal places, with proper rounding.

---

## Module 2: `engine` — Math Script Parser

**Type**: Static library (`epsilon_engine`). Links against `epsilon`. Namespace: `epx::script`.

### 2a. Lexer (`lexer.hpp`, `lexer.cpp`)

**`lexer` class**: Takes a `string_view`, returns `expected<token, token_ec>` on each call.

**Token types** (all in `token` variant):
- `token_integer_literal` — `z'123` (explicit integer)
- `token_real_literal` — `3.14` or `r'3.14` (explicit real)
- `token_id` — alphabetic identifiers like `f`, `sin`, `pi`
- Operators: `+`, `-`, `*`, `/`, `%`, `=`, `.`, `,`, `(`, `)`
- `std::monostate` — no token / drained
- Error: `token_ec::eof` or `token_ec::bad_input`

**Lexing rules**: Whitespace is skipped. `z'` prefix forces integer literal. `r'` prefix or leading digit/dot tries real literal. Letters start identifiers. Single-char operators match directly.

### 2b. Parser (`parser.hpp`, `parser.cpp`)

**Recursive descent parser** with operator precedence. Entry point: `translate(string_view) → expected<mathscript, translate_ec>`.

**AST node hierarchy** (all derive from `details::expr` with `node_kind` tag):
- `val_term` — literal or identifier (holds `variant<token_integer_literal, token_real_literal, token_id>`)
- `binop_expr` — binary operation (`add_expr`, `sub_expr`, `mul_expr`, `div_expr`)
- `paren_expr` — parenthesized sub-expression
- `func_call` — function call with parameter list

**Precedence**: `+`/`-` = 10, `*`/`/` = 20. Left-associative within each level.

**Memory management**: All AST nodes allocated from `std::pmr::monotonic_buffer_resource` (arena). Owned by `mathscript::ctx_`. AST nodes are raw pointers — lifetime is tied to the `mathscript` object.

**`translate()` flow**:
1. Create `lexer` + `tu` (translation unit)
2. `tu::initialize()` — prime first two tokens (lookahead of 1)
3. `tu::parse()` → `parse_expr()` → `parse_expr_with_precedence(min_prec)`:
   - Parse left term via `parse_term()`
   - While next token is binary op with precedence > min: consume op, parse right with higher precedence, build `binop_expr`
4. `parse_term()` dispatches on current token: integer/real/id → `val_term`; `(` → `paren_expr`; id followed by `(` → `func_call`

### 2c. Utils (`utils.hpp`, `utils.cpp`)

**`dump(script) → string`**: Pretty-prints the AST with indentation. Useful for debugging and tests. Handles `val`, `binop_expr`, `paren_expr`, `func_call`.

---

## Module 3: `ut` — Unit Tests

**Type**: Executable (`epsilon_ut`). Links `epsilon`, `epsilon_engine`, and `gtest_main`.

Uses `FetchContent` to pull GoogleTest 1.17.0. Tests discovered via `gtest_discover_tests()`.

Test files:
- `z_tests.cpp` — Integer create, add, sub, mul, div, floor_div, ceil_div, mul_4exp, pow, root (across `sz`, `mz`, `lz` digit sizes)
- `r_tests.cpp` — Real number arithmetic and elementary functions
- `coro_tests.cpp` — Lazy coroutine semantics
- `chars_tests.cpp` — String parse/format round-trips
- `ops_tests.cpp` — Operator overload behavior
- `lexer_tests.cpp` — Token recognition
- `parser_tests.cpp` — AST structure verification via `dump()` output matching

Test type aliases (from `def.hpp`):
- `sz` = `z<vector<uint32_t>>` (standard/default)
- `mz` = `z<vector<uint16_t>>` (medium)
- `lz` = `z<vector<uint64_t>>` (large) — not in `max_digit_type` constraint, used for boundary testing

---

## Dependency Graph

```
epsilon (INTERFACE, header-only)
   ↑
   ├── engine (STATIC lib, links epsilon)
   │
   └──┐
      ut (EXECUTABLE, links epsilon + engine + GTest)
```

No third-party runtime dependencies. GoogleTest is fetched at build time only.

---

## Key Design Decisions

1. **Header-only core**: The `epsilon` math library is header-only for easy integration. All functions are `constexpr`-compatible (though MSVC limits may apply).

2. **Base-4 real arithmetic**: The `r<C>` type uses base B=4 (power-of-two) for the exact real representation, following the Ménissier-Morain paper. This allows efficient bit-shift (`mul_4exp`) instead of general multiplication for scaling.

3. **Container-parameterized integers**: `z<C>` accepts any conforming container, enabling different digit sizes (16-bit, 32-bit, 64-bit) and custom allocators.

4. **Coroutine-based lazy evaluation**: Real number approximations use C++20 coroutines with memoization. This avoids recomputing lower-precision results and naturally composes for complex expressions.

5. **Arena-allocated AST**: The parser uses `pmr::monotonic_buffer_resource` for all AST nodes, giving fast allocation and single-point deallocation.

6. **AI-assisted implementation**: The elementary functions in `r.hpp` (exp, log, sin, cos, arctan, pi, and derived functions) were generated by AI following the algorithms from the Ménissier-Morain paper, as noted in code comments.

---

## Current State & Evolution

Recent commits show active development:
- **Core math** is mature: `z` arithmetic, `r` arithmetic, and all elementary functions are implemented.
- **Script engine** is in early stages: lexer, recursive-descent parser with binary ops, parenthesized expressions, identifiers, and function calls. No evaluator/semantic analysis yet.
- Latest commit (`func_call`) added function-call AST nodes to the parser.
