Experimental project - work in progress

---

![](./docs/imgs/epsilon.svg)

[![Epsilon-CI](https://github.com/tian-lt-personal/epsilon/actions/workflows/epsilon-ci.yml/badge.svg)](https://github.com/tian-lt-personal/epsilon/actions/workflows/epsilon-ci.yml)

An arbitrary precision arithmetic library for **computable real numbers**.

## Highlights

- **Lazy evaluation**: Dynamically adjusts precision as needed, enabling efficient computations and on-demand accuracy.
- **Supported domains**: 
	- **Z**: Big integer numbers.
	- **R**: Computable real numbers, supporting operations with dynamic precision.
- **Algorithmic foundation**: Based on algorithms by V. Menissier-Morain and Donald Knuth, ensuring mathematical rigor and reliability.
- **Modern C++**: Written in modern C++ (C++23), leveraging concepts, coroutines, and advanced type traits for performance and clarity.
- **Cross-platform**: The library can be compiled by clang, g++, and cl (MSVC).

## Features

- Efficient, on-demand computation of real numbers with adjustable precision.
- Operator overloading for natural arithmetic expressions.
- Conversion utilities for parsing and formatting numbers.
- Extensible design for custom digit and container types.

## CLI

A command-line calculator (`cli`) is included. It supports three modes:

**Immediate mode** — inline expression:

```bash
cli "1 + 2 * 3"
cli -p 10 "sin(pi/2)"
```

**File mode** — evaluate each non-empty line of a file (`#` and `//` comments are skipped):

```bash
cli script.math
```

**REPL mode** — interactive session (no arguments):

```bash
cli
```

| Option | Description |
|--------|-------------|
| `-V`, `--version` | Print version and exit |
| `-p N`, `--precision N` | Decimal places for real numbers (default: 50) |
| `--strip-zeros`, `--no-strip-zeros` | Strip / keep trailing zeros |
| `--config-path PATH` | Path to config file (default: `~/.epsilon/config`) |

REPL commands: `.precision N`, `.strip_zeros [on|off]`, `.config`, `.help`, `exit`, `quit`.

Config file (`~/.epsilon/config`): `precision = 50`, `strip_zeros = true`.

Supported expression syntax: real literals (`3.14`), integer literals (`z'42`), operators (`+ - * /`), functions (`sin`, `cos`, `tan`, `exp`, `log`, `ln`, `sqrt`, `arcsin`, `arccos`, `arctan`, `sinh`, `cosh`, `tanh`, `arcsinh`, `arccosh`, `arctanh`, `pow`), and constants (`pi`, `e`).

## References

- V. Menissier-Morain, "Arbitrary Precision Real Arithmetic: design and algorithms"
- Donald Knuth, "The Art of Computer Programming"

### TODOs
- Math expression - parser ✅
- Math expression - evaluation ✅
- CLI application ✅
- Lexer: underscore support in identifiers
- Parser: unary minus
- Parser: multiple statements
- Parser: assignment (`=`)
- Evaluator: user-defined functions with closures
- Evaluator: richer error diagnostics
