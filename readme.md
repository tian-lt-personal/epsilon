Experimental project - work in progress

> ### TODOs
> - (algo.7) Exponential function
> - (algo.8) Logarithm to base $B^\prime$
> - (algo.9) Arctangent
> - (algo.10) Sine function
> - Math expression - lexer/tokenizer
> - Math expression - parser 
> - Math expression - evaluation

---

![](./docs/imgs/epsilon.svg)

# epsilon

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

## References

- V. Menissier-Morain, "Arbitrary Precision Real Arithmetic: design and algorithms"
- Donald Knuth, "The Art of Computer Programming"
