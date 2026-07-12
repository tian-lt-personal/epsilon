// SPDX-License-Identifier: MIT
// Copyright (c) 2026-present Tian Liao

// std
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

// epx
#include "chars.hpp"
#include "machine.hpp"
#include "tmp.hpp"

namespace {

namespace script = epx::script;

struct config {
  unsigned int precision = 50;
  bool strip_zeros = true;
  std::string config_path;
};

struct mode_inline {
  std::string expression;
};

struct mode_file {
  std::string path;
};

struct args {
  config cfg;
  std::optional<std::variant<mode_inline, mode_file>> mode;
};

std::string default_config_path() {
  const char* home = std::getenv("HOME");
  if (home == nullptr) {
    home = std::getenv("USERPROFILE");
  }
  std::string base = home != nullptr ? std::string{home} : std::string{"."};
  return base + "/.epsilon/config";
}

std::string_view trim(std::string_view s) {
  auto is_space = [](char c) { return std::isspace(static_cast<unsigned char>(c)); };

  auto no_leading = s | std::views::drop_while(is_space);

  auto trailing = static_cast<size_t>(std::ranges::distance(
      no_leading | std::views::reverse | std::views::take_while(is_space)));

  auto total = static_cast<size_t>(std::ranges::distance(no_leading));
  auto len = total - trailing;

  if (len == 0) return {};
  return std::string_view{&*no_leading.begin(), len};
}

config load_config(const std::string& path) {
  config cfg;
  cfg.config_path = path;

  std::ifstream file(path);
  if (!file.good()) return cfg;
  file.exceptions(std::ios::failbit | std::ios::badbit);

  std::string line;
  try {
    while (std::getline(file, line)) {
      if (!line.empty() && line.back() == '\r') {
        line.pop_back();
      }

      auto trimmed = trim(line);
      if (trimmed.empty() || trimmed.starts_with('#')) continue;

      auto eq = trimmed.find('=');
      if (eq == std::string_view::npos) continue;

      auto key = trim(trimmed.substr(0, eq));
      auto val = trim(trimmed.substr(eq + 1));

      if (key == "precision") {
        try {
          auto n = std::stoul(std::string{val});
          cfg.precision = static_cast<unsigned int>(n);
        } catch (...) {
          // Malformed value -- keep the default
        }
      } else if (key == "strip_zeros") {
        cfg.strip_zeros = val == "true" || val == "1" || val == "yes";
      }
    }
  } catch (const std::ios_base::failure&) {
    if (!file.eof()) throw;
  }

  return cfg;
}

std::optional<args> parse_args(int argc, char* argv[]) {
  std::optional<unsigned int> precision_override;
  std::optional<std::string> config_path_override;
  std::optional<bool> strip_zeros_override;
  std::optional<std::string> positional;
  bool end_of_options = false;

  for (int i = 1; i < argc; ++i) {
    std::string_view arg = argv[i];

    if (!end_of_options && arg == "--") {
      end_of_options = true;
      continue;
    }

    if (!end_of_options && arg == "-p") {
      if (i + 1 >= argc) return std::nullopt;
      std::string_view val = argv[++i];
      try {
        auto n = std::stoul(std::string{val});
        precision_override = static_cast<unsigned int>(n);
      } catch (...) {
        return std::nullopt;
      }
    } else if (!end_of_options && arg.starts_with("-p") && arg.size() > 2) {
      try {
        auto n = std::stoul(std::string{arg.substr(2)});
        precision_override = static_cast<unsigned int>(n);
      } catch (...) {
        return std::nullopt;
      }
    } else if (!end_of_options && arg == "--precision") {
      if (i + 1 >= argc) return std::nullopt;
      std::string_view val = argv[++i];
      try {
        auto n = std::stoul(std::string{val});
        precision_override = static_cast<unsigned int>(n);
      } catch (...) {
        return std::nullopt;
      }
    } else if (!end_of_options && arg.starts_with("--precision=")) {
      try {
        auto n = std::stoul(std::string{arg.substr(12)});
        precision_override = static_cast<unsigned int>(n);
      } catch (...) {
        return std::nullopt;
      }
    } else if (!end_of_options && arg == "--config-path") {
      if (i + 1 >= argc) return std::nullopt;
      config_path_override = std::string{argv[++i]};
    } else if (!end_of_options && arg == "--strip-zeros") {
      strip_zeros_override = true;
    } else if (!end_of_options && arg == "--no-strip-zeros") {
      strip_zeros_override = false;
    } else if (!end_of_options && arg.starts_with('-')) {
      return std::nullopt;
    } else {
      if (positional.has_value()) return std::nullopt;
      positional = std::string{arg};
    }
  }

  // Resolve and load config
  args result;
  result.cfg = load_config(config_path_override.value_or(default_config_path()));

  // CLI flag overrides config file
  if (precision_override.has_value()) {
    result.cfg.precision = *precision_override;
  }
  if (strip_zeros_override.has_value()) {
    result.cfg.strip_zeros = *strip_zeros_override;
  }

  // Determine mode from positional argument
  if (!positional.has_value()) {
    result.mode = std::nullopt;
  } else {
    std::ifstream test_file(*positional);
    if (test_file.good()) {
      result.mode = mode_file{*positional};
    } else {
      result.mode = mode_inline{*positional};
    }
  }

  return result;
}

void strip_trailing_zeros(std::string& s) {
  auto dot = s.find('.');
  if (dot == std::string::npos) return;

  while (s.back() == '0') s.pop_back();
  if (s.back() == '.') s.pop_back();
}

std::string evaluate(std::string_view expr, unsigned int precision, bool strip_zeros) {
  auto script_res = script::translate(expr);
  if (!script_res.has_value()) return {};

  script::machine m;
  auto exec_res = m.execute(*script_res);
  if (!exec_res.has_value() || exec_res->empty()) return {};

  auto& result = (*exec_res)[0];

  if (std::holds_alternative<script::real>(result)) {
    auto s = epx::to_string(std::move(std::get<script::real>(result)), precision);
    if (strip_zeros) strip_trailing_zeros(s);
    return s;
  }
  if (std::holds_alternative<script::integer>(result)) {
    return epx::to_string(std::get<script::integer>(result));
  }
  return {};
}

std::string evaluate_safe(std::string_view expr, unsigned int precision, bool strip_zeros) {
  try {
    auto result = evaluate(expr, precision, strip_zeros);
    if (result.empty()) {
      std::cerr << "error: could not evaluate expression\n";
    }
    return result;
  } catch (const std::exception& ex) {
    std::cerr << "error: " << ex.what() << '\n';
    return {};
  } catch (...) {
    std::cerr << "error: unknown runtime error\n";
    return {};
  }
}

int process_file(const std::string& path, unsigned int precision, bool strip_zeros) {
  std::ifstream file(path);
  if (!file.good()) {
    std::cerr << "error: cannot open file '" << path << "'\n";
    return EXIT_FAILURE;
  }
  file.exceptions(std::ios::failbit | std::ios::badbit);

  std::string line;
  try {
    while (std::getline(file, line)) {
      if (!line.empty() && line.back() == '\r') {
        line.pop_back();
      }

      auto trimmed = trim(line);
      if (trimmed.empty()) continue;
      if (trimmed.starts_with('#') || trimmed.starts_with("//")) continue;

      auto result = evaluate_safe(line, precision, strip_zeros);
      if (!result.empty()) {
        std::cout << result << '\n';
      }
    }
  } catch (const std::ios_base::failure&) {
    if (file.eof()) return EXIT_SUCCESS;
    std::cerr << "error: read failure in file '" << path << "'\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int run_repl(config& cfg) {
  std::cout << "Epsilon " << EPSILON_VERSION << " -- arbitrary-precision math\n";
  std::cout << "Type an expression or a command (.help for more).\n";

  std::string line;
  for (;;) {
    std::cout << "> " << std::flush;
    if (!std::getline(std::cin, line)) break;

    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }

    auto trimmed = trim(line);
    if (trimmed.empty()) continue;

    // Dot-prefixed REPL commands
    if (trimmed.starts_with('.')) {
      auto cmd = std::string{trimmed};

      if (cmd.starts_with(".precision ")) {
        auto val = trim(std::string_view{cmd}.substr(11));
        try {
          auto n = std::stoul(std::string{val});
          cfg.precision = static_cast<unsigned int>(n);
          std::cout << "precision = " << cfg.precision << '\n';
        } catch (...) {
          std::cerr << "error: invalid precision value\n";
        }
      } else if (cmd.starts_with(".strip_zeros")) {
        auto val = trim(std::string_view{cmd}.substr(12));
        if (val.empty()) {
          cfg.strip_zeros = !cfg.strip_zeros;
        } else {
          cfg.strip_zeros = val == "on" || val == "true" || val == "1";
        }
        std::cout << "strip_zeros = " << (cfg.strip_zeros ? "true" : "false") << '\n';
      } else if (cmd == ".config") {
        std::cout << "config path: " << cfg.config_path << '\n';
        std::cout << "precision:   " << cfg.precision << '\n';
        std::cout << "strip_zeros: " << (cfg.strip_zeros ? "true" : "false") << '\n';
      } else if (cmd == ".help") {
        std::cout << "Commands:\n";
        std::cout << "  .precision N     Set display precision (decimal places)\n";
        std::cout << "  .strip_zeros     Toggle trailing-zero stripping (on|off)\n";
        std::cout << "  .config          Show current configuration\n";
        std::cout << "  .help            Show this help\n";
        std::cout << "  exit, quit     Exit the REPL\n";
        std::cout << "\n";
        std::cout << "Expression syntax:\n";
        std::cout << "  Integers:  z'42\n";
        std::cout << "  Reals:     3.14, .5, r'1.5\n";
        std::cout << "  Operators: + - * /\n";
        std::cout << "  Functions:  sin, cos, tan, exp, log, ln, sqrt,\n";
        std::cout << "              arcsin, arccos, arctan, sinh, cosh, tanh,\n";
        std::cout << "              arcsinh, arccosh, arctanh, pow\n";
        std::cout << "  Constants:  pi, e\n";
        std::cout << "  Examples:   1 + 2 * 3, sin(pi/2), sqrt(2), pow(2, 10)\n";
      } else {
        std::cerr << "error: unknown command '" << cmd << "' (try .help)\n";
      }
      continue;
    }

    // Convenience aliases (non-dot-prefixed)
    if (trimmed == "exit" || trimmed == "quit") break;
    if (trimmed == "help") {
      std::cout << "Type .help for command reference.\n";
      continue;
    }

    auto result = evaluate_safe(line, cfg.precision, cfg.strip_zeros);
    if (!result.empty()) {
      std::cout << result << '\n';
    }
  }

  std::cout << '\n';
  return EXIT_SUCCESS;
}

void print_usage(std::string_view prog_name) {
  std::cerr << "Usage: " << prog_name << " [OPTIONS] [EXPRESSION | FILE]\n";
  std::cerr << "\n";
  std::cerr << "Options:\n";
  std::cerr << "  -p, --precision N   Decimal places for real numbers (default: 50)\n";
  std::cerr << "  --strip-zeros       Strip trailing zeros from real output (default)\n";
  std::cerr << "  --no-strip-zeros    Keep trailing zeros in real output\n";
  std::cerr << "  --config-path PATH  Path to config file (default: ~/.epsilon/config)\n";
  std::cerr << "\n";
  std::cerr << "Modes:\n";
  std::cerr << "  (no positional arg) Interactive REPL session\n";
  std::cerr << "  EXPRESSION           Evaluate and print result\n";
  std::cerr << "  FILE                 Evaluate each non-empty line of the file\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  auto parsed = parse_args(argc, argv);
  if (!parsed.has_value()) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  auto& a = *parsed;

  if (!a.mode.has_value()) {
    return run_repl(a.cfg);
  }

  return std::visit(
      epx::tmp::overloads{
          [&](const mode_inline& m) {
            auto result = evaluate_safe(m.expression, a.cfg.precision, a.cfg.strip_zeros);
            if (result.empty()) return EXIT_FAILURE;
            std::cout << result << '\n';
            return EXIT_SUCCESS;
          },
          [&](const mode_file& m) {
            return process_file(m.path, a.cfg.precision, a.cfg.strip_zeros);
          },
      },
      *a.mode);
}
