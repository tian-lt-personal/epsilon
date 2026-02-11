include(FetchContent)

# The policy uses the download time for timestamp, instead of the timestamp in the archive. This
# allows for proper rebuilds when a projects URL changes.
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
  URL_HASH SHA256=65fab701d9829d38cb77c14acdc431d2108bfdbf8979e40eb8ae567edf10b27c
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE) # For Windows: Prevent overriding the parent project's compiler/linker settings
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE) # Disable GMock
FetchContent_MakeAvailable(googletest)

add_executable(epsilon_ut
  coro_tests.cpp
  n_tests.cpp
  ops_tests.cpp
  z_tests.cpp
)
target_link_libraries(epsilon_ut PRIVATE
    gtest_main
    epsilon
)

if(MSVC)
  target_compile_options(epsilon_ut PRIVATE /W4 /WX)
else()
  target_compile_options(epsilon_ut PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

include(GoogleTest)
gtest_discover_tests(epsilon_ut)
