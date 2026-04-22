add_library(epsilon_engine
  lexer.cpp
  parser.cpp
  utils.cpp
)

target_compile_features(epsilon_engine PUBLIC cxx_std_23)
target_include_directories(epsilon_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(epsilon_engine PUBLIC
  epsilon
)

if(MSVC)
  target_compile_options(epsilon_engine PRIVATE /W4 /WX)
else()
  target_compile_options(epsilon_engine PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
