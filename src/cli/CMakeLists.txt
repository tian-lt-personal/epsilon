add_executable(cli
  main.cpp
)

target_compile_features(cli PRIVATE cxx_std_23)
target_link_libraries(cli PRIVATE epsilon_engine)
target_compile_definitions(cli PRIVATE
  EPSILON_VERSION="${PROJECT_VERSION}"
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
)

if(MSVC)
  target_compile_options(cli PRIVATE /W4 /WX)
else()
  target_compile_options(cli PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
