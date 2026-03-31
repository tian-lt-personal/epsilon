cmake_minimum_required(VERSION 3.20)

project(epsilon VERSION 0.0.1 LANGUAGES CXX)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

add_subdirectory(epsilon)
add_subdirectory(engine)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(ut)
endif()
