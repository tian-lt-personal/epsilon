cmake_minimum_required(VERSION 3.20)

if(DEFINED EPX_VERSION)
  set(_epx_ver "${EPX_VERSION}")
else()
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _epx_ver)
  string(STRIP "${_epx_ver}" _epx_ver)
endif()

project(epsilon VERSION "${_epx_ver}" LANGUAGES CXX)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

add_subdirectory(epsilon)
add_subdirectory(engine)
add_subdirectory(cli)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(ut)
endif()
