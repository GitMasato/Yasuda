if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree")
endif()

cmake_minimum_required(VERSION 3.10)
set(CMAKE_VERBOSE_MAKEFILE True)

project(emgsim LANGUAGES CXX)
include_directories(./include)
include_directories(./include/ThirdParty/clipp/include)
file(GLOB EMGSRC ./src/*.cpp)
add_executable(emgsim ${EMGSRC})

find_package(OpenMP REQUIRED)
target_link_libraries(emgsim PRIVATE stdc++fs OpenMP::OpenMP_CXX)
target_compile_options(emgsim PRIVATE -O3 -Wall -Wextra -Winline)
target_compile_features(emgsim PRIVATE cxx_std_17)
