cmake_minimum_required(VERSION 3.10)
set(CMAKE_VERBOSE_MAKEFILE True)
project(PtclGen LANGUAGES CXX)

set(src
    main.cpp
    Generator.cpp
)

add_executable(PtclGen ${src})
target_link_libraries(PtclGen PRIVATE stdc++fs)
target_compile_options(PtclGen PRIVATE -O3 -Wall -Wextra -Winline)
target_compile_features(PtclGen PRIVATE cxx_std_17)

