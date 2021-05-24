cmake_minimum_required(VERSION 3.10)
set(CMAKE_VERBOSE_MAKEFILE True)
project(DEM LANGUAGES CXX)


set(src
    main.cpp
    Generator.cpp
)

add_executable(GEN ${dem_src})
target_link_libraries(DEM PRIVATE stdc++fs OpenMP::OpenMP_CXX)
target_compile_options(DEM PRIVATE -O3 -Wall -Wextra -Winline)
target_compile_features(DEM PRIVATE cxx_std_17)
