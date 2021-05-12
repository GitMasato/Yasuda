cmake_minimum_required(VERSION 3.10)
set(CMAKE_VERBOSE_MAKEFILE True)
project(DEM LANGUAGES CXX)


set(dem_src 
    main.cpp
    Cell.cpp
    Collision.cpp
    DEM.cpp
    Field.cpp
    LongRange.cpp
    LongRangeParticleList.cpp
    Object.cpp
    OutputAVS.cpp
    OutputGeneral.cpp
    OutputVTK.cpp
    OutputVTK.cpp
    Particle.cpp
)

add_executable(DEM ${dem_src})

find_package(OpenMP REQUIRED)
target_link_libraries(DEM PRIVATE stdc++fs OpenMP::OpenMP_CXX)
target_compile_options(DEM PRIVATE -O3 -Wall -Wextra -Winline)
target_compile_features(DEM PRIVATE cxx_std_17)

