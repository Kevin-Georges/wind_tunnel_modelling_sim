cmake_minimum_required(VERSION 3.16)
project(wind_tunnel_modelling_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Phase 1: command-line physics harness only. No windowing/rendering/UI
# dependencies -- those come back once the LBM solver is validated and
# phase 2 (a windowed app) actually starts.
add_executable(wind_tunnel_physics
    main.cpp
    src/utils/Logger.cpp
    src/geometry/Mesh.cpp
    src/geometry/Transform.cpp
    src/geometry/BoundingBox.cpp
    src/geometry/Voxelizer.cpp
    src/io/StlLoader.cpp
    src/aero/WindTunnel.cpp
    src/aero/LatticeBoltzmannSolver.cpp
)

target_include_directories(wind_tunnel_physics PRIVATE src)
