cmake_minimum_required(VERSION 3.16)
project(SimplexSolver)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(SimplexSolver
        main.cpp
        solver.hpp)
