cmake_minimum_required(VERSION 3.16)

project(1d_inviscid_burgers_step_upwind_conservative LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1d_inviscid_burgers_step_upwind_conservative main.cpp solver.cpp)
