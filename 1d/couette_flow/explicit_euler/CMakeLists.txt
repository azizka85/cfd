cmake_minimum_required(VERSION 3.16)

project(1d_couette_flow_explicit_euler LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1d_couette_flow_explicit_euler main.cpp solver.cpp)
