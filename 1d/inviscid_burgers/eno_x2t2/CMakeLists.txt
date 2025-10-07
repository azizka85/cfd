cmake_minimum_required(VERSION 3.16)

project(1d_inviscid_burgers_eno_x2t2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1d_inviscid_burgers_eno_x2t2 main.cpp solver.cpp)
