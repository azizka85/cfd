cmake_minimum_required(VERSION 3.16)

project(2d_heat_conduction_thermally_insulated_steady_state_sor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(2d_heat_conduction_thermally_insulated_steady_state_sor main.cpp solver.cpp)
