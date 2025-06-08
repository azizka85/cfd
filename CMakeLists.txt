cmake_minimum_required(VERSION 3.16)

project(cfd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(slae)

add_subdirectory(1d/couette_flow/explicit_euler)
add_subdirectory(1d/couette_flow/steady_state)
add_subdirectory(1d/couette_flow/implicit_euler)
add_subdirectory(1d/couette_flow/cranck_nicolson)
add_subdirectory(1d/couette_flow/cranck_nicolson_time_variation)
