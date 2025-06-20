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

add_subdirectory(1d/couette_flow/free_surface/cranck_nicolson_time_variation/1_order_bc)
add_subdirectory(1d/couette_flow/free_surface/cranck_nicolson_time_variation/2_order_bc)

add_subdirectory(1d/couette_flow/gravity/cranck_nicolson_time_variation)

add_subdirectory(1d/couette_flow/high_order_compact/cranck_nicolson_time_variation)

add_subdirectory(1d/heat_conduction/thermally_insulated/cranck_nicolson_time_variation)
add_subdirectory(1d/heat_conduction/thermally_insulated/steady_state)
