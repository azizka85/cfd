cmake_minimum_required(VERSION 3.16)

project(1d_couette_flow_cranck_nicolson_time_variation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1d_couette_flow_cranck_nicolson_time_variation main.cpp solver.cpp)

target_link_libraries(1d_couette_flow_cranck_nicolson_time_variation PRIVATE
    slae
)
