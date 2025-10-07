cmake_minimum_required(VERSION 3.16)

project(1d_linear_advection_tvd_implicit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1d_linear_advection_tvd_implicit main.cpp solver.cpp)

target_link_libraries(1d_linear_advection_tvd_implicit PRIVATE
    slae
)
