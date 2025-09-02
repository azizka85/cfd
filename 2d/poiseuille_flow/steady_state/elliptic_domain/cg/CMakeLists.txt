cmake_minimum_required(VERSION 3.16)

project(2d_poiseuille_flow_steady_state_elliptic_domain_cg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(2d_poiseuille_flow_steady_state_elliptic_domain_cg main.cpp solver.cpp)

target_link_libraries(2d_poiseuille_flow_steady_state_elliptic_domain_cg PRIVATE
    slae
)
