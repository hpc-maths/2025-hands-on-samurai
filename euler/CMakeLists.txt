cmake_minimum_required(VERSION 4.1)
project(euler_tutorial)

find_package(samurai REQUIRED)

set(CMAKE_CXX_STANDARD 20)

add_executable(euler main.cpp)
target_link_libraries(euler PRIVATE samurai::samurai)

add_executable(sod_shock_2d main_sod_2d.cpp)
target_link_libraries(sod_shock_2d PRIVATE samurai::samurai)

add_executable(double_mach_reflection main_double_mach_reflection.cpp)
target_link_libraries(double_mach_reflection PRIVATE samurai::samurai)

add_executable(sedov_blast main_sedov_blast.cpp)
target_link_libraries(sedov_blast PRIVATE samurai::samurai)