cmake_minimum_required(VERSION 4.1)
project(euler_tutorial)

find_package(samurai REQUIRED)

set(CMAKE_CXX_STANDARD 20)
add_executable(euler main.cpp)
target_link_libraries(euler PRIVATE samurai::samurai)