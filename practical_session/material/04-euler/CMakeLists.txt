cmake_minimum_required(VERSION 4.1)
project(euler_tutorial)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai REQUIRED)

add_executable(euler main.cpp)
target_link_libraries(euler PRIVATE samurai::samurai)

