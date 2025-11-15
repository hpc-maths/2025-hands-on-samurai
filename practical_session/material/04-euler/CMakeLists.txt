cmake_minimum_required(VERSION 4.1)
project(euler_tutorial)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai REQUIRED)

add_executable(euler_2d main_2d.cpp)
target_link_libraries(euler_2d PRIVATE samurai::samurai)


