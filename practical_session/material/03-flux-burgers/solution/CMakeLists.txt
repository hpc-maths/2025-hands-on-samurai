cmake_minimum_required(VERSION 3.16)
project(burgers)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai REQUIRED)

add_executable(burgers_1d main_1d.cpp)
target_link_libraries(burgers_1d PRIVATE samurai::samurai)

add_executable(burgers_2d main_2d.cpp)
target_link_libraries(burgers_2d PRIVATE samurai::samurai)

add_executable(burgers_viscous_2d main_viscous_2d.cpp)
target_link_libraries(burgers_viscous_2d PRIVATE samurai::samurai)
