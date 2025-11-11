cmake_minimum_required(VERSION 3.16)
project(flux_burgers)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai REQUIRED)

add_executable(flux_burgers main.cpp)
target_link_libraries(flux_burgers PRIVATE samurai::samurai)
