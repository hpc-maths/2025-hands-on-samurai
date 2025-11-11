cmake_minimum_required(VERSION 3.16)
project(samurai_setup_test)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai REQUIRED)

add_executable(samurai_setup_test main.cpp)
target_link_libraries(samurai_setup_test PRIVATE samurai::samurai)
