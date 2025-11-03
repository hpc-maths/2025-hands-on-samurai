cmake_minimum_required(VERSION 3.16)
project(naive_burgers)

find_package(samurai REQUIRED)

add_executable(naive_burgers main.cpp)
target_link_libraries(naive_burgers PRIVATE samurai::samurai)
