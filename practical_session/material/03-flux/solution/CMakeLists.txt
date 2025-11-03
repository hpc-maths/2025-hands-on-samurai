cmake_minimum_required(VERSION 3.16)
project(burgers)

find_package(samurai REQUIRED)

add_executable(burgers main.cpp)
target_link_libraries(burgers PRIVATE samurai::samurai)

add_executable(burgers_2d main_2d.cpp)
target_link_libraries(burgers_2d PRIVATE samurai::samurai)

add_executable(burgers_viscous_2d main_viscous_2d.cpp)
target_link_libraries(burgers_viscous_2d PRIVATE samurai::samurai)

add_executable(burgers_viscous_vector_2d main_viscous_vector_2d.cpp)
target_link_libraries(burgers_viscous_vector_2d PRIVATE samurai::samurai)
