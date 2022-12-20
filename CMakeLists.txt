cmake_minimum_required (VERSION 3.14)
project(gauss_newton)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src) 

enable_testing ()
add_subdirectory(test)
