cmake_minimum_required (VERSION 3.14)
project(gauss_newton)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(src) 
add_subdirectory(example) 

enable_testing ()
add_subdirectory(test)
