
add_executable(gauss_newton_example
    main.cpp
)

target_link_libraries(gauss_newton_example 
    gauss_newton
)

target_include_directories(gauss_newton_example
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)
