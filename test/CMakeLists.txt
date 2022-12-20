enable_testing()
include(GoogleTest)

add_executable(gauss_newton_test
  test.cpp
)

target_link_libraries(gauss_newton_test
  gauss_newton
  gtest
  gtest_main
)

target_include_directories(gauss_newton_test
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)

gtest_discover_tests(gauss_newton_test)
