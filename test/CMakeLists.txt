include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(gauss_newton_test
  test.cpp
)

target_link_libraries(gauss_newton_test
  gauss_newton
  GTest::gtest_main
)

target_include_directories(gauss_newton_test
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)

include(GoogleTest)
gtest_discover_tests(gauss_newton_test)
