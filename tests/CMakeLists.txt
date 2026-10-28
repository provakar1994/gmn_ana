#----------------------------------------------------------------------------
# Unit tests for the gmnana library.
# Each test is a small executable (tests/test_*.cpp) that links the library and
# returns 0 if all of its checks pass. Run them all from the build directory with:
#   ctest --output-on-failure
# Small input files used by the tests live in tests/data/.

SET(tests
  test_lookup_table
  test_simc_rmc
)

foreach(t ${tests})
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${t} PRIVATE ${PACKAGE})
  target_compile_definitions(${t} PRIVATE GMNANA_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
  add_test(NAME ${t} COMMAND ${t})
endforeach()
