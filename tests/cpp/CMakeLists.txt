add_executable(crispr_gpu_tests
  test_encoding.cpp
  test_index.cpp
  test_gpu_parity.cpp
)

target_link_libraries(crispr_gpu_tests PRIVATE crispr_gpu Catch2::Catch2WithMain)
target_compile_definitions(crispr_gpu_tests PRIVATE TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../data\")
add_test(NAME crispr_gpu_tests COMMAND crispr_gpu_tests)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('pytest') else 1)"
    RESULT_VARIABLE PYTEST_STATUS
  )
  if(PYTEST_STATUS EQUAL 0)
    add_test(NAME crispr_gpu_pytests
      COMMAND ${Python3_EXECUTABLE} -m pytest -q
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/python)
  else()
    message(STATUS "pytest not found; skipping Python tests")
  endif()
endif()
