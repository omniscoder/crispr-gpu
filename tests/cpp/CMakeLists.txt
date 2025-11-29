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
  add_test(NAME crispr_gpu_pytests
    COMMAND ${Python3_EXECUTABLE} -m pytest -q
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/python)
endif()
