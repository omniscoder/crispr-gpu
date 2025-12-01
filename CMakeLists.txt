cmake_minimum_required(VERSION 3.18)
project(crispr-gpu VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(CRISPR_GPU_BUILD_PYTHON "Build Python bindings" ON)
option(CRISPR_GPU_BUILD_CLI "Build CLI binary" ON)
option(CRISPR_GPU_ENABLE_TESTS "Build tests" ON)
option(CRISPR_GPU_ENABLE_CUDA "Enable CUDA kernels" ON)

# Warnings
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Dependencies
include(FetchContent)

# pybind11 for Python bindings
if(CRISPR_GPU_BUILD_PYTHON)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.12.0
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# nlohmann/json for table loading
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

# Catch2 for tests
if(CRISPR_GPU_ENABLE_TESTS)
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.3
  )
  FetchContent_MakeAvailable(catch2)
endif()

set(CRISPR_GPU_BASE_SOURCES
  src/types.cpp
  src/genome_index.cpp
  src/genome_index_build.cpp
  src/fm_index.cpp
  src/engine.cpp
  src/scoring.cpp
  src/version.cpp
)

if(CRISPR_GPU_ENABLE_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  list(APPEND CRISPR_GPU_BASE_SOURCES src/kernels.cu)
else()
  add_definitions(-DCRISPR_GPU_CPU_ONLY)
endif()

add_library(crispr_gpu STATIC ${CRISPR_GPU_BASE_SOURCES})
target_include_directories(crispr_gpu PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(crispr_gpu PRIVATE nlohmann_json::nlohmann_json)
target_compile_definitions(crispr_gpu PUBLIC CRISPR_GPU_DATA_DIR="${CMAKE_INSTALL_PREFIX}/share/crispr_gpu")

if(CRISPR_GPU_ENABLE_CUDA)
  set_target_properties(crispr_gpu PROPERTIES CUDA_STANDARD 14 CUDA_STANDARD_REQUIRED ON)
  target_compile_definitions(crispr_gpu PUBLIC CRISPR_GPU_ENABLE_CUDA)
  target_include_directories(crispr_gpu PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  # Default to a reasonable arch; allow user override via CMAKE_CUDA_ARCHITECTURES
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set_target_properties(crispr_gpu PROPERTIES CUDA_ARCHITECTURES "60;70;75;80")
  endif()
endif()

# CLI
if(CRISPR_GPU_BUILD_CLI)
  add_executable(crispr-gpu src/cli_main.cpp)
  target_link_libraries(crispr-gpu PRIVATE crispr_gpu)
endif()

# Kernel microbenchmark (developer aid)
if(CRISPR_GPU_ENABLE_CUDA)
  add_executable(kernel_microbench benchmarks/kernel_microbench.cu)
  set_target_properties(kernel_microbench PROPERTIES CUDA_STANDARD 14 CUDA_STANDARD_REQUIRED ON)
  target_link_libraries(kernel_microbench PRIVATE CUDA::cudart)
endif()

# Python bindings
if(CRISPR_GPU_BUILD_PYTHON)
  add_subdirectory(python)
endif()

# Tests
if(CRISPR_GPU_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()

# Install rules (minimal for now)
install(TARGETS crispr_gpu
  EXPORT CrisprGpuTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)
install(TARGETS nlohmann_json EXPORT CrisprGpuTargets)
install(DIRECTORY include/ DESTINATION include)

install(EXPORT CrisprGpuTargets
  FILE CrisprGpuTargets.cmake
  NAMESPACE crispr_gpu::
  DESTINATION lib/cmake/crispr_gpu)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/CrisprGpuConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CrisprGpuConfig.cmake
  INSTALL_DESTINATION lib/cmake/crispr_gpu)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/CrisprGpuConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CrisprGpuConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CrisprGpuConfigVersion.cmake
  DESTINATION lib/cmake/crispr_gpu)
