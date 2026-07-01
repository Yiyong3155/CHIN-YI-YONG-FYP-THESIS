cmake_minimum_required(VERSION 3.19)
project(fusx LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python FIRST (required by nanobind)
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find required packages
find_package(MPI REQUIRED)
find_package(DOLFINX REQUIRED)

# Find nanobind
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# Create Python extension module
nanobind_add_module(
  _fusx
  STABLE_ABI
  NB_STATIC
  src/fusx.cpp
  src/forms.c
)

# Link libraries
target_link_libraries(_fusx PRIVATE dolfinx MPI::MPI_CXX)

# Include directories
target_include_directories(_fusx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Install
install(TARGETS _fusx LIBRARY DESTINATION fusx)
