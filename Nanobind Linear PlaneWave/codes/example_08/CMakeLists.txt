cmake_minimum_required(VERSION 3.20)
project(fusx)

# -------------------------------------------------------------------
# Compiler Flags
# -------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "-Ofast -march=native ${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_C_FLAGS   "-Ofast -march=native ${CMAKE_C_FLAGS} -Wall")

# -------------------------------------------------------------------
# Warning if CMake is run manually instead of via scikit-build-core
# -------------------------------------------------------------------
if (NOT SKBUILD)
  message(WARNING "
  This CMake file is intended to be executed via 'scikit-build-core'.

  For users:
    pip install .

  For developers:
    pip install nanobind scikit-build-core[pyproject]
    pip install --no-build-isolation -ve .
  (Optionally add -Ceditable.rebuild=true for auto-rebuild.)
  ")
endif()

# -------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------
if (NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

find_package(Python 3.10
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
)

# -------------------------------------------------------------------
# Build type
# -------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
  )
endif()

# -------------------------------------------------------------------
# Detect nanobind installation directory
# -------------------------------------------------------------------
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT
)

find_package(nanobind CONFIG REQUIRED)

# -------------------------------------------------------------------
# Build Python Extension Module
# -------------------------------------------------------------------
nanobind_add_module(
  _fusx
  # STABLE_ABI   # Uncomment for Python >= 3.12 stable ABI
  src/fusx.cpp
  src/forms.c
)

# Fix: Override nanobind's -Os flag with -O3 for performance
target_compile_options(_fusx PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:-O3>
)
# Remove -Os explicitly
string(REPLACE "-Os" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
string(REPLACE "-Os" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

target_link_libraries(_fusx PRIVATE dolfinx)

# -------------------------------------------------------------------
# Installation for scikit-build-core
# -------------------------------------------------------------------
install(TARGETS _fusx LIBRARY DESTINATION fusx)