# Slab-Pro: Depth-Averaged Material Point Method for Geophysical Mass Flows
# Copyright (C) 2026 ETH Zurich, SLF Davos
#
# Slab-Pro is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Slab-Pro is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Slab-Pro. If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.14)
project(SlabPro VERSION 1.0 LANGUAGES CXX)

# ----------------------
# C++ settings
# ----------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------
# Set default build type if none is specified
# ----------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug or Release" FORCE)
endif()

# ----------------------
# Compiler flags
# ----------------------
# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -fno-inline -fsanitize=address,undefined -fno-omit-frame-pointer")
# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# ----------------------
# Project-wide warnings (INTERFACE target)
# ----------------------
add_library(project_warnings INTERFACE)
target_compile_options(project_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion
)

# ----------------------
# Find core dependencies
# ----------------------
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GDAL REQUIRED)
find_package(nlohmann_json 3.10.5 REQUIRED)

# ----------------------
# Simulation library
# ----------------------
add_library(simulation_lib STATIC
    src/simulation.cpp
    src/tools.cpp
    src/save_function.cpp
    src/complex_terrain.cpp
    src/polygon_tools.cpp
)
target_include_directories(simulation_lib PUBLIC src)
target_link_libraries(simulation_lib
    PUBLIC Eigen3::Eigen
           GDAL::GDAL
           nlohmann_json::nlohmann_json
    PRIVATE project_warnings
)

# ----------------------
# CLI simulation executable
# ----------------------
add_executable(slab-pro src/main.cpp)
target_link_libraries(slab-pro
    PRIVATE simulation_lib
            project_warnings
)

# ----------------------
# Print current build type
# ----------------------
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
