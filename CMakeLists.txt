cmake_minimum_required(VERSION 3.10)  # CMake version check
project(mpi-molecular-simulation)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard
set(CMAKE_CXX_FLAGS "-O3")

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp utils.cpp)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(body3 ${SOURCE_FILES})

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(body3 ${MPI_C_LIBRARIES})