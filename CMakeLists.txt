cmake_minimum_required(VERSION 3.5)

project(MLFORTRAN)

set(VERSION 0.0.1)

enable_language(C Fortran CXX)

include(cmake/fortran.cmake)
include(cmake/c.cmake)
include(cmake/cpp.cmake)

find_package(LAPACK REQUIRED)
find_package(HDF5 COMPONENTS C Fortran REQUIRED)
#find_package(Boost 1.54 REQUIRED)

include_directories(include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src)
include_directories(${HDF5_INCLUDE_DIRS})
set(MLF_INCLUDES ${CMAKE_BINARY_DIR}/MlFortran)
include_directories(${MLF_INCLUDES})


set(MLF_LIB mlfortran)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(utils)
add_subdirectory(examples)

install(DIRECTORY include/ DESTINATION include/MlFortran FILES_MATCHING PATTERN "*.h")

