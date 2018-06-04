cmake_minimum_required(VERSION 3.5)

project(MLFORTRAN)

set(VERSION 0.0.1)

enable_language(C Fortran)

include(cmake/fortran.cmake)
include(cmake/c.cmake)

find_package(LAPACK REQUIRED)
find_package(HDF5 COMPONENTS C Fortran REQUIRED)

include_directories(include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src)
include_directories(${HDF5_INCLUDE_DIRS})

set(MLF_LIB mlfortran)

add_subdirectory(src)
add_subdirectory(tests)

#install(TARGETS ${MLF_LIB} LIBRARY DESTINATION lib)

