cmake_minimum_required(VERSION 3.5)

project(MLFORTRAN)

set(VERSION 0.0.1)

enable_language(C Fortran CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Select the type of build: Debug Release" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(cmake/flags.cmake)
include(cmake/fortran.cmake)
include(cmake/c.cmake)
include(cmake/cpp.cmake)

find_package(LAPACK REQUIRED)
find_package(HDF5 COMPONENTS C Fortran REQUIRED)
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else(OpenMP_FOUND)
  message("OpenMP not found")
endif(OpenMP_FOUND)
find_package(Boost 1.54 REQUIRED)

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
install(DIRECTORY include/ DESTINATION include/MlFortran FILES_MATCHING PATTERN "*.hpp")
install(FILES cmake/MlFortranConfig.cmake DESTINATION share/MlFortran)

