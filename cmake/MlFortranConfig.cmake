get_filename_component(ML_FORTRAN_CM_DIR "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
set(ML_FORTRAN_INCLUDE_DIR ${ML_FORTRAN_CM_DIR}/../../include/MlFortran)
set(ML_FORTRAN_LIB_DIR ${ML_FORTRAN_CM_DIR}/../../lib)
include_directories(${ML_FORTRAN_INCLUDE_DIR})
link_directories(${ML_FORTRAN_LIB_DIR})

if(EXISTS "${ML_FORTRAN_LIB_DIR}/libmlfortran.so")
  message("-- MlFortran library found")
elseif(EXISTS "${ML_FORTRAN_LIB_DIR}/libmlfortran.dylib")
  message("-- MlFortran library found")
else()
  message(SEND_ERROR "Error: MlFortran library not found")
endif()

if(EXISTS "${ML_FORTRAN_INCLUDE_DIR}/mlf_cintf.h")
  message("-- MlFortran headers found")
else()
  message(SEND_ERROR "Error: MlFortran headers not found")
endif()

set(ML_FORTRAN_LIB mlfortran)

