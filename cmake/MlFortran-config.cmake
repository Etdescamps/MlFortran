get_filename_component(ML_FORTRAN_CM_DIR "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
set(ML_FORTRAN_INCLUDE_DIR ${ML_FORTRAN_CM_DIR}/../../include/MlFortran)
set(ML_FORTRAN_LIB_DIR ${ML_FORTRAN_CM_DIR}/../../lib)
include_directories(${ML_FORTRAN_INCLUDE_DIR})
link_directories(${ML_FORTRAN_LIB_DIR})



