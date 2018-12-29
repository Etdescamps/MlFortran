
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -std=f2008ts -fimplicit-none")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${MLF_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${MLF_DEBUG_FLAGS}")

