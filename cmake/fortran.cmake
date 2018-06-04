
set (CMAKE_Fortran_FLAGS "-Wall -std=f2008ts -fimplicit-none")
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -funroll-loops")
set (CMAKE_Fortran_FLAGS_DEBUG "-ggdb -fsanitize=address")
#set (CMAKE_Fortran_FLAGS_DEBUG "-ggdb")

