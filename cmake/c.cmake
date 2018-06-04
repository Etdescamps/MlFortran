

set (CMAKE_C_FLAGS "-Wall -std=c99 -fPIC")
set (CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -funroll-loops")
set (CMAKE_C_FLAGS_DEBUG "-ggdb -fsanitize=address")
#set (CMAKE_C_FLAGS_DEBUG "-ggdb")


