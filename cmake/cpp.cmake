
set (CMAKE_CXX_FLAGS "-Wall -std=c++17 -fPIC")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -funroll-loops")
set (CMAKE_CXX_FLAGS_DEBUG "-ggdb -fsanitize=address")

