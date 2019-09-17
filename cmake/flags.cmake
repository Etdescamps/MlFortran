
#set(MLF_DEBUG_FLAGS "-ggdb -fsanitize=address")
set(MLF_DEBUG_FLAGS "-ggdb")
#set(MLF_RELEASE_FLAGS "-O3 -march=native -ftree-vectorize -ffast-math -funroll-loops")
#set(MLF_RELEASE_FLAGS "-O3 -march=native -ftree-vectorize -funroll-loops")
set(MLF_RELEASE_FLAGS "-O3 -march=native -ftree-vectorize -funroll-loops -fopenmp")
#set(MLF_RELEASE_FLAGS "-O3 -march=native -ftree-vectorize -funroll-loops -fopenmp -ggdb")
#set(MLF_RELEASE_FLAGS "-O3 -march=native -ftree-vectorize -funroll-loops -fopenmp -pg")
#set(MLF_RELEASE_FLAGS "-Og -ggdb")


