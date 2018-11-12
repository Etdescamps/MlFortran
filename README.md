# MlFortran

**MlFortran** is a library designed for the writing and the development of Complex System models.
It is written in a modern object-oriented Fortran (2008) language, and has a C and C++ interface.
Its aim is to ease the development of complex models in lower level language such as C/C++ and Fortran by providing structures and algorithms that can be found in higher level programming language such as Matlab or Python/SciPy.
This library is designed to be sufficiently generalist to handle multiple kind of algorithms:
  - An ODE solver (based on Dormand-Prince RK-5) with advanced constraint and stop handling features.
  - Dimension reduction (linearisation of a function space)
  - Statistical inference (EM algorithm on gaussian mixture model)
  - Non-linear optimisation (currently only (C)MA-ES is implemented)
  - Diverse methods for random number generation and matrix operations
This library use an interface using HDF5 for saving checkpoints and intermediate values.
This interface permits also to get data from the C language.

**This project is in an early stage of development, don't expect this to be stable.**

## License
MlFortran is released under the BSD 3-Clause license.

## Requirements

This project has only been tested in POSIX environments (Ubuntu 18.04 and macOS).
It requires these dependencies:
  - CMake >= 3.5
  - GCC/G++ >= 7.0
  - GFortran >= 7.0
  - LAPACK and BLAS
  - HDF5: a version that contains a Fortran 90/2003 interface compiled with the same
    version of GNU Fortran. (Remove Anaconda from your $PATH if cmake choose the wrong
    version)

## Getting started
After having clone this repository, you can compile the project using CMake:
```console
foo@bar:MlFortran$ mkdir build 
foo@bar:MlFortran$ cd build
foo@bar:build$ cmake ..
foo@bar:build$ make
foo@bar:build$ sudo make install
```

For the instance, only a few demo are present for testing purposes.
You can try out CMA-ES to see if this implementation converges quickly on Rosenbrock function:
```console
foo@bar:build$ tests/test_cmaes
```




