Program sobol_test
  Use mlf_sobol1111
  type(mlf_sobol1111_sampler) :: s
  integer, parameter :: l = 64
  real(8) :: X(l, 5)
  integer :: info, i
  info = s%init(5)
  Do i=1,l
    info = s%sample(X(i,:))
    PRINT *, REAL(X(i,:), KIND = 4)
  End Do
End

