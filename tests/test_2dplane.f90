Program test_2dplane
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_2dplane_integration
  Use test_functions
  IMPLICIT NONE


  type(gauss_fun) :: f
  real(c_double) :: s, s1, kr
  integer :: i, j, k
  type(mlf_2d_h_val) :: D
  s = f%integrateOnPlane(1d-4, 12, 0.5d0, 0.5d0, 0d0, 0d0, 1d0, 1d0)
  PRINT*, "Tolerance: ", 1d-4
  PRINT *, s, (s-0.25d0*ACOS(-1d0))/(0.25d0*ACOS(-1d0)), f%nsteps
  s1 = 0
  k = INT(SQRT(REAL(f%nsteps)))
  kr = REAL(k)
  Do i = 0, k
    Do j = 0, k
      CALL f%getValDer(REAL(i)/kr, REAL(j)/kr, D)
      If((i == 0 .OR. i == k) .AND. (j == 0 .OR. j == k)) Then
        s1 = s1+ 0.25d0*D%val/kr**2
      Else If(i == 0 .OR. i == k .OR. j == 0 .OR. j == k) Then
        s1 = s1+ 0.5d0*D%val/kr**2
      Else
        s1 = s1+ D%val/kr**2
      Endif
    End Do
  End Do
  PRINT *, (k+1)**2, (s1-0.25d0*ACOS(-1d0))/(0.25d0*ACOS(-1d0))
End
