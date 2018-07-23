program test_polynomial
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_cfuns
  integer:: info
  call test_root([0.2d0, 0.3d0])
  call test_root([0.5d0, 0.5d0])
  call test_root([0.2d0, 0.3d0, 0.8d0])
  call test_root([0.2d0, 0.2d0, 0.8d0])
  call test_poly([-1d0, 0d0, 0d0, 1d0])
  call test_poly([1d0, 0d0, 1d0])
contains
  Subroutine test_root(R)
    real(c_double), intent(in) :: R(:)
    real(c_double) :: P(size(R)+1), Q(size(R))
    integer :: nr, i
    print *,"Initial roots:", R
    call mlf_polyFromRoots(R,P)
    print *, "P coefficients: ", P
    forall(i=1:size(R)) Q(i) = mlf_polyVal(P, R(i))
    print *, "P(root) shall be 0: ", Q
    nr = mlf_rootsPoly(P,Q)
    print *, "Estimated roots:", Q(1:nr)
  End Subroutine test_root
  Subroutine test_poly(P)
    real(c_double), intent(in) :: P(:)
    real(c_double) :: Q(size(P)-1)
    integer :: nr
    print *, "Polynomial: ",P
    nr = mlf_rootsPoly(P,Q)
    print *, "Estimated roots:", Q(1:nr)
  End Subroutine test_poly
end

