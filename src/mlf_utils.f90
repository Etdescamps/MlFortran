! Copyright (c) 2017-2018 Etienne Descamps
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation and/or
!    other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Module mlf_utils
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_cfuns
  IMPLICIT NONE
  PRIVATE

  Public :: Xchange, QSortIdx, QSort, Mean, PrintMatrix, mlf_countcl, mlf_reversecl, SetIf
  Public :: mlf_reverseid, mlf_cumulativefromvect, mlf_di_search, mlf_printmatrix_c, mlf_isSorted
  Public :: mlf_solve4DPoly, mlf_rootsPoly, mlf_polyFromRoots, mlf_polyVal

  Interface SetIf
    module procedure mlf_setIfInt
    module procedure mlf_setIfInt64
    module procedure mlf_setIfDouble
  End Interface SetIf
  Interface Xchange
    module procedure mlf_xchangeInt
    module procedure mlf_xchangeInt64
    module procedure mlf_xchangeDouble
    module procedure mlf_xchangeFloat
  End Interface Xchange
  Interface QSortIdx
    module procedure mlf_sortM
    module procedure mlf_sortV
  End Interface QSortIdx
  Interface QSort
    module procedure mlf_sortOutM
    module procedure mlf_sortOutV
  End Interface QSort
  Interface Mean
    module procedure mlf_MeanV
    module procedure mlf_MeanM
    module procedure mlf_MeanMV
  End Interface Mean
  Interface PrintMatrix
    module procedure mlf_printmatrix
    module procedure mlf_printmatrix2
  End Interface PrintMatrix
  real(c_double), public, parameter :: mlf_PI = acos(-1d0)
Contains

  Subroutine mlf_setIfInt(dest, def, val)
    integer(c_int), intent(out) :: dest
    integer(c_int), intent(in) :: def
    integer(c_int), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfInt

  Subroutine mlf_setIfInt64(dest, def, val)
    integer(c_int64_t), intent(out) :: dest
    integer(c_int64_t), intent(in) :: def
    integer(c_int64_t), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfInt64

  Subroutine mlf_setIfDouble(dest, def, val)
    real(c_double), intent(out) :: dest
    real(c_double), intent(in) :: def
    real(c_double), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfDouble


  ! Utility functions for exchanging to variables
  Pure Elemental Subroutine mlf_xchangeInt(a,b)
    integer(c_int32_t), intent(inout) :: a,b
    integer(c_int32_t) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeInt

  ! Utility functions for exchanging to variables
  Pure Elemental Subroutine mlf_xchangeInt64(a,b)
    integer(c_int64_t), intent(inout) :: a,b
    integer(c_int64_t) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeInt64


  Pure Elemental Subroutine mlf_xchangeDouble(a,b)
    real(c_double), intent(inout) :: a,b
    real(c_double) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeDouble

  Pure Elemental Subroutine mlf_xchangeFloat(a,b)
    real(c_float), intent(inout) :: a,b
    real(c_float) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeFloat

  Integer Function rootsQuadratic(R, b, c) result(nr)
    real(c_double), intent(in) :: b, c
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    nr = 0
    delta = 0.25d0*b**2-c
    if(abs(delta) < 1d-12*abs(c)) then
      R(1) = -0.5d0*b
      nr = 1
      RETURN
    endif
    if(delta<0d0) RETURN
    delta = sqrt(delta)
    R(1) = -0.5d0*b-delta
    R(2) = -0.5d0*b+delta
    nr = 2
  End Function rootsQuadratic

  ! Find the real solutions of the equation x^3+p*x+q=0
  Integer Function rootsCardano(R, p, q) result(nr)
    real(c_double), intent(in) :: p, q
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    complex(c_double) :: a
    complex(c_double), parameter :: j = (-0.5d0, 8.660254037844386467869d-1)
    nr = 1 ! At least one solution
    delta = -4d0*p**3-27d0*q**2
    if(abs(delta)<1d-12*abs(p+q)) then ! discriminant = 0
      if(abs(p) < 1d-12) then ! p == q == 0, so only the solution z=0
        R(1) = 0
      else
        R(1) = 3d0*q/p
        R(2) = (-3d0*q)/(2d0*p) ! The double root of the equation
        ! So in this case: x^3+p*x+q=(x-R(1)*(x-R(2))^2
        nr = 2
      endif
      RETURN
    endif
    if(delta<0) then ! One real solution (+ 2 complex solutions)
      delta = sqrt(-delta/27d0)
      R(1) = c_cbrt(0.5d0*(-q+delta))+c_cbrt(0.5d0*(-q-delta))
    else ! Three real solutions
      nr = 3
      a = cmplx(-0.5d0*q, 0.5*sqrt(delta/27d0), kind=c_double)**(1d0/3d0)
      ! Get one complex root of (-q+i*sqrt(delta/27))/2
      ! the other complex roots are a*j and a*j^2
      ! The solutions have the form a*j^k+conjg(a*j^k) = 2*real(a*j^k)
      R(1) = 2d0*REAL(a, kind= c_double)
      a = a*j
      R(2) = 2d0*REAL(a, kind= c_double)
      a = a*j
      R(3) = 2d0*REAL(a, kind= c_double)
    endif
  End Function rootsCardano

  ! Roots of a depressed quartic using Ferrari's method
  Integer Function rootsFerrari(A, p, q, r) result(nr)
    real(c_double), intent(in) :: p, q, r
    real(c_double), intent(out) :: A(:)
    real(c_double) :: u, v, mv, R2(3), m
    integer :: i, j
    nr = 0
    mv = max(abs(p), abs(q), abs(r))
    if(abs(q) < 1d-12*mv) then
      ! Solve x**4+p*x**2+r=0
      j = rootsQuadratic(R2, p, r)
      Do i=1,j
        if(R2(i) < -1d-12*mv) CYCLE
        if(R2(i) < 1d-12*mv) then
          nr = nr+1
          A(nr) = 0d0
        else
          nr = nr+2
          u = sqrt(R2(i))
          A(nr-1:nr) = [-u, u]
        endif
      End Do
      RETURN
    endif
    ASSOCIATE(b => p, c => (0.25d0*p**2-r), d => (1d0/8d0*q**2))
      ! Equation P(x)=x^3+b*x^2+c*x+d is transformed in depressed form:
      ! x is replaced by z = x+b/3, so the equation has the form z^3+p*z+q
      j = rootsCardano(R2, c-(b**2)/3d0, b/27d0*(2d0*b**2-9d0*c)+d)
      R2(1:j) = R2(1:j) - b/3d0 ! x = z-b/3
    END ASSOCIATE
    m = Maxval(R2(1:j))
    if(m<-1d-12*mv) RETURN ! No solutions
    if(m<1d-12*mv) then
      nr = rootsQuadratic(A, 0d0, 0.5*p)
      RETURN
    endif
    u = sqrt(2d0*m)
    v = 0.5d0*q/u
    nr = rootsQuadratic(A, u, 0.5d0*p+m-v)
    nr = nr + rootsQuadratic(A(nr+1:), -u, 0.5d0*p+m+v)
  End Function rootsFerrari

  ! Find real roots of a polynomial
  Integer Function rootsPoly(A, R) result(nr)
    real(c_double), intent(in) :: A(:)
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    integer :: N
    N = size(A)
    nr = 0
    select case (N)
      case(0)
        RETURN
      case(1)
        R(1) = -A(1)
        nr = 1
      case(2) ! Solve analytically the second order equation
        delta = 0.25d0*A(2)**2-A(1)
        if(abs(delta) < 1d-12*abs(A(1))) then
          R(1) = -0.5d0*A(2)
          nr = 1
          RETURN
        endif
        if(delta<0d0) RETURN
        delta = sqrt(delta)
        R(1) = -0.5d0*A(2)-delta
        R(2) = -0.5d0*A(2)+delta
        nr = 2
      case(3) ! Use Cardano's method for 3rd order equations
        ASSOCIATE(b => A(3), c => A(2), d => A(1))
          ! Equation P(x)=x^3+b*x^2+c*x+d is transformed in depressed form:
          ! x is replaced by z = x+b/3, so the equation has the form z^3+p*z+q
          nr = rootsCardano(R, c-(b**2)/3d0, b/27d0*(2d0*b**2-9d0*c)+d)
          R(1:nr) = R(1:nr) - b/3d0 ! x = z-b/3
        END ASSOCIATE
      case(4) ! Use Ferrari's method
        ASSOCIATE(b => A(4), c => A(3), d => A(2), e => A(1))
          nr = rootsFerrari(R, c-3d0/8d0*b**2, 1d0/8d0*b**3-0.5*b*c+d, &
            e-3d0/256d0*b**4-0.25*b*d+1d0/16d0*b**2*c)
          R(1:nr) = R(1:nr) - b/4d0
        END ASSOCIATE
      case default
    end select
  End Function rootsPoly

  ! Interface for rootsPoly, check if the top coefficent of A can be negleted
  Integer Function mlf_rootsPoly(A, R) result(nr)
    real(c_double), intent(in) :: A(:)
    real(c_double), intent(out) :: R(:)
    real(c_double) :: m
    integer :: i, N
    m = maxval(abs(A))
    N = size(A)
    nr = 0
    Do i=N,1,-1
      if(abs(A(i)) > m*1e-12) EXIT
    End Do
    if(i>1) nr = rootsPoly(A(1:i-1)/A(i), R)
  End Function mlf_rootsPoly

  real(c_double) Function mlf_solve4DPoly(t, a0, a1, a2, a3, a4) result(x)
    real(c_double), intent(in) :: t, a0, a1, a2, a3, a4
    real(c_double) :: R(4)
    integer :: nr, i=-1,j
    nr = mlf_rootsPoly([a0, a1, a2, a3, a4], R)
    Do i=1,nr
      if(R(i)>t) EXIT
    End Do
    Do j=i+1,nr
      if(R(j)>t .AND. R(j)<R(i)) i = j
    End Do
    if(i>0) then
      x = R(i)
    else
      x = IEEE_VALUE(x, IEEE_QUIET_NAN)
    endif
  End Function mlf_solve4DPoly

  Pure real(c_double) Function mlf_polyVal(P, x) result(Y)
    real(c_double), intent(in) :: P(:), x
    integer :: N, i
    N =size(P)
    Y = P(N)
    Do i=N-1,1,-1
      Y = x*Y+P(i)
    End Do
  End Function mlf_polyVal

  Subroutine mlf_polyFromRoots(R, P)
    real(c_double), intent(in) :: R(:)
    real(c_double), intent(out) :: P(:)
    real(c_double) :: a
    integer :: k, N
    N = size(R)
    if(N == 1) then
      P(1) = -R(1)
      P(2) = 1d0
      RETURN
    endif
    P(1) = R(1)*R(2)
    P(2) = -R(1)-R(2)
    Do k=3,N
      a = R(k)
      P(k) = P(k-1)-a
      P(2:k-1) = P(1:k-2)-a*P(2:k-1)
      P(1) = -a*P(1)
    End Do
    P(N+1) = 1d0
  End Subroutine mlf_polyFromRoots

  ! Test function written to test QSort function
  logical Function mlf_isSorted(Y) result(T)
    real(c_double) :: Y(:)
    integer :: i
    T=.FALSE.
    do i=1,size(Y)-1
      if(Y(i) > Y(i+1)) RETURN 
    end do
    T=.TRUE.
  End Function mlf_isSorted
  
  ! Quick sort algorithm that uses FORTRAN90 array structures
  ! This method sort the index array idx (Y(idx) will generate a sorted Y array)
  ! Mu is the minimal number of sorted element in the begin of the array
  Subroutine mlf_sortV(Y, idx, mu)
    real(c_double), intent(in), contiguous :: Y(:)
    integer, intent(in), optional :: mu
    integer(c_int), intent(out) :: idx(:)
    integer(c_int) :: N, mu0
    N = size(Y,1)
    mu0 = N
    if(present(mu)) mu0 = mu
    call c_qsort(Y, idx, N, 1, 1, mu0)
  End Subroutine mlf_sortV
  Subroutine mlf_sortM(Y, idx, mu)
    real(c_double), intent(in), contiguous :: Y(:,:)
    integer, intent(in), optional :: mu
    integer(c_int), intent(out) :: idx(:)
    integer(c_int) :: ND, N, mu0
    N = size(Y,2); ND = size(Y,1)
    mu0 = N
    if(present(mu)) mu0 = mu
    call c_qsort(Y, idx, N, ND, ND, mu0)
  End Subroutine mlf_sortM
  
  Function mlf_sortOutV(Y, mu) result(Z)
    real(c_double), intent(in), contiguous :: Y(:)
    integer, intent(in), optional :: mu
    real(c_double), allocatable :: Z(:)
    integer, allocatable :: idx(:)
    integer(c_int) :: N, mu0
    N = size(Y,1)
    mu0 = N
    if(present(mu)) mu0 = mu
    allocate(idx(N))
    call c_qsort(Y, idx, N, 1, 1, mu0)
    Z = Y(idx)
  End Function mlf_sortOutV 
  Function mlf_sortOutM(Y, mu) result(Z)
    real(c_double), intent(in), contiguous :: Y(:,:)
    integer, intent(in), optional :: mu
    real(c_double), allocatable :: Z(:,:)
    integer, allocatable :: idx(:)
    integer(c_int) :: ND, N, mu0
    N = size(Y,2); ND = size(Y,1)
    mu0 = N
    if(present(mu)) mu0 = mu
    allocate(idx(N))
    call c_qsort(Y, idx, N, ND, ND, mu0)
    Z = Y(:,idx)
  End Function mlf_sortOutM

  ! Useful functions for getting a mean value
  real(c_double) function mlf_MeanV(V)
    real(c_double), intent(in) :: V(:)
    mlf_MeanV = sum(V)/real(size(V),kind=8)
  end function mlf_MeanV

  real(c_double) function mlf_MeanM(M)
    real(c_double), intent(in) :: M(:,:)
    mlf_MeanM = sum(M)/real(size(M), kind=8)
  end function mlf_MeanM

  function mlf_MeanMV(M, dim) result(V)
    real(c_double), allocatable :: V(:)
    real(c_double), intent(in) :: M(:,:)
    integer :: dim
    V = sum(M, dim=dim)/real(size(M,dim), kind=8)
  end function mlf_MeanMV

  ! Count class apearance within an integer array
  integer(c_int) function mlf_countcl(id, cnt) result(N)
    integer(c_int), intent(in) :: id(:)
    integer(c_int), intent(out) :: cnt(:)
    integer :: l, u, i, k
    N = 0
    l = lbound(cnt,1)
    u = ubound(cnt,1)
    cnt = 0
    do i=1,size(id,1)
      k = id(i)
      if(k<l .OR. k>u) CYCLE ! Ignore elements outside the array
      cnt(k) = cnt(k) + 1
      N = N+1 ! Return the number of detected elements
    end do
  end function mlf_countcl

  ! Reverse class association array
  ! Use last appearance of the class
  ! Keep the value of cl if cl not present within id
  integer(c_int) function Mlf_reversecl(id, cl) result(N)
    integer(c_int), intent(in) :: id(:)
    integer(c_int), intent(inout) :: cl(:)
    integer :: l, u, i, k
    N = 0
    l = lbound(cl,1)
    u = ubound(cl,1)
    do i=1,size(id,1)
      k = id(i)
      if(k<l .OR. k>u) CYCLE ! Ignore elements outside the array
      cl(k) = i
      N = N+1 ! Return the number of detected elements
    end do
  End Function mlf_reversecl

  Function mlf_reverseid(id) result(cl)
    integer(c_int), intent(in) :: id(:)
    integer(c_int), allocatable :: cl(:)
    integer :: info
    ALLOCATE(cl(size(id)))
    info = mlf_reversecl(id,cl)
  End Function mlf_reverseid
 
  ! Utility function generating cumulative vector
  Function mlf_cumulativefromvect(V) result(W)
    real(c_double), intent(in) :: V(:)
    real(c_double) :: W(size(V)), S
    integer :: i
    S=0
    Do i =1,size(V,1)
      S = S + V(i)
      W(i) = S
    End Do
  End Function mlf_cumulativefromvect

  ! Do a dichotomous search for a value 
  integer(c_int) Function mlf_di_search(V, x) result(k)
    real(c_double), intent(in) :: V(:), x
    integer :: i, j
    i = lbound(V, 1); j = ubound(V,1)
    Do While(j-i>1)
      k = ishft(i+j, -1)
      if(V(k)<x) then
        i = k+1
      else
        j = k
      endif
    End Do
    if(x <= V(i)) then 
      k = i
    else
      k = j
    endif
  End Function mlf_di_search
  ! Print matrix (useful for debug)
  Subroutine mlf_printmatrix(M)
    real(c_double), intent(in) :: M(:,:)
    integer :: i
    Do i=1,size(M,1)
      write(*, "(24f6.2)") M(i,:)
    End Do
  End Subroutine mlf_printmatrix
  ! Print matrix (useful for debug)
  Subroutine mlf_printmatrix2(M, P)
    real(c_double), intent(in) :: M(:,:), P(:,:)
    integer :: i, j
    Do i=1,size(M,1)
      write(*, "(24f6.2)") [([M(i,j), P(i,j)], j=1, size(M,2))]
    End Do
  End Subroutine mlf_printmatrix2

  Subroutine mlf_printmatrix_c(M, nL, nC) bind(C, name="mlf_printMat")
    integer(c_int), value :: nL, nC
    real(c_double), intent(in) :: M(nL, nC)
    call mlf_printmatrix(M)
  End Subroutine mlf_printmatrix_c
End Module mlf_utils

