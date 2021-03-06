! Copyright (c) 2017-2019 Etienne Descamps
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
  Public :: InitOrDefault, FirstEltVal, PresentAndTrue, OutsideBounds
  Public :: mlf_reverseid, mlf_cumulativefromvect, mlf_di_search, mlf_printmatrix_c, mlf_isSorted
  Public :: Median, AddCstr

  Interface FirstEltVal
    module procedure mlf_firstIntVal
    module procedure mlf_firstDoubleVal
  End Interface FirstEltVal

  Interface SetIf
    module procedure mlf_setIfInt
    module procedure mlf_setIfInt64
    module procedure mlf_setIfDouble
  End Interface SetIf
  Interface Xchange
    module procedure mlf_xchangeInt
    module procedure mlf_xchangeInt64
    module procedure mlf_xchangeDouble
    module procedure mlf_xchangeComplex4
    module procedure mlf_xchangeComplex8
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
  Interface InitOrDefault
    module procedure mlf_initOrDefault_bool
    module procedure mlf_initOrDefault_int
    module procedure mlf_initOrDefault_real
    module procedure mlf_initOrDefault_double
    module procedure mlf_initOrDefault_int64
  End Interface InitOrDefault
  Interface Median
    module procedure mlf_median_real8
  End Interface Median
  real(c_double), public, parameter :: mlf_PI = ACOS(-1d0)
Contains
  Real(c_double) Function mlf_median_real8(X, idX) Result(Median)
    real(c_double), intent(in) :: X(:)
    integer, intent(in) :: idX(:)
    integer :: N
    N = SIZE(X)
    If(BTEST(N,0)) Then
      Median = X(idX((N+1)/2))
    Else
      Median = 0.5d0*(X(idX(N/2))+X(idX(N/2+1)))
    Endif
  End Function mlf_median_real8

  Elemental Real(c_double) Function OutsideBounds(X, XMin, XMax) Result(cstr)
    real(c_double), intent(in) :: X, XMin, XMax
    If(X < XMin) Then
      cstr = XMin - X
    Else If(X > XMax) Then
      cstr = XMax - X
    Else
      cstr = 0
    Endif
  End Function OutsideBounds

  Integer Function AddCstr(V, N, i) Result(k)
    integer, intent(inout) :: V(:)
    integer, intent(in) :: N, i
    If(ANY(V(1:N) == i)) Then
      k = N
    Else
      V(N+1) = i
      k = N+1
    Endif
  End Function AddCstr

  Logical Function PresentAndTrue(x) Result(y)
    logical, optional :: x
    y = .FALSE.
    If(PRESENT(x)) y = x
  End Function PresentAndTrue

  Pure Subroutine mlf_initOrDefault_bool(var, def, opt)
    logical, intent(out) :: var
    logical, intent(in) :: def
    logical, intent(in), optional :: opt
    If(PRESENT(opt)) Then
      var = opt
    Else
      var = def
    Endif
  End Subroutine mlf_initOrDefault_bool

  Pure Subroutine mlf_initOrDefault_int(var, def, opt)
    integer, intent(out) :: var
    integer, intent(in) :: def
    integer, intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_int

  Pure Subroutine mlf_initOrDefault_int64(var, def, opt)
    integer(8), intent(out) :: var
    integer(8), intent(in) :: def
    integer(8), intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_int64

  Pure Subroutine mlf_initOrDefault_real(var, def, opt)
    real, intent(out) :: var
    real, intent(in) :: def
    real, intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_real

  Pure Subroutine mlf_initOrDefault_double(var, def, opt)
    real(8), intent(out) :: var
    real(8), intent(in) :: def
    real(8), intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_double

  Pure Subroutine mlf_setIfInt(dest, def, val)
    integer(c_int), intent(out) :: dest
    integer(c_int), intent(in) :: def
    integer(c_int), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfInt

  Pure Subroutine mlf_setIfInt64(dest, def, val)
    integer(c_int64_t), intent(out) :: dest
    integer(c_int64_t), intent(in) :: def
    integer(c_int64_t), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfInt64

  Pure Subroutine mlf_setIfDouble(dest, def, val)
    real(c_double), intent(out) :: dest
    real(c_double), intent(in) :: def
    real(c_double), intent(in), optional :: val
    if(present(val)) then
      dest = val
    else
      dest = def
    endif
  End Subroutine mlf_setIfDouble

  Pure Integer Function mlf_firstIntVal(array, val) Result(k)
    integer(c_int32_t), intent(in) :: array(:), val
    integer :: i, N
    N = size(array)
    k = -1 ! Not found
    Do i = 1,N
      if(array(i) == val) then
        k = i
        RETURN
      endif
    End Do
  End Function mlf_firstIntVal

  Pure Integer Function mlf_firstDoubleVal(array, val) Result(k)
    real(c_double), intent(in) :: array(:), val
    integer :: i, N
    N = size(array)
    k = -1 ! Not found
    Do i = 1,N
      if(array(i) == val) then
        k = i
        RETURN
      endif
    End Do
  End Function mlf_firstDoubleVal

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

  Pure Elemental Subroutine mlf_xchangeComplex4(a,b)
    complex(kind = 4), intent(inout) :: a,b
    complex(kind = 4) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeComplex4

  Pure Elemental Subroutine mlf_xchangeComplex8(a,b)
    complex(kind = 8), intent(inout) :: a,b
    complex(kind = 8) :: k
    k = a
    a = b
    b = k
  End Subroutine mlf_xchangeComplex8


  ! Test function written to test QSort function
  Pure Logical Function mlf_isSorted(Y) result(T)
    real(c_double), intent(in) :: Y(:)
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
    integer(c_int), allocatable :: idx(:)
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
    integer(c_int), allocatable :: idx(:)
    integer(c_int) :: ND, N, mu0
    N = size(Y,2); ND = size(Y,1)
    mu0 = N
    if(present(mu)) mu0 = mu
    allocate(idx(N))
    call c_qsort(Y, idx, N, ND, ND, mu0)
    Z = Y(:,idx)
  End Function mlf_sortOutM

  ! Useful functions for getting a mean value
  Pure Real(c_double) function mlf_MeanV(V)
    real(c_double), intent(in) :: V(:)
    mlf_MeanV = sum(V)/real(size(V),kind=8)
  end function mlf_MeanV

  Pure Real(c_double) Function mlf_MeanM(M)
    real(c_double), intent(in) :: M(:,:)
    mlf_MeanM = sum(M)/real(size(M), kind=8)
  End Function mlf_MeanM

  Pure Function mlf_MeanMV(M, dim) Result(V)
    real(c_double), allocatable :: V(:)
    real(c_double), intent(in) :: M(:,:)
    integer, intent(in) :: dim
    V = SUM(M, DIM=dim)/REAL(SIZE(M,dim), KIND=8)
  End Function mlf_MeanMV

  ! Count class apearance within an integer array
  Integer(c_int) Function mlf_countcl(id, cnt, l) result(N)
    integer(c_int), intent(in) :: id(:), l
    integer(c_int), intent(out) :: cnt(l:)
    integer :: u, i, k
    N = 0
    u = UBOUND(cnt,1)
    cnt = 0
    Do i=1,SIZE(id,1)
      k = id(i)
      If(k<l .OR. k>u) CYCLE ! Ignore elements outside the array
      cnt(k) = cnt(k) + 1
      N = N+1 ! Return the number of detected elements
    End Do
  End Function mlf_countcl

  ! Reverse class association array
  ! Use last appearance of the class
  ! Keep the value of cl if cl not present within id
  Integer(c_int) Function mlf_reversecl(id, cl, l) Result(N)
    integer(c_int), intent(in) :: id(:), l
    integer(c_int), intent(inout) :: cl(l:)
    integer :: u, i, k
    N = 0
    u = UBOUND(cl,1)
    Do i=1,SIZE(id,1)
      k = id(i)
      If(k<l .OR. k>u) CYCLE ! Ignore elements outside the array
      cl(k) = i
      N = N+1 ! Return the number of detected elements
    End Do
  End Function mlf_reversecl

  Function mlf_reverseid(id) Result(cl)
    integer(c_int), intent(in) :: id(:)
    integer(c_int), allocatable :: cl(:)
    integer :: info
    ALLOCATE(cl(SIZE(id)))
    info = mlf_reversecl(id,cl,1)
  End Function mlf_reverseid
 
  ! Utility function generating cumulative vector
  Pure Function mlf_cumulativefromvect(V) result(W)
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
  Pure Integer(c_int) Function mlf_di_search(V, x) result(k)
    real(c_double), intent(in) :: V(:), x
    integer :: i, j
    i = 1; j = SIZE(V)
    Do While(j-i>1)
      k = (i+j)/2
      If(V(k) < x) Then
        i = k
      Else
        j = k
      Endif
    End Do
    If(x < V(j)) Then 
      k = i
    Else
      k = j
    Endif
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

