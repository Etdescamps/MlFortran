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
  Public :: InitOrDefault
  Public :: mlf_reverseid, mlf_cumulativefromvect, mlf_di_search, mlf_printmatrix_c, mlf_isSorted

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
  Interface InitOrDefault
    module procedure mlf_initOrDefault_bool
    module procedure mlf_initOrDefault_int
    module procedure mlf_initOrDefault_real
    module procedure mlf_initOrDefault_double
    module procedure mlf_initOrDefault_int64
  End Interface InitOrDefault
  real(c_double), public, parameter :: mlf_PI = acos(-1d0)
Contains

  Subroutine mlf_initOrDefault_bool(var, def, opt)
    logical, intent(out) :: var
    logical, intent(in) :: def
    logical, intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_bool

  Subroutine mlf_initOrDefault_int(var, def, opt)
    integer, intent(out) :: var
    integer, intent(in) :: def
    integer, intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_int

  Subroutine mlf_initOrDefault_int64(var, def, opt)
    integer(8), intent(out) :: var
    integer(8), intent(in) :: def
    integer(8), intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_int64

  Subroutine mlf_initOrDefault_real(var, def, opt)
    real, intent(out) :: var
    real, intent(in) :: def
    real, intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_real

  Subroutine mlf_initOrDefault_double(var, def, opt)
    real(8), intent(out) :: var
    real(8), intent(in) :: def
    real(8), intent(in), optional :: opt
    if(present(opt)) then
      var = opt
    else
      var = def
    endif
  End Subroutine mlf_initOrDefault_double

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

