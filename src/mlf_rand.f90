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

Module mlf_rand
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Public :: RandSign, RandN, RandD, mlf_rand_class, randPerm, randInt, mlf_randNArrayC
  Public :: Rand3DSurf, RandFromIArr, mlf_randN2
  
  Type, Public, Abstract :: mlf_1DRealSampler
  Contains
    procedure(mlf_1DRealSampler_fun), deferred :: random
    procedure :: sample => mlf_1DRealSampler_sample
  End Type mlf_1DRealSampler

  Abstract Interface
    Function mlf_1DRealSampler_fun(this)
      use iso_c_binding
      import :: mlf_1DRealSampler
      class(mlf_1DRealSampler), intent(inout) :: this
      real(c_double) :: mlf_1DRealSampler_fun
    End Function mlf_1DRealSampler_fun
  End Interface

  Interface RandSign  
    module procedure mlf_RandSignV
    module procedure mlf_RandSignM
  End Interface RandSign
  Interface RandN
    module procedure mlf_randNVect
    module procedure mlf_randNMatrix
  End Interface RandN
  Interface RandD
    module procedure mlf_randD0
    module procedure mlf_randD1
  End Interface RandD
  Interface Rand3DSurf
    module procedure Rand3DSurf_float
    module procedure Rand3DSurf_double
  End Interface Rand3DSurf

Contains
  Subroutine mlf_1DRealSampler_sample(this, X)
    class(mlf_1DRealSampler), intent(inout) :: this
    real(c_double), intent(out) :: X(:)
    integer :: i
    Do i = LBOUND(X,1), UBOUND(X,1)
      X(i) = this%random()
    End Do
  End Subroutine mlf_1DRealSampler_sample
  
  Subroutine mlf_RandSignV(V)
    real(c_double), intent(out) :: V(:)
    call RANDOM_NUMBER(V)
    where (V<0.5)
      V=-1
    elsewhere
      V=1
    end where
  End Subroutine mlf_RandSignV

  Subroutine mlf_RandSignM(M)
    real(c_double), intent(out) :: M(:,:)
    call RANDOM_NUMBER(M)
    where (M<0.5)
      M=-1
    elsewhere
      M=1
    end where
  End Subroutine mlf_RandSignM

  ! Utility function for generating random double number in [0,1[
  real(c_double) Function mlf_randD0()
    call RANDOM_NUMBER(mlf_randD0)
  End Function mlf_randD0
  
  ! Utility function for generating random double number in [0,A[
  real(c_double) Function mlf_randD1(A) result(R)
    real(c_double), intent(in) :: A
    call RANDOM_NUMBER(R)
    R = R * A
  End Function mlf_randD1
  
  ! Utility function for generating random integer (from 0 to N-1)
  ! (if N is huge, it is better to use a more precise generator)
  integer(c_int) Function randInt(N)
    integer(c_int), intent(in) :: N
    real(c_double) :: r
    call RANDOM_NUMBER(r)
    randInt = MIN(floor(r*N, kind=c_int), N-1)
  End Function randInt

  ! Random point on surface using Marsaglia's method
  Subroutine Rand3DSurf_double(X, r, X0)
    real(c_double), intent(in) :: r
    real(c_double), intent(in), optional :: X0(3)
    real(c_double), intent(out) :: X(3)
    real(c_double) :: U(2), W
 10 CALL RANDOM_NUMBER(U)
    U = 2.0*U-1.0
    W = DOT_PRODUCT(U,U)
    If(W > 1.0d0) GOTO 10
    if(PRESENT(X0)) then
      X(1:2) = X0(1:2) + 2d0*r*U*sqrt(1d0-W)
      X(3) = X0(3) + r*(1d0-2d0*W)
    else
      X(1:2) = 2d0*r*U*sqrt(1d0-W)
      X(3) = r*(1d0-2d0*W)
    endif
  End Subroutine Rand3DSurf_double

  ! Random point on surface using Marsaglia's method
  Subroutine Rand3DSurf_float(X, r, X0)
    real(c_float), intent(in) :: r
    real(c_float), intent(in), optional :: X0(3)
    real(c_float), intent(out) :: X(3)
    real(c_float) :: U(2), W
 10 CALL RANDOM_NUMBER(U)
    U = 2.0*U-1.0
    W = DOT_PRODUCT(U,U)
    If(W > 1.0d0) GOTO 10
    if(PRESENT(X0)) then
      X(1:2) = X0(1:2) + 2.0*r*U*sqrt(1.0-W)
      X(3) = X0(3) + r*(1.0-2.0*W)
    else
      X(1:2) = 2.0*r*U*sqrt(1.0-W)
      X(3) = r*(1.0-2.0*W)
    endif
  End Subroutine Rand3DSurf_float

  integer(c_int) Function RandFromIArr(X) Result(N)
    integer(c_int) :: X(:)
    integer :: i
    If(SIZE(X) == 1) Then
      N = X(LBOUND(X,1))
      RETURN
    Endif
    i = randInt(SIZE(X))+LBOUND(X,1)
    N = X(i)
  End Function RandFromIArr


  ! Generate random class from a vector proportional to the cumulative probabilies
  Subroutine mlf_rand_class(W, id) 
    real(c_double), intent(in) :: W(:)
    integer(c_int), intent(out) :: id(:)
    integer :: i
    real(c_double) :: r
    Do i=1,size(id)
      call RANDOM_NUMBER(r)
      id(i) = mlf_di_search(W, r*W(size(W)))
    End Do
  End Subroutine mlf_rand_class

  ! Generate a random permutation of an integer array
  Subroutine randPerm(idx)
    integer(c_int), intent(inout) :: idx(:)
    integer :: i, j, k, N
    N = size(idx)
    Do i=1,N-1
      k = idx(i)
      j = randInt(N+1-i)+i
      idx(i) = idx(j)
      idx(j) = k
    End Do
  End Subroutine randPerm

  ! Generate two random normal variables using Box-Muller(Marsaglia) method
  Subroutine mlf_randN2(R)
    real(c_double), intent(out) :: R(2)
    real(c_double) :: rsq
    ! GFORTRAN uses the xorshift 1024 method as RNG.
    ! This method is thread-safe with recent versions of GFORTRAN.
    ! This shall be sufficient for this kind of problems.
 10 CALL RANDOM_NUMBER(R)
    R = 2.0*R-1.0
    rsq = DOT_PRODUCT(R,R)
    If(rsq > 1.0d0 .OR. rsq == 0d0) GOTO 10 ! Outside the circle (probability=pi/4)
    R = R*sqrt(-2.0*log(rsq)/rsq)
  end subroutine mlf_randN2

  ! Works only with contiguous array
  Subroutine mlf_randNArrayC(V, N) bind(C, name="mlf_randN")
    real(c_double), intent(out) :: V(*)
    real(c_double) :: r(2)
    integer(c_int), value :: N
    integer :: i
    Do i=1,N,2
      call mlf_randN2(V(i:i+1))
    End do
    if(btest(N,0)) then ! If N is odd
      call mlf_randN2(r)
      V(N) = r(1)
    endif
  End Subroutine mlf_randNArrayC

  ! Use previous method for generating an arbitrary size array of random variables
  Subroutine mlf_randNVect(V)
    real(c_double), intent(out) :: V(:)
    real(c_double) :: r(2)
    integer :: i
    Do i=lbound(V,1),ubound(V,1)-1,2
      call mlf_randN2(r)
      V(i:i+1) = r(:)
    End Do
    if(i<ubound(V,1)) then
      call mlf_randN2(r)
      V(i+1) = r(1)
    endif
  End Subroutine mlf_randNVect

  Subroutine mlf_randNMatrix(M, sigma, X0, C12)
    real(c_double), intent(out) :: M(:,:)
    real(c_double), optional, intent(in) :: sigma, X0(:), C12(:,:)
    integer :: i
    Do i=lbound(M,2),ubound(M,2)
      call mlf_randNVect(M(:,i))
    End Do
    if (present(sigma)) then
      M = sigma*M
    endif
    if(PRESENT(C12)) M = matmul(C12,M)
    if(PRESENT(X0)) then
      Do i=lbound(M,2),ubound(M,2)
        M(:,i) = X0 + M(:,i)
      End Do
    endif
  End Subroutine mlf_randNMatrix

End Module mlf_rand

