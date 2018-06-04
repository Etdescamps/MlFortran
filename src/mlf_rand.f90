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
Contains
  Subroutine mlf_RandSignV(V)
    real(c_double), intent(out) :: V(:)
    call random_number(V)
    where (V<0.5)
      V=-1
    elsewhere
      V=1
    end where
  End Subroutine mlf_RandSignV

  Subroutine mlf_RandSignM(M)
    real(c_double), intent(out) :: M(:,:)
    call random_number(M)
    where (M<0.5)
      M=-1
    elsewhere
      M=1
    end where
  End Subroutine mlf_RandSignM

  ! Utility function for generating random double number in [0,1[
  real(c_double) function mlf_randD0()
    call random_number(mlf_randD0)
  end function mlf_randD0
  
  ! Utility function for generating random double number in [0,A[
  real(c_double) function mlf_randD1(A) result(R)
    real(c_double), intent(in) :: A
    call random_number(R)
    R = R * A
  end function mlf_randD1
  
  ! Utility function for generating random integer
  ! (if N is huge, it is better to use a more precise generator)
  integer(c_int) function randInt(N)
    integer(c_int), intent(in) :: N
    real(c_double) :: r
    call random_number(r)
    randInt = floor(r*N, kind=c_int)
  end function randInt

  ! Generate random class from a vector proportional to the cumulative probabilies
  subroutine mlf_rand_class(W, id) 
    real(c_double), intent(in) :: W(:)
    integer(c_int), intent(out) :: id(:)
    integer :: i
    real(c_double) :: r
    do i=1,size(id)
      call random_number(r)
      id(i) = mlf_di_search(W, r*W(size(W)))
    end do
  end subroutine mlf_rand_class

  ! Generate a random permutation of an integer array
  subroutine randPerm(idx)
    integer(c_int), intent(inout) :: idx(:)
    integer :: i, j, k, N
    N = size(idx)
    do i=1,N-1
      k = idx(i)
      j = randInt(N+1-i)+i
      idx(i) = idx(j)
      idx(j) = k
    end do
  end subroutine randPerm
  ! Generate two random normal variables using Box-Muller(Marsaglia) method
  subroutine mlf_randN2(R)
    real(c_double), intent(out) :: R(2)
    real(c_double) :: rsq
    do
      ! GFORTRAN uses the xorshift 1024 method as RNG.
      ! This method is thread-safe with recent versions of GFORTRAN.
      ! This shall be sufficient for this kind of problems.
      call random_number(R)
      R=2.0*R-1.0
      rsq = dot_product(R,R)
      if(rsq > 0.0 .AND. rsq <= 1.0) EXIT ! Outside the circle (probability=pi/4)
    end do
    R = R*sqrt(-2.0*log(rsq)/rsq)
  end subroutine mlf_randN2

  ! Works only with contiguous array
  Subroutine mlf_randNArrayC(V, N) bind(C, name="mlf_randN")
    real(c_double), intent(out) :: V(*)
    real(c_double) :: r(2)
    integer(c_int), value :: N
    integer :: i
    do i=1,N,2
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
    do i=lbound(V,1),ubound(V,1)-1,2
      call mlf_randN2(r)
      V(i:i+1) = r(:)
    end do
    if(i<ubound(V,1)) then
      call mlf_randN2(r)
      V(i+1) = r(1)
    endif
  End Subroutine mlf_randNVect

  Subroutine mlf_randNMatrix(M, sigma, X0, C12)
    real(c_double), intent(out) :: M(:,:)
    real(c_double), optional, intent(in) :: sigma, X0(:), C12(:,:)
    integer :: i
    do i=lbound(M,2),ubound(M,2)
      call mlf_randNVect(M(:,i))
    end do
    if (present(sigma)) then
      M = sigma*M
    endif
    if(present(C12)) M = matmul(C12,M)
    if(present(X0)) then
      do i=lbound(M,2),ubound(M,2)
        M(:,i) = X0 + M(:,i)
      end do
    endif
  End Subroutine mlf_randNMatrix

End Module mlf_rand

