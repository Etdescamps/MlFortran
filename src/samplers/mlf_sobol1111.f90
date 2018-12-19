! Copyright (c) 2018 Etienne Descamps
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

Module mlf_sobol1111
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_sobol1111_table
  IMPLICIT NONE
  PRIVATE
  
  ! Generator of quasirandom Sobol sequences based on:
  ! Algorithm 659: Implementing Sobol's Quasirandom Sequence Generator
  ! from Paul Bratley and Bennett Fox
  ! ACM Transactions on Mathematical Software, Volume 14, Number 1, March 1988
  !
  ! Using improvements from Stephen Joe and Frances Kuo described in
  ! Remark on Algorithm 659: Impelmenting Sobol's Quasirandom
  ! Sequence Generator
  ! ACM Transations on Mathematical Software, Vol 29, No. 1, March 2003
  ! 
  ! More details on this method can be found at the address:
  ! http://web.science.unsw.edu.au/~fkuo/sobol/
  !
  ! The init procedure of the module has been inspired by the method i8_sobol,
  ! from the file sobol.f90 that can be found here:
  ! https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html

  Integer, Parameter :: s_dimMax = 1111
  Integer, Parameter :: s_logMax = 62
  Integer(8), Parameter :: s_mask = SHIFTR(-1_8,2)

  Type, Public :: mlf_sobol1111_sampler
    integer :: numDim
    integer(8) :: seed
    integer(8), allocatable :: table(:,:), last(:)
  Contains
    procedure :: init => mlf_sobol1111_init
    procedure :: mlf_sobol1111_sampleReal, mlf_sobol1111_sampleInt
    procedure :: mlf_sobol1111_sampleMatrix
    generic :: sample => mlf_sobol1111_sampleReal, mlf_sobol1111_sampleInt, &
      mlf_sobol1111_sampleMatrix
  End Type mlf_sobol1111_sampler
Contains

  Subroutine SetVect1(V, S, B)
    integer(8), intent(inout) :: V(:)
    integer, intent(in) :: B
    integer(1), dimension(B:), intent(in) :: S
    integer :: N
    N = SIZE(V)
    If(B > N) RETURN
    V(B:) = INT(S(B:N), KIND = 8)
  End Subroutine SetVect1

  Subroutine SetVect2(V, S, B)
    integer(8), intent(inout) :: V(:)
    integer, intent(in) :: B
    integer(2), dimension(B:), intent(in) :: S
    integer :: N
    N = SIZE(V)
    If(B > N) RETURN
    V(B:) = INT(S(B:N), KIND = 8)
  End Subroutine SetVect2
 
  Integer Function mlf_sobol1111_init(this, numDim) Result(info)
    class(mlf_sobol1111_sampler), intent(inout) :: this
    integer, intent(in) :: numDim
    integer(8) :: newV
    integer :: i, j, m, k
    integer(2) :: poly
    info = -1
    If(numDim > s_dimMax .OR. numDim <= 0) RETURN
    ALLOCATE(this%table(numDim, s_logMax))
    ASSOCIATE(V => this%table)
      V = 0
      v(1,:) = 1
      V(2:,1) = 1
      CALL SetVect1(V(:,2), Sobol1111_V2, LBOUND(Sobol1111_V2, 1))
      CALL SetVect1(V(:,3), Sobol1111_V3, LBOUND(Sobol1111_V3, 1))
      CALL SetVect1(V(:,4), Sobol1111_V4, LBOUND(Sobol1111_V4, 1))
      CALL SetVect1(V(:,5), Sobol1111_V5, LBOUND(Sobol1111_V5, 1))
      CALL SetVect1(V(:,6), Sobol1111_V6, LBOUND(Sobol1111_V6, 1))
      CALL SetVect1(V(:,7), Sobol1111_V7, LBOUND(Sobol1111_V7, 1))
      CALL SetVect2(V(:,8), Sobol1111_V8, LBOUND(Sobol1111_V8, 1))
      CALL SetVect2(V(:,9), Sobol1111_V9, LBOUND(Sobol1111_V9, 1))
      CALL SetVect2(V(:,10), Sobol1111_V10, LBOUND(Sobol1111_V10, 1))
      CALL SetVect2(V(:,11), Sobol1111_V11, LBOUND(Sobol1111_V11, 1))
      CALL SetVect2(V(:,12), Sobol1111_V12, LBOUND(Sobol1111_V12, 1))
      CALL SetVect2(V(:,13), Sobol1111_V13, LBOUND(Sobol1111_V13, 1))
      Do i = 1, numDim
        poly = Sobol1111_Poly(i)
        m = BIT_SIZE(poly)-LEADZ(poly)-1
        Do j = m+1, s_logMax
          newV = V(i, j-m)
          Do k = 1, m
            If(BTEST(poly, m-k)) newV = IEOR(newV, SHIFTL(V(i, j-k), k))
          End Do
          V(i,j) = newV
        End Do
      End Do
      FORALL(i=1:s_logMax-1) V(:,i) = SHIFTL(V(:,i), s_logMax-i)
    END ASSOCIATE
    info = 0
  End Function mlf_sobol1111_init

  Elemental Real(c_double) Function ConvToReal(X) Result(Y)
    integer(8), intent(in) :: X
    integer(8), parameter :: Z = TRANSFER(1d0, 1_8)
    Y = TRANSFER(IOR(Z, SHIFTR(X, 10)), 1d0)-1d0
  End Function ConvToReal

  Subroutine Sobol1111NextSample(this)
    class(mlf_sobol1111_sampler), intent(inout) :: this
    integer :: l, numDim
    numDim = SIZE(this%table, 1)
    If(.NOT. ALLOCATED(this%last)) Then
      ALLOCATE(this%last(numDim))
      this%last = 0
      this%seed = 0
    Else
      l = TRAILZ(NOT(this%seed))+1
      this%last = IEOR(this%last, this%table(:,l))
      this%seed = IAND(this%seed+1, s_mask)
    Endif
  End Subroutine Sobol1111NextSample

  Integer Function mlf_sobol1111_sampleInt(this, X) Result(info)
    class(mlf_sobol1111_sampler), intent(inout) :: this
    integer(c_int64_t), intent(out) :: X(:)
    CALL Sobol1111NextSample(this)
    X = SHIFTL(this%last, 1)
    info = 0
  End Function mlf_sobol1111_sampleInt

  Integer Function mlf_sobol1111_sampleReal(this, X) Result(info)
    class(mlf_sobol1111_sampler), intent(inout) :: this
    real(c_double), intent(out) :: X(:)
    CALL Sobol1111NextSample(this)
    X = ConvToReal(this%last)
    info = 0
  End Function mlf_sobol1111_sampleReal

  Integer Function mlf_sobol1111_sampleMatrix(this, X) Result(info)
    class(mlf_sobol1111_sampler), intent(inout) :: this
    real(c_double), intent(out) :: X(:,:)
    integer :: i
    Do i = 1, SIZE(X,2)
      info = mlf_sobol1111_sampleReal(this, X(:,i))
    End Do
    info = 0
  End Function mlf_sobol1111_sampleMatrix
End Module mlf_sobol1111

