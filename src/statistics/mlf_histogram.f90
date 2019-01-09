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

Module mlf_histogram
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Type, Public :: mlf_histogram_int
    real(c_double), allocatable :: X(:), W(:)
    integer(c_int64_t), allocatable :: V(:)
  Contains
    procedure :: mlf_histogram_int_init_withX, mlf_histogram_int_init_abN
    procedure :: get      => mlf_histogram_int_get
    procedure :: addPoint => mlf_histogram_int_addPoint
    procedure :: addVectSorted  => mlf_histogram_int_addVectSorted
    generic :: init => mlf_histogram_int_init_withX, mlf_histogram_int_init_abN
  End Type mlf_histogram_int
Contains

  Integer Function mlf_histogram_int_init_abN(this, a, b, N, with_W) Result(info)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: a, b
    logical, intent(in), optional :: with_W
    integer, intent(in) :: N
    integer :: i
    info = 0
    ALLOCATE(this%X(N), this%V(N+1))
    FORALL(i=1:N) this%X(i) = a+(b-a)*REAL(i-1,8)/REAL(N-1,8)
    this%V = 0
    If(PresentAndTrue(with_W)) Then
      ALLOCATE(this%W(N+1))
      this%W = 0
    Endif
  End Function mlf_histogram_int_init_abN

  Integer Function mlf_histogram_int_init_withX(this, X, with_W) Result(info)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    logical, intent(in), optional :: with_W
    integer :: N
    N = SIZE(X)
    info = 0
    this%X = X
    ALLOCATE(this%V(N+1))
    this%V = 0
    If(PresentAndTrue(with_W)) Then
      ALLOCATE(this%W(N+1))
      this%W = 0
    Endif
  End Function mlf_histogram_int_init_withX

  Integer Function mlf_histogram_int_get(this, W, X) Result(info)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(out) :: W(:)
    real(c_double), intent(out), optional :: X(:)
    integer :: N
    real(c_double) :: S
    info = -1
    N = SIZE(W)
    If(SIZE(this%V) /= N) RETURN
    S = SUM(this%V)
    W = this%V/S
    If(PRESENT(X)) Then
      If(SIZE(X) /= N) RETURN
      X(1) = this%X(1)
      X(2:N-1) = 0.5d0*(this%X(1:N-2)+this%X(2:N-1))
      X(N) = this%X(N-1)
      If(ALLOCATED(this%W)) Where(this%V /= 0) X = this%W/this%V
    Endif
    info = 0
  End Function mlf_histogram_int_get

  Subroutine mlf_histogram_int_addVectSorted(this, X)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    integer :: i, j
    real(c_double) :: V
    j = 1; V = this%X(j)
    Do i = 1,SIZE(X)
      Do While(X(i) >= V)
        If(j >= SIZE(this%X)) Then
          V = HUGE(1d0)
          EXIT
        Else
          j = j+1
          V = this%X(j)
        Endif
      End Do
      this%V(j) = this%V(j) + 1
      If(ALLOCATED(this%W)) this%W(j) = this%W(j) + X(i)
    End Do
  End Subroutine mlf_histogram_int_addVectSorted

  Subroutine mlf_histogram_int_addPoint(this, x, N)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: x
    integer(c_int64_t), intent(in), optional :: N
    integer :: i
    If(x < this%X(1)) Then
      i = 1
    Else
      i = mlf_di_search(this%X, x) + 1
    Endif
    ASSOCIATE(V => this%V, W => this%W)
      If(PRESENT(N)) Then
        V(i) = V(i) + N
        If(ALLOCATED(this%W)) W(i) = W(i) + N*x
      Else
        V(i) = V(i) + 1
        If(ALLOCATED(this%W)) W(i) = W(i) + x
      Endif
    END ASSOCIATE
  End Subroutine mlf_histogram_int_addPoint
End Module mlf_histogram

