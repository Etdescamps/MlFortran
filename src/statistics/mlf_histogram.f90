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
    procedure :: init     => mlf_histogram_int_init
    procedure :: get      => mlf_histogram_int_get
    procedure :: addPoint => mlf_histogram_int_addPoint
  End Type mlf_histogram_int
Contains
  Integer Function mlf_histogram_int_init(this, X) Result(info)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    integer :: N
    N = SIZE(X)
    info = 0
    this%X = X
    ALLOCATE(this%V(N+1), this%W(N+1))
    this%V = 0; this%W = 0
  End Function mlf_histogram_int_init

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
      Where(this%V /= 0) X = this%W/this%V
    Endif
    info = 0
  End Function mlf_histogram_int_get

  Subroutine mlf_histogram_int_addPoint(this, x, N)
    class(mlf_histogram_int), intent(inout) :: this
    real(c_double), intent(in) :: x
    integer(c_int64_t), intent(in) :: N
    integer :: i
    If(x < this%X(1)) Then
      i = 1
    Else
      i = mlf_di_search(this%X, x) + 1
    Endif
    ASSOCIATE(V => this%V, W => this%W)
      V(i) = V(i) + N
      W(i) = W(i) + N*x
    END ASSOCIATE
  End Subroutine mlf_histogram_int_addPoint
End Module mlf_histogram

