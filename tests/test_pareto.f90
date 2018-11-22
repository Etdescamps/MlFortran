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

Program test_pareto
  use ieee_arithmetic
  use iso_c_binding
  use mlf_utils
  use mlf_matrix
  use mlf_rand
  use mlf_pareto
  use mlf_intf
  use test_common
  implicit none
  integer :: ND, NX, info

  info = mlf_init()

  If(COMMAND_ARGUMENT_COUNT()<2) Then
    PRINT *, "Error missing arguments: ND NX"
    STOP
  Endif
  ND = GetIntParameter(1)
  NX = GetIntParameter(2)
  PRINT *,"ND: ", ND, " NX: ", NX
  CALL testPareto()

  info = mlf_init()
Contains

  Subroutine testPareto()
    real(c_double), allocatable :: M(:,:), M2(:,:)
    integer(c_int), allocatable :: np(:), sel(:), idx(:)
    integer :: N, Mu,i
    ALLOCATE(idx(NX), M(ND, NX), M2(ND,NX))
    call RandN(M)
    np = mlf_getPareto(M, idx)
    N = size(np)
    PRINT *, np
    Do i=1,NX
      M2(:,i) = M(:,idx(i))
    End Do
    PRINT *, test_non_dominance(M2(:,:np(1)))
    PRINT *, test_dominance(M2(:,:np(1)), M2(:,np(1)+1:np(2)))
    Mu = NX/2
    ALLOCATE(sel(Mu))
    CALL mlf_paretoSelect(M, Mu, idx, np)
  End Subroutine testPareto

  logical Function is_dominated(P, V) result(r)
    real(c_double), intent(in) :: P(:,:), V(:)
    integer :: i, N
    N = size(P,2)
    Do i=1,N
      If(mlf_pareto_dominate(P(:,i), V) == 1) Then
        r = .TRUE.
        RETURN
      Endif
    End Do
    r = .FALSE.
  End Function is_dominated

  logical Function test_non_dominance(V) result(r)
    real(c_double), intent(in) :: V(:,:)
    integer :: i, j, N, k
    N = size(V,2)
    Do i=1,N-1
      Do j=i+1,N
        k = mlf_pareto_dominate(V(:,i), V(:,j))
        If(k == 1 .OR. k == -1) Then
          r = .FALSE.
          RETURN
        Endif
      End Do
    End Do
    r = .TRUE.
  End Function test_non_dominance

  logical Function test_dominance(P, V) result(r)
    real(c_double), intent(in) :: P(:,:), V(:,:)
    integer :: i, N
    N = size(V,2)
    Do i=1,N
      If(.NOT. is_dominated(P, V(:,1))) Then
        r = .FALSE.
        RETURN
      Endif
    End Do
    r = .TRUE.
  End Function test_dominance
End Program test_pareto

