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

Module mlf_sumtree
  Use ieee_arithmetic
  Use iso_c_binding
  IMPLICIT NONE
  PRIVATE

  Integer, Parameter :: hp_nAlloc0 = 4096

  Type, Public :: mlf_sumtree_double
    integer :: nIndiv
    real(c_double), allocatable :: V(:,:,:)
  Contains
    procedure :: init => sumtree_double_init
    procedure :: getTop => sumtree_double_getTop
    procedure :: getIndiv => sumtree_double_getIndiv
    procedure :: getAction => sumtree_double_getAction
    procedure :: realloc => sumtree_double_realloc
    procedure :: updateIndiv => sumtree_double_updateIndiv
    procedure :: setNIndiv => sumtree_double_setNIndiv
  End Type mlf_sumtree_double
Contains
  Subroutine sumtree_double_setNIndiv(this, nIndiv)
    class(mlf_sumtree_double), intent(inout) :: this
    integer, intent(in) :: nIndiv
    this%nIndiv = nIndiv
  End Subroutine sumtree_double_setNIndiv

  Integer Function sumtree_double_init(this, dimVect, nActions, nAlloc) Result(info)
    class(mlf_sumtree_double), intent(inout) :: this
    integer, intent(in) :: dimVect, nActions
    integer, intent(in), optional :: nAlloc
    integer :: nAlloc0
    info = 0
    nAlloc0 = hp_nAlloc0
    If(PRESENT(nAlloc)) Then
      Do While(nAlloc0 < nAlloc)
        nAlloc0 = 2*nAlloc0
      End Do
    Endif
    ! Deallocate the array if it does not fit the parameters
    If(ALLOCATED(this%V)) DEALLOCATE(this%V)
    ! Allocate the array if it is not allocated or deallocated
    If(.NOT.ALLOCATED(this%V)) Then
      ALLOCATE(this%V(dimVect, nActions, 0:nAlloc0-1), STAT = info)
    Endif
    this%V = 0
    this%nIndiv = 0
  End Function sumtree_double_init

  Integer Function sumtree_double_realloc(this, nRealloc) Result(info)
    class(mlf_sumtree_double), intent(inout) :: this
    integer, intent(in) :: nRealloc
    real(c_double), allocatable :: V(:,:,:)
    integer :: sV(3)
    sV = SHAPE(this%V)
    If(this%nIndiv > nRealloc) Then
      info = -1
      RETURN
    Endif
    ALLOCATE(V(sV(1), sV(2), 0:nRealloc-1), stat = info)
    If(info /= 0) RETURN
    V(:,:,0:this%nIndiv-1) = this%V(:,:,0:this%nIndiv-1)
    V(:,:,this%nIndiv:) = 0
    CALL MOVE_ALLOC(V, this%V)
  End Function sumtree_double_realloc

  Integer Function sumtree_double_getAction(this, vSel, val, sumV) Result(idAction)
    class(mlf_sumtree_double), intent(inout) :: this
    real(c_double), intent(in) :: vSel(:)
    real(c_double), intent(out) :: val
    real(c_double), intent(out), optional :: sumV
    real(c_double) :: X(SIZE(this%V, 2)), r
    X = MATMUL(vSel, this%V(:,:,0))
    CALL RANDOM_NUMBER(r)
    If(PRESENT(sumV)) Then
      sumV = SUM(x)
      r = r*sumV
    Else
      r = r*SUM(X)
    Endif
    Do idAction = 1, SIZE(this%V, 2)
      r = r - X(idAction)
      If(r <= 0) EXIT
    End Do
    idAction = MIN(SIZE(this%V, 2), idAction)
    val = X(idAction)
  End Function sumtree_double_getAction

  Real(c_double) Function sumtree_double_getTop(this, vSel, idAction, vActions) Result(val)
    class(mlf_sumtree_double), intent(inout) :: this
    real(c_double), intent(in) :: vSel(:)
    integer, intent(in), optional :: idAction
    real(c_double), intent(out), optional :: vActions(:)
    If(this%nIndiv == 0) Then
      val = 0
      If(PRESENT(vActions)) vActions = 0
      RETURN
    Endif
    If(PRESENT(vActions)) Then
      vActions = MATMUL(vSel, this%V(:,:,0))
      val = SUM(vActions)
    Else If(PRESENT(idAction)) Then
      val = DOT_PRODUCT(this%V(:,idAction,0), vSel)
    Else
      val = SUM(MATMUL(vSel, this%V(:,:,0)))
    Endif
  End Function sumtree_double_getTop

  Integer Function sumtree_double_updateIndiv(this, id, DV) Result(info)
    class(mlf_sumtree_double), intent(inout) :: this
    integer, intent(in) :: id
    real(c_double), intent(in) :: DV(:,:)
    integer :: j, k, nRealloc
    If(id >= this%nIndiv) Then
      nRealloc = SIZE(this%V, 3)
      If(id >= nRealloc) Then
        Do While(nRealloc <= id)
          nRealloc = nRealloc*2
        End Do
        info = this%realloc(nRealloc)
        If(info /= 0) RETURN
      Endif
      this%nIndiv = id
    Endif
    info = 0
    j = 1
    k = id-1
    this%V(:,:,0) = this%V(:,:,0) + DV
    Do While(k>0)
      If(BTEST(k,0)) Then
        this%V(:,:,j) = this%V(:,:,j) + DV
        j = ISHFT(j,1) + 1
      Else
        j = ISHFT(j,1)
      Endif
      k = ISHFT(k,-1)
    End Do
  End Function sumtree_double_updateIndiv

  Integer Function sumtree_double_getIndiv(this, vSel, idAction, val) Result(k)
    class(mlf_sumtree_double), intent(inout) :: this
    real(c_double), intent(in) :: vSel(:)
    real(c_double), intent(inout) :: val
    integer, intent(in) :: idAction
    integer :: j, u
    real(c_double) :: x, z, r
    j = 1; k = 0; u = 1
    x = val
    Do While(k+u < this%nIndiv)
      CALL RANDOM_NUMBER(r)
      z = DOT_PRODUCT(this%V(:,idAction,j), vSel)
      If(r*x < z) Then
        x = z
        j = ISHFT(j,1) + 1
        k = k + u
      Else
        x = x - z
        j = ISHFT(j,1)
      Endif
      u = ISHFT(u,1)
    End Do
    val = x
    k = k + 1
  End Function sumtree_double_getIndiv
End Module mlf_sumtree

