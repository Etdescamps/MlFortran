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

Module mlf_gridmodel
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  IMPLICIT NONE
  PRIVATE

  Type, Public, extends(mlf_obj_model) :: mlf_2dgrid
    class(mlf_approx_model), pointer :: funmodel
    real(c_double), pointer :: grid(:,:,:)
    real(c_double), pointer :: XMin, XMax, YMin, YMax
  Contains
    procedure :: initF => mlf_GridModelInit
  End Type mlf_2dgrid

  Type, Public, extends(mlf_approx_model) :: mlf_2dgrid_model
    class(mlf_2dgrid), pointer :: top
  Contains
    procedure :: getValue => mlf_GridModelFunValue
    procedure :: getProj => mlf_GridModelGetProjection
  End Type mlf_2dgrid_model
Contains
  
  ! Model initilisator
  Subroutine mlf_model_funbasis_init(this, top)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_2dgrid), intent(in), target :: top
    ALLOCATE(mlf_2dgrid_model :: this)
    select type(this)
      class is (mlf_2dgrid_model)
        this%top => top
    end select
  End Subroutine mlf_model_funbasis_init

  ! C wrapper for init function
  type(c_ptr) Function c_2dgridmodel_init(cmodel, XMin, XMax, YMin, YMax, nW, nX0, nY0) &
      bind(C, name="mlf_2dgridModelInit")
    type(c_ptr), value :: cmodel
    class(mlf_model), pointer :: model
    type(mlf_2dgrid), pointer :: x
    class (*), pointer :: obj
    real(c_double), value :: XMin, XMax, YMin, YMax
    integer(c_int), value :: nW, nX0, nY0
    integer :: info
    c_2dgridmodel_init = C_NULL_PTR
    model => mlf_getModel(cmodel)
    if(.NOT. associated(model)) RETURN
    select type(model)
      class is (mlf_approx_model)
        ALLOCATE(x)
        info = x%initF(model, XMin, XMax, YMin, YMax, nW, nX0, nY0)
        if(info < 0) RETURN
        obj => x
        c_2dgridmodel_init = c_allocate(obj)
    end select
  End Function c_2dgridmodel_init

  Integer Function mlf_GridModelInit(this, fmodel, XMin, XMax, YMin, YMax, nW, nX0, nY0, data_handler) result(info)
    class(mlf_2dgrid), intent(inout) :: this
    class(mlf_approx_model), target :: fmodel
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: XMin, XMax, YMin, YMax
    integer(c_int), intent(in), optional :: nX0, nY0, nW
    real(c_double), allocatable :: X(:,:,:)
    type(mlf_rsc_numFields) :: numFields = mlf_rsc_numFields(0,4,1)
    integer :: i
    integer(c_int64_t) :: ndGrid(3)
    ndGrid = -1
    if(present(nW)) ndGrid(1) = nW
    if(present(nX0)) ndGrid(2) = nX0
    if(present(nY0)) ndGrid(3) = nY0
    info = mlf_arr_init(this, numFields, data_handler)
    if(info < 0) RETURN
    call this%addRPar(numFields, this%XMin, "XMin")
    call this%addRPar(numFields, this%XMax, "XMax")
    call this%addRPar(numFields, this%YMin, "YMin")
    call this%addRPar(numFields, this%YMax, "YMax")
    info = this%add_RMatrix3D(numFields, ndGrid, this%grid, C_CHAR_"grid", data_handler = data_handler)
    if(info < 0) RETURN
    this%funmodel => fmodel
    if(.NOT. present(data_handler)) then
      ALLOCATE(X(2,nX0,nY0))
      this%XMin = XMin; this%XMax = XMax; this%YMin = YMin; this%YMax = YMax;
      forall(i=1:nX0) X(1, i, :) = Xmin+real(i-1)/real(nX0-1)*(Xmax-Xmin)
      forall(i=1:nY0) X(2, :, i) = Ymin+real(i-1)/real(nY0-1)*(Ymax-Ymin)
      do i = 1,nY0
        info = fmodel%getProj(X(:,:,i), this%grid(:,:,i))
      end do
    endif
    call mlf_model_funbasis_init(this%model, this)
  End Function mlf_GridModelInit

  Function mlf_GridModelFunValue(this, W, x) result(Y)
    class(mlf_2dgrid_model), intent(in) :: this
    real(c_double), intent(in) :: W(:), x
    real(c_double) :: Y
    Y = this%top%funmodel%getValue(W, x)
  End Function mlf_GridModelFunValue

  integer Function mlf_GridModelGetProjection(this, Y, W, Aerror) result(info)
    class(mlf_2dgrid_model), intent(in), target :: this
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    real(c_double) :: vX, vY, ax, ay, dX, dY, diX, diY
    integer :: nXG, nYG, nIn, i, j, k
    info = -1
    if(size(Y,1) /= 2) RETURN
    nIn = size(Y,2)
    ASSOCIATE(grid => this%top%grid, XMin => this%top%XMin, XMax => this%top%XMax, &
        YMin => this%top%YMin, YMax => this%top%YMax)
      nXG = size(grid,2); nYG = size(grid,3)
      dX = (XMax-XMin)/real(nXG-1, kind = 8)
      dY = (YMax-YMin)/real(nYG-1, kind = 8)
      diX = 1d0/dX; diY = 1d0/dY
      do k = 1, nIn
        vX = Y(1,k); vY = Y(2,k)
        i = floor((vX-XMin)*diX)+1
        j = floor((vY-YMin)*diY)+1
        ax = ((vX-XMin)-(i-1)*dX)*diX
        ay = ((vY-YMin)-(j-1)*dY)*diY
        W(:,k) = (1d0-ay)*((1d0-ax)*grid(:,i,j)+ax*grid(:,i+1,j)) &
               + ay*((1d0-ax)*grid(:,i,j+1)+ax*grid(:,i+1,j+1))
      end do
    END ASSOCIATE
    ! TODO: add error estimation
    if(present(Aerror)) Aerror = 0
    info = 0
  End Function mlf_GridModelGetProjection
End Module mlf_gridmodel

