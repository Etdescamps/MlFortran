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
  Use mlf_fun_intf
  Use mlf_rsc_array
  Use mlf_models
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_obj_model) :: mlf_2dgrid
    class(mlf_reduce_model), pointer :: funmodel
    real(c_float), pointer :: grid(:,:,:)
    real(c_double), pointer :: XMin, XMax, YMin, YMax
  Contains
    procedure :: initF => mlf_GridModelInit
  End Type mlf_2dgrid

  Type, Public, Extends(mlf_b_reduce_model) :: mlf_2dgrid_model
    class(mlf_2dgrid), pointer :: top
  Contains
    procedure :: getProjSingleFloat => mlf_GridModelGetProjectionSingle
    procedure :: getProjSingle => mlf_GridModelGetProjectionSingleDouble
    procedure :: getProjMultFloat => mlf_GridModelGetProjection
    procedure :: getBounds => mlf_GridModelGetBounds
  End Type mlf_2dgrid_model

  Type, Public, Extends(mlf_2dgrid_model) :: mlf_2dgrid_fun_model
    class(mlf_bijection_fun), pointer :: fun
    real(c_double) :: XMin, XMax, YMin, YMax
  Contains
    procedure :: getProjSingleFloat => mlf_fun_GridModelGetProjectionSingle
    procedure :: getProjSingle => mlf_fun_GridModelGetProjectionSingleDouble
    procedure :: getProjMultFloat => mlf_fun_GridModelGetProjection
    procedure :: getBounds => mlf_fun_GridModelGetBounds
  End Type mlf_2dgrid_fun_model

Contains
  
  ! Model initilisator
  Subroutine mlf_model_funbasis_init(this, top, fun, XMin, XMax, YMin, YMax)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_2dgrid), intent(in), target :: top
    class(mlf_bijection_fun), intent(in), target, optional :: fun
    real(c_double), intent(in), optional :: XMin, XMax, YMin, YMax
    If(PRESENT(fun)) Then
      ALLOCATE(mlf_2dgrid_fun_model :: this)
    Else
      ALLOCATE(mlf_2dgrid_model :: this)
    Endif
    Select Type(this)
    Class is (mlf_2dgrid_fun_model)
      this%top => top
      this%fun => fun
      If(PRESENT(XMin)) this%XMin = XMin
      If(PRESENT(XMax)) this%XMax = XMax
      If(PRESENT(YMin)) this%YMin = YMin
      If(PRESENT(YMax)) this%YMax = YMax
    Class is (mlf_2dgrid_model)
      this%top => top
    End Select
  End Subroutine mlf_model_funbasis_init

  ! C wrapper for init function
  Type(c_ptr) Function c_2dgridmodel_init(cmodel, XMin, XMax, YMin, YMax, nW, nX0, nY0) &
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
    If(.NOT. ASSOCIATED(model)) RETURN
    Select Type(model)
      Class is (mlf_reduce_model)
        ALLOCATE(x)
        info = x%initF(model, XMin, XMax, YMin, YMax, nW, nX0, nY0)
        If(info < 0) RETURN
        obj => x
        c_2dgridmodel_init = c_allocate(obj)
    End Select
  End Function c_2dgridmodel_init

  Integer Function mlf_GridModelInit(this, fmodel, XMin, XMax, YMin, YMax, nW, nX0, nY0, fun, nFastX, data_handler) result(info)
    class(mlf_2dgrid), intent(inout) :: this
    class(mlf_reduce_model), target :: fmodel
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: XMin, XMax, YMin, YMax
    integer(c_int), intent(in), optional :: nX0, nY0, nW, nFastX
    real(c_double), allocatable :: X(:,:,:), G(:,:)
    class(mlf_bijection_fun), intent(in), target, optional :: fun
    type(mlf_rsc_numFields) :: numFields
    real(c_double) :: XYMin(2), XYMax(2), XTemp(2)
    integer :: i, j
    integer(c_int64_t) :: ndGrid(3)
    numFields = mlf_rsc_numFields(0,4,1)
    ndGrid = -1
    If(PRESENT(nW)) ndGrid(1) = nW
    If(PRESENT(nX0)) ndGrid(2) = nX0
    If(PRESENT(nY0)) ndGrid(3) = nY0
    info = mlf_arr_init(this, numFields, data_handler)
    If(info < 0) RETURN
    CALL this%addRPar(numFields, this%XMin, "XMin")
    CALL this%addRPar(numFields, this%XMax, "XMax")
    CALL this%addRPar(numFields, this%YMin, "YMin")
    CALL this%addRPar(numFields, this%YMax, "YMax")
    info = this%add_FMatrix3D(numFields, ndGrid, this%grid, C_CHAR_"grid", data_handler = data_handler)
    If(info < 0) RETURN
    this%funmodel => fmodel
    If(.NOT. PRESENT(data_handler)) Then
      ALLOCATE(X(2,nX0,nY0), G(nW, nX0))
      If(PRESENT(fun)) Then
        CALL fun%proj([XMin, YMin], XYMin)
        CALL fun%proj([XMax, YMax], XYMax)
        this%XMin = XYMin(1); this%YMin = XYMin(2)
        this%XMax = XYMax(1); this%YMax = XYMax(2)
      Else
        this%XMin = XMin; this%XMax = XMax; this%YMin = YMin; this%YMax = YMax
      Endif
      FORALL(i=1:nX0) X(1, i, :) = this%Xmin+REAL(i-1)/REAL(nX0-1)*(this%Xmax-this%Xmin)
      FORALL(i=1:nY0) X(2, :, i) = this%Ymin+REAL(i-1)/REAL(nY0-1)*(this%Ymax-this%Ymin)
      If(PRESENT(fun)) Then
        Do i=1,nX0
          Do j=1,nY0
            CALL fun%invert(X(:,i,j), XTemp)
            X(:,i,j)=XTemp
          End Do
        End Do
      Endif
      Do i = 1,nY0
        If(PRESENT(nFastX)) Then
          Select Type(fmodel)
          Class is (mlf_fast_approx_linear)
            info = fmodel%getFastProj(X(:,:,i), nFastX, G)
          Class Default
            info = fmodel%getProj(X(:,:,i), G)
          End Select
        Else
          info = fmodel%getProj(X(:,:,i), G)
        Endif
        this%grid(:,:,i) = REAL(G,4)
      End Do
    Endif
    CALL mlf_model_funbasis_init(this%model, this, fun, XMin, XMax, YMin, YMax)
  End Function mlf_GridModelInit
  
  Subroutine mlf_GridModelGetBounds(this, XMin, XMax)
    class(mlf_2dgrid_model), intent(in), target :: this
    real(c_double), intent(out) :: XMin(:), XMax(:)
    XMin(1) = this%top%XMin
    XMin(2) = this%top%YMin
    XMax(1) = this%top%XMax
    XMax(2) = this%top%YMax
  End Subroutine mlf_GridModelGetBounds

  Subroutine mlf_fun_GridModelGetBounds(this, XMin, XMax)
    class(mlf_2dgrid_fun_model), intent(in), target :: this
    real(c_double), intent(out) :: XMin(:), XMax(:)
    XMin(1) = this%XMin
    XMin(2) = this%YMin
    XMax(1) = this%XMax
    XMax(2) = this%YMax
  End Subroutine mlf_fun_GridModelGetBounds


  Integer Function mlf_GridModelGetProjectionSingleDouble(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_model), intent(in), target :: this
    real(c_double), intent(in) :: Y(:)
    real(c_double), intent(out) :: W(:)
    real(c_double), optional, intent(out) :: Aerror(:)
    real(c_float) :: W0(SIZE(W))
    info = this%getProj(REAL(Y, c_float), W0)
    W = REAL(W0, c_double)
    If(PRESENT(Aerror)) Aerror = 0
  End Function mlf_GridModelGetProjectionSingleDouble

  Integer Function mlf_fun_GridModelGetProjectionSingleDouble(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_fun_model), intent(in), target :: this
    real(c_double), intent(in) :: Y(:)
    real(c_double), intent(out) :: W(:)
    real(c_double), optional, intent(out) :: Aerror(:)
    real(c_double) :: Y0(SIZE(Y))
    CALL this%fun%proj(Y, Y0)
    info = mlf_GridModelGetProjectionSingleDouble(this, Y0, W, Aerror)
  End Function mlf_fun_GridModelGetProjectionSingleDouble

  Integer Function mlf_GridModelGetProjectionSingle(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:)
    real(c_float), intent(out) :: W(:)
    real(c_float), optional, intent(out) :: Aerror(:)
    real(c_double) :: vX, vY, ax, ay, dX, dY, diX, diY
    integer :: nXG, nYG, i, j
    info = -1
    If(SIZE(Y,1) /= 2) RETURN
    ASSOCIATE(grid => this%top%grid, XMin => this%top%XMin, XMax => this%top%XMax, &
        YMin => this%top%YMin, YMax => this%top%YMax)
      nXG = SIZE(grid,2); nYG = SIZE(grid,3)
      dX = (XMax-XMin)/REAL(nXG-1, KIND = 4)
      dY = (YMax-YMin)/REAL(nYG-1, KIND = 4)
      diX = 1d0/dX; diY = 1d0/dY
      vX = Y(1); vY = Y(2)
      i = FLOOR((vX-XMin)*diX)+1
      j = FLOOR((vY-YMin)*diY)+1
      ax = ((vX-XMin)-(i-1)*dX)*diX
      ay = ((vY-YMin)-(j-1)*dY)*diY
      W(:) = (1d0-ay)*((1d0-ax)*grid(:,i,j)+ax*grid(:,i+1,j)) &
             + ay*((1d0-ax)*grid(:,i,j+1)+ax*grid(:,i+1,j+1))
    END ASSOCIATE
    ! TODO: add error estimation
    If(PRESENT(Aerror)) Aerror = 0
    info = 0
  End Function mlf_GridModelGetProjectionSingle

  Integer Function mlf_fun_GridModelGetProjectionSingle(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_fun_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:)
    real(c_float), intent(out) :: W(:)
    real(c_float), optional, intent(out) :: Aerror(:)
    real(c_double) :: Y0(SIZE(Y))
    CALL this%fun%proj(REAL(Y, KIND=8), Y0)
    info = mlf_GridModelGetProjectionSingle(this, REAL(Y0, KIND=4), W, Aerror)
  End Function mlf_fun_GridModelGetProjectionSingle

  Integer Function mlf_GridModelGetProjection(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:,:)
    real(c_float), intent(out) :: W(:,:)
    real(c_float), optional, intent(out) :: Aerror(:,:)
    real(c_double) :: vX, vY, ax, ay, dX, dY, diX, diY
    integer :: nXG, nYG, nIn, i, j, k
    info = -1
    If(SIZE(Y,1) /= 2) RETURN
    nIn = SIZE(Y,2)
    ASSOCIATE(grid => this%top%grid, XMin => this%top%XMin, XMax => this%top%XMax, &
        YMin => this%top%YMin, YMax => this%top%YMax)
      nXG = SIZE(grid,2); nYG = SIZE(grid,3)
      dX = (XMax-XMin)/REAL(nXG-1, KIND = 8)
      dY = (YMax-YMin)/REAL(nYG-1, KIND = 8)
      diX = 1d0/dX; diY = 1d0/dY
      Do k = 1, nIn
        vX = Y(1,k); vY = Y(2,k)
        i = FLOOR((vX-XMin)*diX)+1
        j = FLOOR((vY-YMin)*diY)+1
        ax = ((vX-XMin)-(i-1)*dX)*diX
        ay = ((vY-YMin)-(j-1)*dY)*diY
        W(:,k) = (1d0-ay)*((1d0-ax)*grid(:,i,j)+ax*grid(:,i+1,j)) &
               + ay*((1d0-ax)*grid(:,i,j+1)+ax*grid(:,i+1,j+1))
      End Do
    END ASSOCIATE
    ! TODO: add error estimation
    If(PRESENT(Aerror)) Aerror = 0
    info = 0
  End Function mlf_GridModelGetProjection

  Integer Function mlf_fun_GridModelGetProjection(this, Y, W, Aerror) Result(info)
    class(mlf_2dgrid_fun_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:,:)
    real(c_float), intent(out) :: W(:,:)
    real(c_float), optional, intent(out) :: Aerror(:,:)
    real(c_float), allocatable :: X(:,:)
    real(c_double) :: Y0(SIZE(Y,1))
    integer :: i
    ALLOCATE(X, MOLD=Y)
    Do i=1,SIZE(Y,2)
      CALL this%fun%proj(REAL(Y(:,i), KIND=8), Y0)
      X(:,i) = REAL(Y0, KIND=4)
    End Do
    info = mlf_GridModelGetProjection(this, X, W, Aerror)
  End Function mlf_fun_GridModelGetProjection
End Module mlf_gridmodel

