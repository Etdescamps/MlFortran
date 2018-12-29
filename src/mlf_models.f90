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

Module mlf_models
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_getClass, mlf_getProba, mlf_getNumClasses, mlf_getValue, mlf_getProj
  Public :: mlf_getModel

  Type, Public :: mlf_model
  End Type mlf_model
  
  Type, Public, Abstract, Extends(mlf_arr_obj) :: mlf_obj_model
    class(mlf_model), allocatable :: model
  End Type mlf_obj_model

  Type, Public, abstract, extends(mlf_model) :: mlf_class_model
  Contains
    procedure (mlf_model_getClass), deferred :: getClass
  End Type mlf_class_model
   
  Type, Public, Abstract, Extends(mlf_class_model) :: mlf_class_proba_model
  Contains
    procedure (mlf_model_getProba), deferred :: getProba
    procedure (mlf_model_getNumClasses), deferred :: getNumClasses
    procedure :: getClass => mlf_model_getClass_proba
  End Type mlf_class_proba_model

  Type, Public, Abstract, Extends(mlf_model) :: mlf_reduce_model
  Contains
    procedure (mlf_model_getProjSingle), deferred :: getProjSingle
    procedure :: getProjSingleFloat => mlf_model_getProjSingleFloat
    procedure :: getProjMult => mlf_model_getProjMult
    procedure :: getProjMultFloat => mlf_model_getProjMultFloat
    generic :: getProj => getProjSingle, getProjMult, getProjSingleFloat, &
      getProjMultFloat
  End Type mlf_reduce_model

  Type, Public, Abstract, Extends(mlf_reduce_model) :: mlf_b_reduce_model
  Contains
    procedure (mlf_model_reduce_getBounds), deferred :: getBounds
  End Type mlf_b_reduce_model

  Type, Public, Abstract, Extends(mlf_reduce_model) :: mlf_approx_model
  Contains
    procedure (mlf_model_getValue), deferred :: getValue
    procedure (mlf_model_getValueBounds), deferred :: getValueBounds
  End Type mlf_approx_model

  Type, Public, Abstract, Extends(mlf_approx_model) :: mlf_approx_linear
  Contains
    procedure :: getValue => mlf_approx_linear_getValue
    procedure (mlf_model_getValueBasis), deferred :: getValueBasis
  End Type mlf_approx_linear

  Type, Public, Abstract, Extends(mlf_approx_linear) :: mlf_fast_approx_linear
  Contains
    procedure (mlf_model_getFastProj), deferred :: getFastProj
  End Type mlf_fast_approx_linear

  Abstract Interface
    Integer Function mlf_model_getClass(this, X, Cl)
      use iso_c_binding
      import :: mlf_class_model
      class(mlf_class_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      integer(c_int), intent(out) :: Cl(:)
    End Function mlf_model_getClass

    Integer Function mlf_model_getNumClasses(this)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
    End Function mlf_model_getNumClasses

    Integer Function mlf_model_getProba(this, X, Proba)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      real(c_double), intent(out) :: Proba(:,:)
    End Function mlf_model_getProba

    Subroutine mlf_model_reduce_getBounds(this, XMin, XMax)
      use iso_c_binding
      import :: mlf_b_reduce_model
      class(mlf_b_reduce_model), intent(in), target :: this
      real(c_double), intent(out) :: XMin(:), XMax(:)
    End Subroutine mlf_model_reduce_getBounds

    Integer Function mlf_model_getProjSingle(this, Y, W, Aerror)
      use iso_c_binding
      import :: mlf_reduce_model
      class(mlf_reduce_model), intent(in), target :: this
      real(c_double), intent(in) :: Y(:)
      real(c_double), intent(out) :: W(:)
      real(c_double), optional, intent(out) :: Aerror(:)
    End Function mlf_model_getProjSingle

    Function mlf_model_getValue(this, W, x) Result(Y)
      use iso_c_binding
      import :: mlf_approx_model
      class(mlf_approx_model), intent(in) :: this
      real(c_double), intent(in) :: W(:), x
      real(c_double) :: Y
    End Function mlf_model_getValue

    Integer Function mlf_model_getValueBasis(this, x, Y) Result(info)
      use iso_c_binding
      import :: mlf_approx_linear
      class(mlf_approx_linear), intent(in) :: this
      real(c_double), intent(in) :: x
      real(c_double), intent(out) :: Y(:)
    End Function mlf_model_getValueBasis

    Integer Function mlf_model_getFastProj(this, P, nX, W) Result(info)
      use iso_c_binding
      import :: mlf_fast_approx_linear
      class(mlf_fast_approx_linear), intent(in) :: this
      real(c_double), intent(in) :: P(:,:)
      integer, intent(in) :: nX
      real(c_double), intent(out) :: W(:,:)
    End Function mlf_model_getFastProj

    Subroutine mlf_model_getValueBounds(this, xMin, xMax)
      use iso_c_binding
      import :: mlf_approx_model
      class(mlf_approx_model), intent(in) :: this
      real(c_double), intent(out), optional :: xMin, xMax
    End Subroutine mlf_model_getValueBounds
  End Interface
Contains
  Integer Function mlf_model_getProjSingleFloat(this, Y, W, Aerror) Result(info)
    class(mlf_reduce_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:)
    real(c_float), intent(out) :: W(:)
    real(c_float), optional, intent(out) :: Aerror(:)
    real(c_double) :: Y0(SIZE(Y))
    real(c_double) :: W0(SIZE(W))
    integer :: l
    Y0 = REAL(Y, c_double)
    If(PRESENT(Aerror)) Then
      l = SIZE(Aerror)
      BLOCK
        real(c_double) :: Aerror0(l)
        info = this%getProj(Y0, W0, Aerror0)
        Aerror = REAL(Aerror0, c_float)
      END BLOCK
    Else
      info = this%getProj(Y0, W0)
    Endif
    W = REAL(W0, c_float)
  End Function mlf_model_getProjSingleFloat

  Integer Function mlf_model_getProjMult(this, Y, W, Aerror) Result(info)
    class(mlf_reduce_model), intent(in), target :: this
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    integer :: i
    info = 0
    Do i=1,size(Y,2)
      info = this%getProj(Y(:,i), W(:,i), Aerror(:,i))
    End Do
  End Function mlf_model_getProjMult

  Integer Function mlf_model_getProjMultFloat(this, Y, W, Aerror) Result(info)
    class(mlf_reduce_model), intent(in), target :: this
    real(c_float), intent(in) :: Y(:,:)
    real(c_float), intent(out) :: W(:,:)
    real(c_float), optional, intent(out) :: Aerror(:,:)
    real(c_double) :: Y0(SIZE(Y,1))
    real(c_double) :: W0(SIZE(W,1))
    integer :: i, l
    info = 0
    If(PRESENT(Aerror)) Then
      l = SIZE(Aerror, 1)
      BLOCK
        real(c_double) :: Aerror0(l)
        Do i=1,SIZE(Y,2)
          Y0 = REAL(Y(:,i), c_double)
          info = this%getProj(Y0, W0, Aerror0)
          W(:,i) = REAL(W0, c_float)
          Aerror(:,i) = REAL(Aerror0, c_float)
        End Do
      END BLOCK
    Else
      Do i=1,SIZE(Y,2)
        Y0 = REAL(Y(:,i), c_double)
        info = this%getProj(Y0, W0)
        W(:,i) = REAL(W0, c_float)
      End Do
    Endif
  End Function mlf_model_getProjMultFloat

  Real(c_double) Function mlf_approx_linear_getValue(this, W, x) Result(Y)
    class(mlf_approx_linear), intent(in) :: this
    real(c_double), intent(in) :: W(:), x
    real(c_double) :: Z(SIZE(W))
    integer :: info
    info = this%getValueBasis(x, Z)
    If(info == 0) Then
      Y = DOT_PRODUCT(W,Z)
    Else
      Y = IEEE_VALUE(0d0, ieee_quiet_nan)
    Endif
  End Function mlf_approx_linear_getValue

  Integer Function mlf_model_getClass_proba(this, X, Cl) Result(info)
    class(mlf_class_proba_model), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    integer(c_int), intent(out) :: Cl(:)
    real(c_double), allocatable :: Proba(:,:)
    integer :: nX, nC
    nX = SIZE(X,2); nC = this%getNumClasses()
    ALLOCATE(Proba(nX, nC))
    info = this%getProba(X, Proba)
    cl = MAXLOC(Proba, DIM=2)
  End Function mlf_model_getClass_proba

  Function mlf_getModel(cptr) result(model)
    type(c_ptr), value :: cptr
    class(*), pointer :: obj
    class(mlf_model), pointer :: model
    model => NULL()
    obj => mlf_getrawfromc(cptr)
    If(.NOT. associated(obj)) RETURN
    Select Type(obj)
    Class is (mlf_obj_model)
      If(allocated(obj%model)) model => obj%model
    Class is (mlf_model)
      model => obj
    End Select
  End Function mlf_getModel

  Integer(c_int) Function mlf_getClass(model, X, Cl) Result(info)
    class(mlf_model), intent(inout) :: model
    integer(c_int), intent(out) :: Cl(:)
    real(c_double), intent(in) :: X(:,:)
    info = -1
    Select Type(model)
    Class is (mlf_class_model)
      info = model%getClass(X, Cl)
    End Select
  End Function mlf_getClass

  Integer(c_int) Function c_getClass(cptr, cX, cCl, nX, nIn) Result(info) Bind(C, name="mlf_getClass")
    type(c_ptr), value :: cptr, cX, cCl
    integer(c_int), value :: nX, nIn
    real(c_double), pointer :: X(:,:)
    class(mlf_model), pointer :: model
    integer(c_int), pointer :: Cl(:)
    info = -1
    model => mlf_getModel(cptr)
    If(.NOT. ASSOCIATED(model)) RETURN
    CALL C_F_POINTER(cX, X, [nX, nIn])
    CALL C_F_POINTER(cCl, Cl, [nIn])
    info = mlf_getClass(model, X, Cl)
  End Function c_getClass

  Integer(c_int) Function mlf_getProba(model, X, Cl) Result(info)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(out) :: Cl(:,:)
    real(c_double), intent(in) :: X(:,:)
    info = -1
    Select Type(model)
    Class is (mlf_class_proba_model)
      info = model%getProba(X, Cl)
    End Select
  End Function mlf_getProba

  Integer(c_int) Function c_getProba(cptr, cX, cCl, nX, nIn, nCl) Result(info) Bind(C, name="mlf_getProba")
    type(c_ptr), value :: cptr, cX, cCl
    integer(c_int), value :: nX, nIn, nCl
    class(mlf_model), pointer :: model
    real(c_double), pointer :: X(:,:), Cl(:,:)
    info = -1
    model => mlf_getModel(cptr)
    If(.NOT. ASSOCIATED(model)) RETURN
    CALL C_F_POINTER(cX, X, [nX, nIn])
    CALL C_F_POINTER(cCl, Cl, [nCl, nIn])
    info = mlf_getProba(model, X, Cl)
  End Function c_getProba

  Integer(c_int) Function mlf_getNumClasses(model) Result(info)
    class(mlf_model), intent(inout) :: model
    info = -1
    Select Type(model)
    Class is (mlf_class_proba_model)
      info = model%getNumClasses()
    End Select
  End Function mlf_getNumClasses

  Integer(c_int) Function c_getNumClasses(cptr) Result(info) Bind(C, name="mlf_getNumClasses")
    type(c_ptr), value :: cptr
    class(mlf_model), pointer :: model
    info = -1
    model => mlf_getModel(cptr)
    If(.NOT. associated(model)) RETURN
    info = mlf_getNumClasses(model)
  End Function c_getNumClasses

  Integer(c_int) Function mlf_getProj(model, Y, W, Aerror) Result(info)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    info = -1
    Select Type(model)
    Class is (mlf_reduce_model)
      info = model%getProj(Y,W,Aerror)
    End Select
  End Function mlf_getProj

  Integer(c_int) Function mlf_getProj_f(model, Y, W, Aerror) Result(info)
    class(mlf_model), intent(inout) :: model
    real(c_float), intent(in) :: Y(:,:)
    real(c_float), intent(out) :: W(:,:)
    real(c_float), optional, intent(out) :: Aerror(:,:)
    info = -1
    Select Type(model)
    Class is (mlf_reduce_model)
      info = model%getProj(Y,W,Aerror)
    End Select
  End Function mlf_getProj_f

  Integer(c_int) Function c_getProj(cptr, cY, cW, nIn, nDimIn, nDimOut) Result(info) Bind(C, name="mlf_getProj")
    type(c_ptr), value :: cptr, cY, cW
    integer(c_int), value :: nIn
    class(mlf_model), pointer :: model
    real(c_double), pointer :: Y(:,:), W(:,:)
    integer(c_int), value :: nDimIn, nDimOut
    info = -1
    model => mlf_getModel(cptr)
    If(.NOT. ASSOCIATED(model)) RETURN
    CALL C_F_POINTER(cY, Y, [nDimIn, nIn])
    CALL C_F_POINTER(cW, W, [nDimOut, nIn])
    info = mlf_getProj(model, Y, W)
  End Function c_getProj

  Integer(c_int) Function c_getProj_f(cptr, cY, cW, nIn, nDimIn, nDimOut) Result(info) Bind(C, name="mlf_getProj_f")
    type(c_ptr), value :: cptr, cY, cW
    integer(c_int), value :: nIn
    class(mlf_model), pointer :: model
    real(c_float), pointer :: Y(:,:), W(:,:)
    integer(c_int), value :: nDimIn, nDimOut
    info = -1
    model => mlf_getModel(cptr)
    If(.NOT. ASSOCIATED(model)) RETURN
    CALL C_F_POINTER(cY, Y, [nDimIn, nIn])
    CALL C_F_POINTER(cW, W, [nDimOut, nIn])
    info = mlf_getProj_f(model, Y, W)
  End Function c_getProj_f


  Real(c_double) Function mlf_getValue(model, W, t) Result(Y)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(in) :: W(:), t
    Y = IEEE_VALUE(Y, IEEE_QUIET_NAN)
    Select type(model)
    Class is (mlf_approx_model)
      Y = model%getValue(W,t)
    End Select
  End Function mlf_getValue

  Real(c_double) Function c_getValue(cptr, cW, t, nDimOut) Result(Y) Bind(C, name="mlf_getValue")
    type(c_ptr), value :: cptr, cW
    real(c_double), value :: t
    class(mlf_model), pointer :: model
    real(c_double), pointer :: W(:)
    integer(c_int), value :: nDimOut
    Y = IEEE_VALUE(Y, IEEE_QUIET_NAN)
    model => mlf_getModel(cptr)
    If(.NOT. ASSOCIATED(model)) RETURN
    CALL C_F_POINTER(cW, W, [nDimOut])
    Y = mlf_getValue(model, W, t)
  End Function c_getValue
End Module mlf_models

