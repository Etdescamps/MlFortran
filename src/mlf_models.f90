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

  Type, public :: mlf_model
  End Type mlf_model
  
  Type, public, abstract, extends(mlf_arr_obj) :: mlf_obj_model
    class(mlf_model), allocatable :: model
  End Type mlf_obj_model

  Type, public, abstract, extends(mlf_model) :: mlf_class_model
  Contains
    procedure (mlf_model_getClass), deferred :: getClass
  End Type mlf_class_model
   
  Type, public, abstract, extends(mlf_class_model) :: mlf_class_proba_model
  Contains
    procedure (mlf_model_getProba), deferred :: getProba
    procedure (mlf_model_getNumClasses), deferred :: getNumClasses
    procedure :: getClass => mlf_model_getClass_proba
  End Type mlf_class_proba_model

  Type, public, abstract, extends(mlf_model) :: mlf_reduce_model
  Contains
    procedure (mlf_model_getProjSingle), deferred :: getProjSingle
    procedure :: getProjMult => mlf_model_getProjMult
    generic :: getProj => getProjSingle, getProjMult
  End Type mlf_reduce_model

  Type, public, abstract, extends(mlf_reduce_model) :: mlf_approx_model
  Contains
    procedure (mlf_model_getValue), deferred :: getValue
  End Type mlf_approx_model

  Type, public, abstract, extends(mlf_approx_model) :: mlf_approx_linear
  Contains
    procedure :: getValue => mlf_approx_linear_getValue
    procedure (mlf_model_getValueBasis), deferred :: getValueBasis
  End Type mlf_approx_linear


  Abstract Interface
    integer Function mlf_model_getClass(this, X, Cl)
      use iso_c_binding
      import :: mlf_class_model
      class(mlf_class_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      integer(c_int), intent(out) :: Cl(:)
    End Function mlf_model_getClass

    integer Function mlf_model_getNumClasses(this)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
    End Function mlf_model_getNumClasses

    integer Function mlf_model_getProba(this, X, Proba)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      real(c_double), intent(out) :: Proba(:,:)
    End Function mlf_model_getProba

    integer Function mlf_model_getProjSingle(this, Y, W, Aerror)
      use iso_c_binding
      import :: mlf_reduce_model
      class(mlf_reduce_model), intent(in), target :: this
      real(c_double), intent(in) :: Y(:)
      real(c_double), intent(out) :: W(:)
      real(c_double), optional, intent(out) :: Aerror(:)
    End Function mlf_model_getProjSingle

    Function mlf_model_getValue(this, W, x) result(Y)
      use iso_c_binding
      import :: mlf_approx_model
      class(mlf_approx_model), intent(in) :: this
      real(c_double), intent(in) :: W(:), x
      real(c_double) :: Y
    End Function mlf_model_getValue

    Integer Function mlf_model_getValueBasis(this, x, Y) result(info)
      use iso_c_binding
      import :: mlf_approx_linear
      class(mlf_approx_linear), intent(in) :: this
      real(c_double), intent(in) :: x
      real(c_double), intent(out) :: Y(:)
    End Function mlf_model_getValueBasis
  End Interface
Contains
  integer Function mlf_model_getProjMult(this, Y, W, Aerror) result(info)
    class(mlf_reduce_model), intent(in), target :: this
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    integer :: i
    Do i=1,size(Y,2)
      info = this%getProj(Y(:,i), W(:,i), Aerror(:,i))
    End Do
  End Function mlf_model_getProjMult


  real(c_double) Function mlf_approx_linear_getValue(this, W, x) result(Y)
    class(mlf_approx_linear), intent(in) :: this
    real(c_double), intent(in) :: W(:), x
    real(c_double) :: Z(size(W))
    integer :: info
    info = this%getValueBasis(x, Z)
    Y = dot_product(W,Z)
  End Function mlf_approx_linear_getValue

  integer Function mlf_model_getClass_proba(this, X, Cl) result(info)
    class(mlf_class_proba_model), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    integer(c_int), intent(out) :: Cl(:)
    real(c_double), allocatable :: Proba(:,:)
    integer :: nX, nC
    nX = size(X,2); nC = this%getNumClasses()
    allocate(Proba(nX, nC))
    info = this%getProba(X, Proba)
    cl = maxloc(Proba, dim=2)
  End Function mlf_model_getClass_proba

  Function mlf_getModel(cptr) result(model)
    type(c_ptr), value :: cptr
    class(*), pointer :: obj
    class(mlf_model), pointer :: model
    model => NULL()
    obj => mlf_getrawfromc(cptr)
    if(.NOT. associated(obj)) RETURN
    select type(obj)
      class is (mlf_obj_model)
        if(allocated(obj%model)) model => obj%model
      class is (mlf_model)
        model => obj
    end select
  End Function mlf_getModel

  integer(c_int) Function mlf_getClass(model, X, Cl) result(info)
    class(mlf_model), intent(inout) :: model
    integer(c_int), intent(out) :: Cl(:)
    real(c_double), intent(in) :: X(:,:)
    info = -1
    select type(model)
      class is (mlf_class_model)
        info = model%getClass(X, Cl)
    end select
  End Function mlf_getClass

  integer(c_int) Function c_getClass(cptr, cX, cCl, nX, nIn) result(info) bind(C, name="mlf_getClass")
    type(c_ptr), value :: cptr, cX, cCl
    integer(c_int), value :: nX, nIn
    real(c_double), pointer :: X(:,:)
    class(mlf_model), pointer :: model
    integer(c_int), pointer :: Cl(:)
    info = -1
    model => mlf_getModel(cptr)
    if(.NOT. associated(model)) RETURN
    call C_F_POINTER(cX, X, [nX, nIn])
    call C_F_POINTER(cCl, Cl, [nIn])
    info = mlf_getClass(model, X, Cl)
  End Function c_getClass

  integer(c_int) Function mlf_getProba(model, X, Cl) result(info)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(out) :: Cl(:,:)
    real(c_double), intent(in) :: X(:,:)
    info = -1
    select type(model)
      class is (mlf_class_proba_model)
        info = model%getProba(X, Cl)
    end select
  End Function mlf_getProba

  integer(c_int) Function c_getProba(cptr, cX, cCl, nX, nIn, nCl) result(info) bind(C, name="mlf_getProba")
    type(c_ptr), value :: cptr, cX, cCl
    integer(c_int), value :: nX, nIn, nCl
    class(mlf_model), pointer :: model
    real(c_double), pointer :: X(:,:), Cl(:,:)
    info = -1
    model => mlf_getModel(cptr)
    if(.NOT. associated(model)) RETURN
    call C_F_POINTER(cX, X, [nX, nIn])
    call C_F_POINTER(cCl, Cl, [nCl, nIn])
    info = mlf_getProba(model, X, Cl)
  End Function c_getProba

  integer(c_int) Function mlf_getNumClasses(model) result(info)
    class(mlf_model), intent(inout) :: model
    info = -1
    select type(model)
      class is (mlf_class_proba_model)
        info = model%getNumClasses()
    end select
  End Function mlf_getNumClasses

  integer(c_int) Function c_getNumClasses(cptr) result(info) bind(C, name="mlf_getNumClasses")
    type(c_ptr), value :: cptr
    class(mlf_model), pointer :: model
    info = -1
    model => mlf_getModel(cptr)
    if(.NOT. associated(model)) RETURN
    info = mlf_getNumClasses(model)
  End Function c_getNumClasses

  integer(c_int) Function mlf_getProj(model, Y, W, Aerror) result(info)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    info = -1
    select type(model)
      class is (mlf_reduce_model)
        info = model%getProj(Y,W,Aerror)
    end select
  End Function mlf_getProj

  integer(c_int) Function c_getProj(cptr, cY, cW, nIn, nDimIn, nDimOut) result(info) bind(C, name="mlf_getProj")
    type(c_ptr), value :: cptr, cY, cW
    integer(c_int), value :: nIn
    class(mlf_model), pointer :: model
    real(c_double), pointer :: Y(:,:), W(:,:)
    integer(c_int), value :: nDimIn, nDimOut
    info = -1
    model => mlf_getModel(cptr)
    if(.NOT. associated(model)) RETURN
    call C_F_POINTER(cY, Y, [nDimIn, nIn])
    call C_F_POINTER(cW, W, [nDimOut, nIn])
    info = mlf_getProj(model, Y, W)
  End Function c_getProj

  real(c_double) Function mlf_getValue(model, W, t) result(Y)
    class(mlf_model), intent(inout) :: model
    real(c_double), intent(in) :: W(:), t
    Y = IEEE_VALUE(Y, IEEE_QUIET_NAN)
    select type(model)
      class is (mlf_approx_model)
        Y = model%getValue(W,t)
    end select
  End Function mlf_getValue

  real(c_double) Function c_getValue(cptr, cW, t, nDimOut) result(Y) bind(C, name="mlf_getValue")
    type(c_ptr), value :: cptr, cW
    real(c_double), value :: t
    class(mlf_model), pointer :: model
    real(c_double), pointer :: W(:)
    integer(c_int), value :: nDimOut
    Y = IEEE_VALUE(Y, IEEE_QUIET_NAN)
    model => mlf_getModel(cptr)
    if(.NOT. associated(model)) RETURN
    call C_F_POINTER(cW, W, [nDimOut])
    Y = mlf_getValue(model, W, t)
  End Function c_getValue


End Module mlf_models

